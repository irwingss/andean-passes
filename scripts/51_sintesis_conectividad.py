"""51_sintesis_conectividad.py - Synthesis of the present climatic connectivity (PRD C8 and C9), without inference (D21).
Reads data/interim/censo/passes_andes.csv, the pass values written by 49_passes_corriente.py and the tile tables of
50_conectividad_teselas.py for one run, and writes:
  data/interim/conectividad/<run>/normalizada.npy   normalized current, array [configuration, reading, variable, pass]
  outputs/database/APD_passes.csv.gz                 the census with, for tropical passes, the median and range over the 19
                                                   variables of the normalized current in each reading and configuration
  outputs/tables/Table1_top_passes.csv, Table2_control_passes.csv
  outputs/supp/tables/S*.csv                       census, validation, envelopes, profiles, variables and sensitivity
  data/interim/conectividad/<run>/grupos/<group>_<reading>_<tag>.tif   median normalized current of each cell over the
                                                   temperature (temp, BIO1 to BIO11) or precipitation (prec, BIO12 to
                                                   BIO19) variables, per configuration (tag as in the tile maps)
python -P scripts/51_sintesis_conectividad.py principal
Configurations: q 0.1 with 200 km tiles (main), q 0.05 and 0.2 with 200 km, and 100 and 400 km with q 0.1. A statement
about a pass holds when it holds in all five (C8). Passes in the same 1 km cell share their values; rankings and
correlations keep one pass per cell, the deepest. The normalized current is left undefined where the current under
uniform resistance (corriente / normalizada) is below 1e-6 A, where the ratio of two near-zero sums is unstable.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer

RUN = sys.argv[1]
D = f"data/interim/conectividad/{RUN}"
CFG = [(0.1, 200), (0.05, 200), (0.2, 200), (0.1, 100), (0.1, 400)]      # first one is the main run
CFG_LAB = ["main", "q 0.05", "q 0.2", "100 km", "400 km"]
LECT = ["baja", "alta"]                                                  # transverse (lowland), longitudinal (montane)
LECT_LAB = {"baja": "transverse", "alta": "longitudinal"}
VARS = [f"bio{v:02d}" for v in range(1, 20)]
LN, LS = 11.0, -23.5                                                     # reported tropics (D1)
TAB, SUP, DB = "outputs/tables", "outputs/supp/tables", "outputs/database"
for d in (TAB, SUP, DB):
    os.makedirs(d, exist_ok=True)

passes = pd.read_csv("data/interim/censo/passes_andes.csv", dtype={"hybas_term_a": str, "hybas_term_b": str})
N = len(passes)
trop = ((passes.lat <= LN) & (passes.lat >= LS)).to_numpy()

# normalized current per configuration, reading, variable and pass (cached: reading the 63 million rows takes minutes)
NPY = f"{D}/normalizada.npy"
U_MIN = 1e-6      # below this current under uniform resistance (of the 1 A per tile) the ratio is numerically unstable
if not os.path.exists(NPY) or os.path.getmtime(NPY) < os.path.getmtime(f"{D}/passes_corriente.csv.gz"):
    A = np.full((len(CFG), 2, 19, N), np.nan, np.float32)
    ci = {c: i for i, c in enumerate(CFG)}
    for ch in pd.read_csv(f"{D}/passes_corriente.csv.gz", chunksize=5_000_000,
                          dtype={"pass_id": np.int32, "var": str, "lectura": str, "q": np.float64,
                                 "lado_km": np.int32, "corriente": np.float32, "normalizada": np.float32}):
        k = np.array([ci[(q, s)] for q, s in zip(ch.q.round(3), ch.lado_km)])
        A[k, (ch.lectura == "alta").to_numpy().astype(int), ch["var"].str[3:].astype(int).to_numpy() - 1,
          ch.pass_id.to_numpy()] = np.where(ch.corriente / ch.normalizada >= U_MIN, ch.normalizada, np.nan)
    np.save(NPY, A)
A = np.load(NPY)
with np.errstate(all="ignore"):
    MED = np.nanmedian(A, axis=2)                      # [cfg, reading, pass], median over the 19 variables
    MIN, MAX = np.nanmin(A, axis=2), np.nanmax(A, axis=2)
NV = np.isfinite(A).sum(axis=2)
MED[NV < 19] = np.nan                                  # a pass needs the 19 variables
ROB_HI = np.all(MED > 1, axis=0)                       # concentrates current in all five configurations
ROB_LO = np.all(MED < 1, axis=0)                       # diverts current in all five

# 1 km cell of every pass, and one pass per cell (the deepest) for rankings and correlations
with rasterio.open(f"{D}/insumos/elev.tif") as r:
    to = Transformer.from_crs("EPSG:4326", r.crs, always_xy=True)
    rr, cc = rasterio.transform.rowcol(r.transform, *to.transform(passes.lon.values, passes.lat.values))
    cell = np.asarray(rr).astype(np.int64) * r.width + np.asarray(cc)
order = np.lexsort((-passes.depth_m.to_numpy(), cell))
first = np.zeros(N, bool)
first[order[np.r_[True, cell[order][1:] != cell[order][:-1]]]] = True
UNI = first & trop

# database: census plus connectivity of tropical passes
db = passes.drop(columns=[c for c in passes.columns if c.startswith("gateway") or c in ("isolation_km_window", "block",
                                                                                         "core")])
db.insert(0, "apd_id", [f"APD{i + 1:06d}" for i in range(N)])
db["tropics"] = trop.astype(int)
db["cell_1km"] = np.where(trop, cell, -1)
for c, lab in enumerate(CFG_LAB):
    tag = "" if c == 0 else "_" + lab.replace(" ", "").replace(".", "")
    for j, l in enumerate(LECT):
        p = LECT_LAB[l][:5]
        db[f"{p}_med{tag}"] = MED[c, j].round(4)
        if c == 0:
            db[f"{p}_min"], db[f"{p}_max"] = MIN[c, j].round(4), MAX[c, j].round(4)
for j, l in enumerate(LECT):
    p = LECT_LAB[l][:5]
    db[f"{p}_robust"] = np.select([ROB_HI[j], ROB_LO[j]], ["above 1", "below 1"], "depends")
    db.loc[np.isnan(MED[:, j]).any(0), f"{p}_robust"] = ""
db.to_csv(f"{DB}/APD_passes.csv.gz", index=False, float_format="%.6g")

# ---- census (S1) and validation (S2 to S5) --------------------------------------------------------------------
t = passes.assign(zone=np.where(trop, "tropics", "rest of the Andes"), d100=passes.depth_m >= 100)
s1 = (t.groupby(["zone", "gmba", "rank"]).agg(passes=("z_m", "size"), depth_100m=("d100", "sum"),
                                              median_z_m=("z_m", "median"), median_depth_m=("depth_m", "median"))
      .reset_index())
tot = t.groupby(["zone", "rank"]).agg(passes=("z_m", "size"), depth_100m=("d100", "sum"), median_z_m=("z_m", "median"),
                                      median_depth_m=("depth_m", "median")).reset_index().assign(gmba="all units")
s1 = pd.concat([s1, tot]).sort_values(["zone", "gmba", "rank"], ascending=[False, True, True])
s1.round(1).to_csv(f"{SUP}/S1_census_by_unit.csv", index=False)
V = "data/interim/censo/validacion"
pd.read_csv(f"{V}/O5_4_pasos_control.csv").to_csv(f"{SUP}/S2_control_passes_census.csv", index=False)
o1 = pd.read_csv(f"{V}/O5_1_proporciones.csv").rename(columns={"Unnamed: 0": "group"})
o2 = pd.read_csv(f"{V}/O5_2_identidad.csv")
o1.merge(o2, on="group", how="outer").round(3).to_csv(f"{SUP}/S3_osm_nodes_identity.csv", index=False)
pd.read_csv(f"{V}/O5_3_alturas_referencia.csv").round(2).to_csv(f"{SUP}/S4_reference_heights.csv", index=False)
pd.read_csv(f"{V}/O5_3_error_por_zona.csv").rename(columns={"Unnamed: 0": "reference"}).round(2).to_csv(
    f"{SUP}/S5_height_error_by_zone.csv", index=False)
pd.read_csv(f"{V}/O5_1_sin_pass_1km.csv").to_csv(f"{SUP}/S6_osm_nodes_without_pass.csv", index=False)

# ---- montane limit and envelopes (S7) -------------------------------------------------------------------------
I = f"{D}/insumos"
t2 = pd.read_csv(f"{I}/t2000.csv")
env = pd.read_csv(f"{I}/envolventes.csv")
env = env[np.isclose(env.lat * 2, np.round(env.lat * 2))]                    # every 0.5 degree, to keep it readable
env = env.pivot_table(index=["var", "lat"], columns="gremio", values=["p5", "p95"]).reset_index()
env.columns = ["var", "lat", "p5_montane", "p5_lowland", "p95_montane", "p95_lowland"]
env = env.merge(t2, on="lat", how="left").merge(pd.read_csv(f"{I}/margenes.csv"), on="var")
env["var"] = env["var"].map(lambda v: f"bio{v:02d}")
env = env[(env.lat <= LN) & (env.lat >= LS)].sort_values(["var", "lat"], ascending=[True, False])
env[["var", "lat", "t2000", "iqr", "p5_lowland", "p95_lowland", "p5_montane", "p95_montane"]].round(3).to_csv(
    f"{SUP}/S7_envelopes.csv", index=False)

# ---- tiles: latitudinal profiles over every tile (S8) ------------------------------------------------------------
tt = []
for f in sorted(glob.glob(f"{D}/teselas/teselas_s*_bio??*.csv")):
    b = os.path.basename(f)[:-4].split("_")
    side, var = int(b[1][1:]), b[2]
    q = float(b[3][1:]) if len(b) > 3 else 0.1
    tt.append(pd.read_csv(f).assign(var=var, q=q, lado_km=side))
tt = pd.concat(tt)
tt = tt[tt.resistencia.str.endswith("env") | tt.resistencia.str.contains("_env_q")]
tt["cfg"] = [CFG_LAB[CFG.index((q, s))] for q, s in zip(tt.q.round(3), tt.lado_km)]
tt["reading"] = tt.lectura.map(LECT_LAB)
tt["band"] = np.floor(tt.lat / 0.5) * 0.5 + 0.25
pr = (tt[(tt.lat <= LN) & (tt.lat >= LS)].groupby(["cfg", "reading", "var", "band"]).conductancia_rel
      .agg(tiles="size", median="median", p25=lambda x: x.quantile(0.25), p75=lambda x: x.quantile(0.75))
      .reset_index().sort_values(["cfg", "reading", "var", "band"], ascending=[True, True, True, False]))
pr.round(4).to_csv(f"{SUP}/S8_tile_conductance_profiles.csv", index=False)

# ---- tiles that cross the continental divide (S17, S18) and the zone conductance of every cell -------------------
# A tile crosses the divide when continental passes lie in its central half along the direction of the reading, that is,
# between its two walls (C5 puts the walls in the outer quarters). A profile of every tile mixes in tiles that never
# touch the cordillera; these are the ones that measure how hard it is to cross it.
with rasterio.open(f"{D}/insumos/elev.tif") as r:
    PROF, RH, RW = r.profile, r.height, r.width
    crr, ccc = rasterio.transform.rowcol(r.transform, *to.transform(passes.lon[passes["rank"] == "continental"].values,
                                                                    passes.lat[passes["rank"] == "continental"].values))
div = np.zeros((RH, RW), np.int32)
div[np.clip(crr, 0, RH - 1), np.clip(ccc, 0, RW - 1)] = 1
INT = np.pad(div.cumsum(0).cumsum(1), ((1, 0), (1, 0)))


def n_passes(r0, r1, c0, c1):
    """Continental passes inside each rectangle, from the integral image."""
    r0, r1, c0, c1 = (np.clip(v, 0, m) for v, m in ((r0, RH), (r1, RH), (c0, RW), (c1, RW)))
    return INT[r1, c1] - INT[r0, c1] - INT[r1, c0] + INT[r0, c0]


S = tt.lado_km.to_numpy()
r0, c0 = tt.row.to_numpy() - S // 2, tt.col.to_numpy() - S // 2
low = (tt.lectura == "baja").to_numpy()                                  # transverse: walls west and east
tt["crosses_divide"] = np.where(low, n_passes(r0, r0 + S, c0 + S // 4, c0 + S - S // 4),
                                n_passes(r0 + S // 4, r0 + S - S // 4, c0, c0 + S)) > 0
GRP = {"temperature": VARS[:11], "precipitation": VARS[11:]}
dv = tt[tt.crosses_divide & (tt.lat <= LN) & (tt.lat >= LS)]
w = dv.pivot_table(index=["cfg", "reading", "band", "row", "col"], columns="var", values="conductancia_rel")
s17 = []
for g, vs in GRP.items():
    a = w[vs].to_numpy()
    s17.append(pd.DataFrame({"group": g, "median": np.median(a, axis=1), "most_limiting": a.min(axis=1),
                             "least_limiting": a.max(axis=1)}, index=w.index))
s17 = pd.concat(s17).groupby(["cfg", "reading", "group", "band"]).agg(
    tiles=("median", "size"), median=("median", "median"), most_limiting=("most_limiting", "median"),
    least_limiting=("least_limiting", "median")).reset_index()
s17["cfg"] = pd.Categorical(s17.cfg, CFG_LAB)
s17.sort_values(["cfg", "reading", "group", "band"], ascending=[True, True, True, False]).round(4).to_csv(
    f"{SUP}/S17_divide_tile_profiles.csv", index=False)
s18 = []                                                                 # which variable limits its group (main run)
m = w.loc["main"]
for g, vs in GRP.items():
    lim = pd.Series(np.array(vs)[m[vs].to_numpy().argmin(axis=1)], index=m.index, name="var")
    s = lim.groupby(level="reading").value_counts(normalize=True).rename("share_of_tiles").reset_index()
    s18.append(s.rename(columns={"index": "var", 0: "var"}).assign(group=g))
pd.concat(s18)[["reading", "group", "var", "share_of_tiles"]].sort_values(
    ["reading", "group", "share_of_tiles"], ascending=[True, True, False]).round(4).to_csv(
    f"{SUP}/S18_limiting_variable.csv", index=False)
# zone conductance of each 1 km cell: per tile the median relative conductance over the temperature variables (BIO1
# to BIO11, the measure of the temperature line of Fig. 3c and 3d), then the mean over the tiles of that direction
# that cover the cell (C7), each weighted by a Hann window centred on the tile, so the value changes gradually instead
# of in steps of S/4 and a cell takes most of its value from the tiles centred near it
PROF.update(count=len(CFG), dtype="float32", nodata=np.nan, compress="deflate")
for lect in LECT:
    z = []
    for q, side in CFG:
        t = tt[(tt.lectura == lect) & (tt.lado_km == side) & (tt.q.round(3) == q) & tt["var"].isin(GRP["temperature"])]
        t = t.pivot_table(index=["row", "col"], columns="var", values="conductancia_rel").dropna().median(axis=1)
        h = np.sin(np.pi * (np.arange(side) + 0.5) / side) ** 2              # Hann window, as in 50
        win = np.outer(h, h)
        su, cnt = np.zeros((RH, RW)), np.zeros((RH, RW))
        for (rw, cl), g in t.items():
            a, b = rw - side // 2, cl - side // 2
            ww = win[max(-a, 0):RH - a, max(-b, 0):RW - b]
            su[max(a, 0):a + side, max(b, 0):b + side] += g * ww
            cnt[max(a, 0):a + side, max(b, 0):b + side] += ww
        z.append(np.where(cnt > 0, su / np.where(cnt > 0, cnt, 1), np.nan))
    with rasterio.open(f"{D}/zona_temperatura_{lect}.tif", "w", **PROF) as o:
        o.write(np.stack(z).astype("float32"))
# normalized current of each cell summarized by group, as the passes are summarized over the 19 variables: the median
# over the variables of the group, undefined where any of them is (uniform current below U_MIN)
os.makedirs(f"{D}/grupos", exist_ok=True)
for q, side in CFG:
    tag = f"s{side}_env" + ("" if q == 0.1 else f"_q{q:g}")
    for lect in LECT:
        for g, vs in GRP.items():
            st = None
            for k, v in enumerate(vs):
                with rasterio.open(f"{D}/teselas/{v}_{lect}_{tag}_normalizada.tif") as r:
                    a, gp = r.read(1), r.profile
                with rasterio.open(f"{D}/teselas/{v}_{lect}_{tag}_corriente.tif") as r:
                    c = r.read(1)
                with np.errstate(all="ignore"):
                    a = np.where((a != -9999) & (c != -9999) & (c / a >= U_MIN), a, np.nan).astype(np.float32)
                st = np.empty((len(vs),) + a.shape, np.float32) if st is None else st
                st[k] = a
            gp.update(dtype="float32", nodata=np.nan, compress="deflate")
            with rasterio.open(f"{D}/grupos/{g[:4]}_{lect}_{tag}.tif", "w", **gp) as o:
                o.write(np.median(st, axis=0), 1)
    print("group current written:", tag, flush=True)
print("zone conductance written; tiles that cross the divide:",
      tt[tt.cfg == "main"].groupby("reading").crosses_divide.sum().to_dict())
# ---- bands of local maxima of the divide profiles (S19), numbered as in Fig. 3 ---------------------------------------
# Each band is a run of 0.5° bins of S17 (main configuration) around a local extreme of the temperature median: seven
# maxima of the transverse profile (1 to 7, openings across the cordillera) and three minima of the longitudinal profile
# (8 to 10, barriers along it), north to south. Edges are bin edges.
BANDS19 = [("transverse", 8.5, 8.0), ("transverse", 6.5, 5.5), ("transverse", 3.5, 3.0), ("transverse", 1.5, 1.0),
           ("transverse", -4.0, -6.5), ("transverse", -13.0, -14.0), ("transverse", -16.0, -16.5),
           ("longitudinal", 8.5, 7.0), ("longitudinal", 6.0, 2.0), ("longitudinal", -3.5, -6.0)]
m17 = s17[s17.cfg == "main"]
ctl19 = pd.read_csv("data/refs/pasos_control.csv")
cont19 = db[(db["rank"] == "continental") & (db.depth_m >= 100)]
s19 = []
for k, (rd, n, s) in enumerate(BANDS19, 1):
    r = {"no": k, "reading": rd, "north": n, "south": s, "width_deg": n - s, "bins": int((n - s) / 0.5)}
    for g in ("temperature", "precipitation"):
        q = m17[(m17.reading == rd) & (m17.group == g)]
        inb, nb = q[(q.band < n) & (q.band > s)], q[((q.band > n) & (q.band < n + 1)) | ((q.band < s) & (q.band > s - 1))]
        if g == "temperature":
            r["tiles_min"], r["tiles_max"] = int(inb.tiles.min()), int(inb.tiles.max())
        r[f"{g[:4]}_median"] = inb["median"].mean()
        r[f"{g[:4]}_extreme"] = inb["median"].max() if rd == "transverse" else inb["median"].min()
        r[f"{g[:4]}_most_limiting"], r[f"{g[:4]}_least_limiting"] = inb.most_limiting.mean(), inb.least_limiting.mean()
        r[f"{g[:4]}_median_neighbours"] = nb["median"].mean()          # the 1° on either side of the band
    p = cont19[(cont19.lat <= n) & (cont19.lat >= s)].sort_values("z_m").head(1)
    r.update({"lowest_pass": p.apd_id.iat[0], "lowest_lat": p.lat.iat[0], "lowest_lon": p.lon.iat[0],
              "lowest_z_m": p.z_m.iat[0], "lowest_depth_m": p.depth_m.iat[0]} if len(p) else {})
    r["control_passes"] = "; ".join(ctl19[(ctl19.lat <= n) & (ctl19.lat >= s)].name)
    s19.append(r)
pd.DataFrame(s19).round(3).to_csv(f"{SUP}/S19_profile_bands.csv", index=False)

# ---- control of Huancabamba and the Tachira depression (S9, S20), on the tiles that cross the divide -------------
# each tile is summarized by the median of its group of variables, as in S17; a band is the median over its tiles
gt = tt[tt.crosses_divide].pivot_table(index=["cfg", "reading", "row", "col", "lat", "lon", "band"], columns="var",
                                       values="conductancia_rel").reset_index()
for g, vs in GRP.items():
    gt[g] = np.median(gt[vs].to_numpy(), axis=1)                 # nan if a variable of the group is missing
gl = gt.melt(id_vars=["cfg", "reading", "lat", "lon", "band"], value_vars=list(GRP), var_name="group",
             value_name="conductance").dropna(subset=["conductance"])
REG = {"Huancabamba": lambda x: x[(x.lat <= -2) & (x.lat >= -10)],
       "Tachira": lambda x: x[(x.lon >= -74) & (x.lon <= -70) & (x.lat >= 5)]}
ZONE = {"Huancabamba": {"zone": (-7.5, -4.0), "north": (-4.0, -2.0), "south": (-9.5, -7.5)},
        "Tachira": {"zone": (7.5, 8.0), "south": (5.0, 7.5)}}
c9, s20 = [], []
for lab, f in REG.items():
    b = f(gl).groupby(["reading", "group", "cfg", "band"]).conductance.agg(tiles="size", median="median").reset_index()
    lo, hi = ZONE[lab]["zone"]
    b["in_zone"] = (b.band > lo) & (b.band < hi)
    c9.append(b.assign(control=lab))
    for (rd, g, c), x in b.groupby(["reading", "group", "cfg"]):
        r = {"reading": rd, "control": lab, "group": g, "cfg": c}
        for k, (a0, a1) in ZONE[lab].items():
            r[k] = x[(x.band > a0) & (x.band < a1)]["median"].mean()
        r["lowest_band"] = x.band.iat[x["median"].to_numpy().argmin()]
        r["lowest_band_in_zone"] = bool(x.in_zone.iat[x["median"].to_numpy().argmin()])
        r["zone_below_neighbours"] = r["zone"] < min(r[k] for k in ZONE[lab] if k != "zone")
        s20.append(r)
c9 = pd.concat(c9)
c9 = c9[c9.reading == "longitudinal"].drop(columns="reading")      # S9 is the control of the Methods
c9["cfg"] = pd.Categorical(c9.cfg, CFG_LAB)
c9[["control", "group", "cfg", "band", "tiles", "median", "in_zone"]].sort_values(
    ["control", "group", "cfg", "band"], ascending=[True, False, True, False]).round(3).to_csv(
    f"{SUP}/S9_control_huancabamba_tachira.csv", index=False)
s20 = pd.DataFrame(s20)
s20["cfg"] = pd.Categorical(s20.cfg, CFG_LAB)
s20.sort_values(["reading", "control", "group", "cfg"], ascending=[True, True, False, True]).round(3).to_csv(
    f"{SUP}/S20_control_summary.csv", index=False)

# every variable, main run: mean of the band medians inside each control zone (as in C9, tiles that cross the divide)
# and in the 2 degrees next to it on each side where tiles exist (S15)
s15 = []
for rd in ("longitudinal", "transverse"):
    m = tt[tt.crosses_divide & (tt.cfg == "main") & (tt.reading == rd)]
    for lab, s, parts in (("Huancabamba", m, {"zone": (-7.5, -4.0), "north": (-4.0, -2.0), "south": (-9.5, -7.5)}),
                          ("Tachira", m[(m.lon >= -74) & (m.lon <= -70)], {"zone": (7.5, 8.0), "south": (5.0, 7.5)})):
        b = s.groupby(["var", "band"]).conductancia_rel.median().reset_index()
        for k, (lo, hi) in parts.items():
            x = b[(b.band > lo) & (b.band < hi)].groupby("var").conductancia_rel.mean()
            s15 += [(rd, lab, v, k, x.get(v, np.nan)) for v in VARS]
s15 = pd.DataFrame(s15, columns=["reading", "control", "var", "part", "conductance"]).pivot_table(
    index=["reading", "control", "var"], columns="part", values="conductance").reset_index()
s15["lowest_in_zone"] = s15["zone"] < s15[["north", "south"]].min(axis=1)
s15.round(3).to_csv(f"{SUP}/S15_controls_by_variable.csv", index=False)

# ---- which variables change the result (S10, S11): Spearman between variables, passes and tiles, main run --------
rows = []
for j, l in enumerate(LECT):
    x = pd.DataFrame(A[0, j][:, UNI].T, columns=VARS).dropna()
    rc = x.rank().corr()
    tl = tt[(tt.cfg == "main") & (tt.lectura == l)].pivot_table(index=["row", "col"], columns="var",
                                                               values="conductancia_rel").dropna()
    rt = tl.rank().corr()
    for a in VARS:
        for b in VARS:
            rows.append((LECT_LAB[l], a, b, rc.loc[a, b], rt.loc[a, b]))
    print(f"{LECT_LAB[l]}: {len(x):,} cells with 19 variables, {len(tl):,} tiles")
sp = pd.DataFrame(rows, columns=["reading", "var_a", "var_b", "rho_passes", "rho_tiles"])
sp.round(3).to_csv(f"{SUP}/S10_spearman_variables.csv", index=False)
# per variable: median rho with the other 18, and rho with the pass median over the 19
s11 = []
for j, l in enumerate(LECT):
    med = MED[0, j][UNI]
    for k, v in enumerate(VARS):
        x = A[0, j, k][UNI]
        ok = np.isfinite(x) & np.isfinite(med)
        o = sp[(sp.reading == LECT_LAB[l]) & (sp.var_a == v) & (sp.var_b != v)]
        s11.append((LECT_LAB[l], v, o.rho_passes.median(), o.rho_passes.min(), o.rho_tiles.median(),
                    pd.Series(x[ok]).rank().corr(pd.Series(med[ok]).rank())))
pd.DataFrame(s11, columns=["reading", "var", "median_rho_other_vars_passes", "min_rho_other_vars_passes",
                           "median_rho_other_vars_tiles", "rho_with_pass_median"]).round(3).to_csv(
    f"{SUP}/S11_variable_agreement.csv", index=False)

# ---- sensitivity of the pass medians (S12) -----------------------------------------------------------------------
s12 = []
for j, l in enumerate(LECT):
    base = MED[0, j][UNI]
    top = base >= np.nanquantile(base, 0.9)
    for c in range(1, len(CFG)):
        x = MED[c, j][UNI]
        ok = np.isfinite(base) & np.isfinite(x)
        s12.append((LECT_LAB[l], CFG_LAB[c], ok.sum(), np.isfinite(base).sum() - ok.sum(),
                    pd.Series(base[ok]).rank().corr(pd.Series(x[ok]).rank()),
                    np.mean((x[ok & top] >= np.nanquantile(x, 0.9))),
                    np.mean((base[ok] > 1) == (x[ok] > 1))))
pd.DataFrame(s12, columns=["reading", "configuration", "cells_both", "cells_lost", "rho_with_main",
                           "top_decile_kept", "same_side_of_1"]).round(3).to_csv(f"{SUP}/S12_sensitivity.csv",
                                                                                  index=False)

# ---- passes by unit and latitude (S13, S14) -----------------------------------------------------------------------
u = db[UNI].copy()
s13 = []
for l in LECT:
    p = LECT_LAB[l][:5]
    x = u[u[f"{p}_med"].notna()]
    for g in (x.groupby("gmba"), x.assign(gmba="all units").groupby("gmba")):
        s13.append(g.agg(cells=(f"{p}_med", "size"), median=(f"{p}_med", "median"),
                         p90=(f"{p}_med", lambda x: x.quantile(0.9)),
                         n_robust_above_1=(f"{p}_robust", lambda x: (x == "above 1").sum()),
                         n_robust_below_1=(f"{p}_robust", lambda x: (x == "below 1").sum()),
                         robust_above_1=(f"{p}_robust", lambda x: (x == "above 1").mean()),
                         robust_below_1=(f"{p}_robust", lambda x: (x == "below 1").mean()))
                   .reset_index().assign(reading=LECT_LAB[l]))
pd.concat(s13)[["reading", "gmba", "cells", "median", "p90", "n_robust_above_1", "n_robust_below_1",
                 "robust_above_1", "robust_below_1"]].round(3).to_csv(
    f"{SUP}/S13_by_unit.csv", index=False)
u["band"] = np.floor(u.lat) + 0.5
s14 = []
for l in LECT:
    p = LECT_LAB[l][:5]
    g = u[u[f"{p}_med"].notna()].groupby("band")
    s14.append(g.agg(cells=(f"{p}_med", "size"), median=(f"{p}_med", "median"),
                     p90=(f"{p}_med", lambda x: x.quantile(0.9)),
                     robust_above_1=(f"{p}_robust", lambda x: (x == "above 1").mean()))
               .reset_index().assign(reading=LECT_LAB[l]))
pd.concat(s14)[["reading", "band", "cells", "median", "p90", "robust_above_1"]].sort_values(
    ["reading", "band"], ascending=[False, False]).round(3).to_csv(f"{SUP}/S14_by_latitude.csv", index=False)

# ---- passes by height band (S16): the two readings against the height of the pass -------------------------------------
u["z_band"] = (np.floor(u.z_m / 500) * 500).clip(lower=0, upper=5500)          # 5500 m holds every pass above it
s16 = []
for l in LECT:
    p = LECT_LAB[l][:5]
    g = u[u[f"{p}_med"].notna()].groupby("z_band")
    s16.append(g.agg(cells=(f"{p}_med", "size"), median=(f"{p}_med", "median"),
                     p10=(f"{p}_med", lambda x: x.quantile(0.1)), p90=(f"{p}_med", lambda x: x.quantile(0.9)),
                     robust_above_1=(f"{p}_robust", lambda x: (x == "above 1").mean()),
                     robust_below_1=(f"{p}_robust", lambda x: (x == "below 1").mean()))
               .reset_index().assign(reading=LECT_LAB[l]))
pd.concat(s16)[["reading", "z_band", "cells", "median", "p10", "p90", "robust_above_1", "robust_below_1"]].round(3)    .to_csv(f"{SUP}/S16_by_height.csv", index=False)

# ---- Table 1: per unit and reading, the cell with the highest median that stay above 1 in all five ---------
t1 = []
for l in LECT:
    p = LECT_LAB[l][:5]
    x = u[u[f"{p}_robust"] == "above 1"]
    for g, d in x.groupby("gmba"):
        d = d.nlargest(1, f"{p}_med")
        cfg_passes = [f"{p}_med"] + [f"{p}_med_" + c.replace(" ", "").replace(".", "") for c in CFG_LAB[1:]]
        t1.append(d.assign(reading=LECT_LAB[l], med=d[f"{p}_med"], vmin=d[f"{p}_min"], vmax=d[f"{p}_max"],
                           cfg_lo=d[cfg_passes].min(axis=1), cfg_hi=d[cfg_passes].max(axis=1)))
t1 = pd.concat(t1)[["reading", "gmba", "apd_id", "lat", "lon", "z_m", "depth_m", "rank", "med", "vmin", "vmax",
                    "cfg_lo", "cfg_hi"]]
t1.round({"lat": 3, "lon": 3, "z_m": 0, "depth_m": 0, "med": 2, "vmin": 2, "vmax": 2, "cfg_lo": 2,
          "cfg_hi": 2}).to_csv(f"{TAB}/Table1_top_passes.csv", index=False)

# ---- Table 2: control passes of the tropics, nearest pass ------------------------------------------------------------
ctl = pd.read_csv("data/refs/pasos_control.csv")
ctl = ctl[(ctl.in_area == 1) & (ctl.lat <= LN) & (ctl.lat >= LS)]


def hav(la, lo, LA, LO):
    la, lo, LA, LO = map(np.radians, (la, lo, LA, LO))
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.sin((LA - la) / 2) ** 2 + np.cos(la) * np.cos(LA) * np.sin((LO - lo) / 2) ** 2))


t2r = []
for _, c in ctl.iterrows():
    dk = hav(c.lat, c.lon, passes.lat.values, passes.lon.values)
    i = int(np.argmin(dk))
    r = {"no": len(t2r) + 1, "pass": c["name"], "system": c.system, "lat": c.lat, "lon": c.lon, "apd_id": db.apd_id[i],
         "census_km": dk[i], "z_m": passes.z_m[i], "depth_m": passes.depth_m[i], "rank": passes["rank"][i]}
    for j, l in enumerate(LECT):
        p = LECT_LAB[l][:5]
        r.update({f"{p}_med": MED[0, j, i], f"{p}_min": MIN[0, j, i], f"{p}_max": MAX[0, j, i],
                  f"{p}_cfg_lo": np.nanmin(MED[:, j, i]), f"{p}_cfg_hi": np.nanmax(MED[:, j, i])})
    t2r.append(r)
pd.DataFrame(t2r).round(2).to_csv(f"{TAB}/Table2_control_passes.csv", index=False)

# ---- numbers quoted in the text ------------------------------------------------------------------------------------
print(f"passes in the census {N:,}; tropics {trop.sum():,}; cells {UNI.sum():,}")
for j, l in enumerate(LECT):
    v = MED[0, j][UNI]
    print(f"{LECT_LAB[l]}: cells with 19 vars {np.isfinite(v).sum():,}, median {np.nanmedian(v):.3f}, "
          f"robust above 1 {ROB_HI[j][UNI].sum():,}, below 1 {ROB_LO[j][UNI].sum():,}")
print("->", TAB, SUP, DB)
