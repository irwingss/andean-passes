"""51_sintesis_conectividad.py - Synthesis of the present climatic connectivity (PRD C8 and C9), without inference (D21).
Reads data/interim/censo/cols_andes.csv, the col values written by 49_cols_corriente.py and the tile tables of
50_conectividad_teselas.py for one run, and writes:
  data/interim/conectividad/<run>/normalizada.npy   normalized current, array [configuration, reading, variable, col]
  outputs/database/APD_cols.csv.gz                 the census with, for tropical cols, the median and range over the 19
                                                   variables of the normalized current in each reading and configuration
  outputs/tables/Table1_top_cols.csv, Table2_control_passes.csv
  outputs/supp/tables/S*.csv                       census, validation, envelopes, profiles, variables and sensitivity
python -P scripts/51_sintesis_conectividad.py principal
Configurations: q 0.1 with 200 km tiles (main), q 0.05 and 0.2 with 200 km, and 100 and 400 km with q 0.1. A statement
about a col holds when it holds in all five (C8). Cols in the same 1 km cell share their values; rankings and
correlations keep one col per cell, the deepest. The normalized current is left undefined where the current under
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

cols = pd.read_csv("data/interim/censo/cols_andes.csv", dtype={"hybas_term_a": str, "hybas_term_b": str})
N = len(cols)
trop = ((cols.lat <= LN) & (cols.lat >= LS)).to_numpy()

# normalized current per configuration, reading, variable and col (cached: reading the 63 million rows takes minutes)
NPY = f"{D}/normalizada.npy"
U_MIN = 1e-6      # below this current under uniform resistance (of the 1 A per tile) the ratio is numerically unstable
if not os.path.exists(NPY) or os.path.getmtime(NPY) < os.path.getmtime(f"{D}/cols_corriente.csv.gz"):
    A = np.full((len(CFG), 2, 19, N), np.nan, np.float32)
    ci = {c: i for i, c in enumerate(CFG)}
    for ch in pd.read_csv(f"{D}/cols_corriente.csv.gz", chunksize=5_000_000,
                          dtype={"col_id": np.int32, "var": str, "lectura": str, "q": np.float64,
                                 "lado_km": np.int32, "corriente": np.float32, "normalizada": np.float32}):
        k = np.array([ci[(q, s)] for q, s in zip(ch.q.round(3), ch.lado_km)])
        A[k, (ch.lectura == "alta").to_numpy().astype(int), ch["var"].str[3:].astype(int).to_numpy() - 1,
          ch.col_id.to_numpy()] = np.where(ch.corriente / ch.normalizada >= U_MIN, ch.normalizada, np.nan)
    np.save(NPY, A)
A = np.load(NPY)
with np.errstate(all="ignore"):
    MED = np.nanmedian(A, axis=2)                      # [cfg, reading, col], median over the 19 variables
    MIN, MAX = np.nanmin(A, axis=2), np.nanmax(A, axis=2)
NV = np.isfinite(A).sum(axis=2)
MED[NV < 19] = np.nan                                  # a col needs the 19 variables
ROB_HI = np.all(MED > 1, axis=0)                       # concentrates current in all five configurations
ROB_LO = np.all(MED < 1, axis=0)                       # diverts current in all five

# 1 km cell of every col, and one col per cell (the deepest) for rankings and correlations
with rasterio.open(f"{D}/insumos/elev.tif") as r:
    to = Transformer.from_crs("EPSG:4326", r.crs, always_xy=True)
    rr, cc = rasterio.transform.rowcol(r.transform, *to.transform(cols.lon.values, cols.lat.values))
    cell = np.asarray(rr).astype(np.int64) * r.width + np.asarray(cc)
order = np.lexsort((-cols.depth_m.to_numpy(), cell))
first = np.zeros(N, bool)
first[order[np.r_[True, cell[order][1:] != cell[order][:-1]]]] = True
UNI = first & trop

# database: census plus connectivity of tropical cols
db = cols.drop(columns=[c for c in cols.columns if c.startswith("gateway") or c in ("isolation_km_window", "block",
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
db.to_csv(f"{DB}/APD_cols.csv.gz", index=False, float_format="%.6g")

# ---- census (S1) and validation (S2 to S5) --------------------------------------------------------------------
t = cols.assign(zone=np.where(trop, "tropics", "rest of the Andes"), d100=cols.depth_m >= 100)
s1 = (t.groupby(["zone", "gmba", "rank"]).agg(cols=("z_m", "size"), depth_100m=("d100", "sum"),
                                              median_z_m=("z_m", "median"), median_depth_m=("depth_m", "median"))
      .reset_index())
tot = t.groupby(["zone", "rank"]).agg(cols=("z_m", "size"), depth_100m=("d100", "sum"), median_z_m=("z_m", "median"),
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
pd.read_csv(f"{V}/O5_1_sin_col_1km.csv").to_csv(f"{SUP}/S6_osm_nodes_without_col.csv", index=False)

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

# ---- tiles: profiles (S8) and the control of Huancabamba and the Tachira depression (S9) -----------------------
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
c9 = []
m = tt[(tt["var"] == "bio01") & (tt.lectura == "alta")]
for lab, s, zone in (("Huancabamba", m[(m.lat <= -2) & (m.lat >= -10)], (-7.5, -4)),
                     ("Tachira", m[(m.lon >= -74) & (m.lon <= -70) & (m.lat >= 5)], (7.5, 8))):
    p = s.groupby(["cfg", "band"]).conductancia_rel.agg(tiles="size", median="median").reset_index()
    c9.append(p.assign(control=lab, in_zone=(p.band >= zone[0]) & (p.band <= zone[1])))
c9 = pd.concat(c9)
c9["cfg"] = pd.Categorical(c9.cfg, CFG_LAB)
c9.sort_values(["control", "cfg", "band"], ascending=[True, True, False]).round(3).to_csv(
    f"{SUP}/S9_control_huancabamba_tachira.csv", index=False)

# every variable, main run: mean of the band medians inside each control zone and in the latitudes next to it (S15)
s15 = []
for rd in ("longitudinal", "transverse"):
    m = tt[(tt.cfg == "main") & (tt.reading == rd)]
    for lab, s, parts in (("Huancabamba", m, {"zone": (-6.0, -4.0), "north": (-4.0, -2.0), "south": (-8.5, -6.5)}),
                          ("Tachira", m[(m.lon >= -74) & (m.lon <= -70)], {"zone": (7.5, 8.0), "south": (5.0, 7.5)})):
        b = s.groupby(["var", "band"]).conductancia_rel.median().reset_index()
        for k, (lo, hi) in parts.items():
            x = b[(b.band > lo) & (b.band < hi)].groupby("var").conductancia_rel.mean()
            s15 += [(rd, lab, v, k, x.get(v, np.nan)) for v in VARS]
s15 = pd.DataFrame(s15, columns=["reading", "control", "var", "part", "conductance"]).pivot_table(
    index=["reading", "control", "var"], columns="part", values="conductance").reset_index()
s15["lowest_in_zone"] = s15["zone"] < s15[["north", "south"]].min(axis=1)
s15.round(3).to_csv(f"{SUP}/S15_controls_by_variable.csv", index=False)

# ---- which variables change the result (S10, S11): Spearman between variables, cols and tiles, main run --------
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
sp = pd.DataFrame(rows, columns=["reading", "var_a", "var_b", "rho_cols", "rho_tiles"])
sp.round(3).to_csv(f"{SUP}/S10_spearman_variables.csv", index=False)
# per variable: median rho with the other 18, and rho with the col median over the 19
s11 = []
for j, l in enumerate(LECT):
    med = MED[0, j][UNI]
    for k, v in enumerate(VARS):
        x = A[0, j, k][UNI]
        ok = np.isfinite(x) & np.isfinite(med)
        o = sp[(sp.reading == LECT_LAB[l]) & (sp.var_a == v) & (sp.var_b != v)]
        s11.append((LECT_LAB[l], v, o.rho_cols.median(), o.rho_cols.min(), o.rho_tiles.median(),
                    pd.Series(x[ok]).rank().corr(pd.Series(med[ok]).rank())))
pd.DataFrame(s11, columns=["reading", "var", "median_rho_other_vars_cols", "min_rho_other_vars_cols",
                           "median_rho_other_vars_tiles", "rho_with_col_median"]).round(3).to_csv(
    f"{SUP}/S11_variable_agreement.csv", index=False)

# ---- sensitivity of the col medians (S12) -----------------------------------------------------------------------
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

# ---- cols by unit and latitude (S13, S14) -----------------------------------------------------------------------
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

# ---- Table 1: per unit and reading, the cell with the highest median that stay above 1 in all five ---------
t1 = []
for l in LECT:
    p = LECT_LAB[l][:5]
    x = u[u[f"{p}_robust"] == "above 1"]
    for g, d in x.groupby("gmba"):
        d = d.nlargest(1, f"{p}_med")
        cfg_cols = [f"{p}_med"] + [f"{p}_med_" + c.replace(" ", "").replace(".", "") for c in CFG_LAB[1:]]
        t1.append(d.assign(reading=LECT_LAB[l], med=d[f"{p}_med"], vmin=d[f"{p}_min"], vmax=d[f"{p}_max"],
                           cfg_lo=d[cfg_cols].min(axis=1), cfg_hi=d[cfg_cols].max(axis=1)))
t1 = pd.concat(t1)[["reading", "gmba", "apd_id", "lat", "lon", "z_m", "depth_m", "rank", "med", "vmin", "vmax",
                    "cfg_lo", "cfg_hi"]]
t1.round({"lat": 3, "lon": 3, "z_m": 0, "depth_m": 0, "med": 2, "vmin": 2, "vmax": 2, "cfg_lo": 2,
          "cfg_hi": 2}).to_csv(f"{TAB}/Table1_top_cols.csv", index=False)

# ---- Table 2: control passes of the tropics, nearest col ------------------------------------------------------------
ctl = pd.read_csv("data/refs/pasos_control.csv")
ctl = ctl[(ctl.in_area == 1) & (ctl.lat <= LN) & (ctl.lat >= LS)]


def hav(la, lo, LA, LO):
    la, lo, LA, LO = map(np.radians, (la, lo, LA, LO))
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.sin((LA - la) / 2) ** 2 + np.cos(la) * np.cos(LA) * np.sin((LO - lo) / 2) ** 2))


t2r = []
for _, c in ctl.iterrows():
    dk = hav(c.lat, c.lon, cols.lat.values, cols.lon.values)
    i = int(np.argmin(dk))
    r = {"no": len(t2r) + 1, "pass": c["name"], "system": c.system, "lat": c.lat, "lon": c.lon, "apd_id": db.apd_id[i],
         "col_km": dk[i], "z_m": cols.z_m[i], "depth_m": cols.depth_m[i], "rank": cols["rank"][i]}
    for j, l in enumerate(LECT):
        p = LECT_LAB[l][:5]
        r.update({f"{p}_med": MED[0, j, i], f"{p}_min": MIN[0, j, i], f"{p}_max": MAX[0, j, i],
                  f"{p}_cfg_lo": np.nanmin(MED[:, j, i]), f"{p}_cfg_hi": np.nanmax(MED[:, j, i])})
    t2r.append(r)
pd.DataFrame(t2r).round(2).to_csv(f"{TAB}/Table2_control_passes.csv", index=False)

# ---- numbers quoted in the text ------------------------------------------------------------------------------------
print(f"cols in the census {N:,}; tropics {trop.sum():,}; cells {UNI.sum():,}")
for j, l in enumerate(LECT):
    v = MED[0, j][UNI]
    print(f"{LECT_LAB[l]}: cells with 19 vars {np.isfinite(v).sum():,}, median {np.nanmedian(v):.3f}, "
          f"robust above 1 {ROB_HI[j][UNI].sum():,}, below 1 {ROB_LO[j][UNI].sum():,}")
print("->", TAB, SUP, DB)
