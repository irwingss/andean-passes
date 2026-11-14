"""57_barreras_pasos.py - The pass of the census through which each barrier of 56_barreras_sintesis.py is crossed.
A barrier that cuts the continental divide is crossed along the crest by its lowest continental pass: the lowest
continental pass of the census inside the barrier (barreras_id.tif) is the pass of the barrier, and the barrier should
cut the montane belt from the first montane limit above the height of that pass. A barrier without continental passes
lies between ranges, a valley that the continental divide does not cross; for it, the floor of the gap is the cell of
highest resistance along the least-cost path of its links (55_barreras_grafo.py), and the census pass with a depth of
100 m or more that most links point to, nearest to their floors, is given with the median height of the floors.
For every barrier: support at each montane limit, the lowest limit at which seven or more variables (of 19) support it,
and the first limit above the height of its continental pass.
Writes outputs/supp/tables/S26_barrier_passes.csv (one row per barrier of Table S24) and
outputs/tables/Table3_barrier_passes.csv (the barriers of the design that seven or more variables support in some
configuration of the design).
python -P scripts/57_barreras_pasos.py principal
"""
import sys

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from scipy.spatial import cKDTree

RUN = sys.argv[1]
B = f"data/interim/conectividad/{RUN}/barreras"
SUP = "outputs/supp/tables"
HS = [("h500", 500), ("h1000", 1000), ("h1500", 1500), ("main", 2000), ("h2500", 2500), ("h3000", 3000),
      ("h3500", 3500), ("h4000", 4000), ("h4500", 4500)]
DESIGN = ["main", "q0.05", "q0.2", "a100", "a1000", "h1500", "h2500"]
NMIN = 7                                                    # variables that make a barrier well supported
CTL_KM = 3                                                  # a control pass within this distance names the pass
with rasterio.open(f"{B}/barreras_id.tif") as s:
    U, TR, CRS = s.read(1), s.transform, s.crs
with rasterio.open(f"data/interim/conectividad/{RUN}/insumos/elev.tif") as s:
    Z = s.read(1, masked=True).filled(np.nan)
to_m = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)


def cells(lon, lat):
    r, c = rasterio.transform.rowcol(TR, *to_m.transform(np.asarray(lon), np.asarray(lat)))
    return np.clip(np.asarray(r), 0, U.shape[0] - 1), np.clip(np.asarray(c), 0, U.shape[1] - 1)


apd = pd.read_csv("outputs/database/APD_passes.csv.gz",
                  usecols=["apd_id", "lat", "lon", "z_m", "depth_m", "rank", "gmba", "tropics"])
apd = apd[apd.tropics == 1].reset_index(drop=True)
apd["barrier_id"] = U[cells(apd.lon, apd.lat)]
deep = apd[apd.depth_m >= 100].reset_index(drop=True)
tree = cKDTree(np.column_stack(to_m.transform(deep.lon.values, deep.lat.values)))
ctl = pd.read_csv("data/refs/pasos_control.csv", encoding="utf-8")
ctree = cKDTree(np.column_stack(to_m.transform(ctl.lon.values, ctl.lat.values)))

s24 = pd.read_csv(f"{SUP}/S24_barriers.csv")
lk = pd.read_csv(f"{SUP}/S25_barrier_links.csv")
lk = lk[lk.barrier_id > 0].copy()
d, i = tree.query(np.column_stack(to_m.transform(lk.gap_lon.values, lk.gap_lat.values)))
lk["pass"], lk["pass_km"] = i, d / 1000
lk["gap_z_m"] = Z[cells(lk.gap_lon, lk.gap_lat)]


def control(p):
    cd, ci = ctree.query(to_m.transform(p.lon, p.lat))
    return ctl.name.iat[ci] if cd / 1000 <= CTL_KM else ""


rows = []
for _, b in s24.iterrows():
    r = {"no": b.no, "barrier_id": b.barrier_id, "defined_by": b.defined_by, "lat": b.lat, "lon": b.lon,
         "nearest_hazzi": b.nearest_hazzi, "hazzi_km": b.hazzi_km, "n_max": b.n_max}
    for c, h in HS:
        r[f"n_{h}"] = b[f"n_{c}"]
    on = [h for c, h in HS if b[f"n_{c}"] >= NMIN]
    r["lowest_limit_nmin"] = min(on) if on else np.nan
    inside = apd[apd.barrier_id == b.barrier_id]
    cont = inside[inside["rank"] == "continental"]
    r["passes_inside"], r["continental_inside"] = len(inside), len(cont)
    if len(cont):
        p = cont.loc[cont.z_m.idxmin()]
        r["kind"] = "divide"
        above = [h for _, h in HS if h > p.z_m]
        r["first_limit_above_pass"] = min(above) if above else np.nan
    else:
        e = lk[lk.barrier_id == b.barrier_id]
        if not len(e):
            rows.append(r | {"kind": "valley"})
            continue
        cnt = e["pass"].value_counts()
        p = deep.loc[cnt.index[0]]
        m = e["pass"] == cnt.index[0]
        r.update({"kind": "valley", "share_links": cnt.iloc[0] / len(e), "gap_km_to_pass": e.pass_km[m].median(),
                  "gap_z_m": e.gap_z_m[m].median()})
    r.update({"pass_id": p.apd_id, "pass_lat": p.lat, "pass_lon": p.lon, "pass_z_m": p.z_m, "pass_depth_m": p.depth_m,
              "pass_rank": p["rank"], "pass_unit": p.gmba, "control_pass": control(p)})
    rows.append(r)
out = pd.DataFrame(rows)
out.round(3).to_csv(f"{SUP}/S26_barrier_passes.csv", index=False, encoding="utf-8")

design_max = s24.set_index("no").loc[out.no, [f"n_{c}" for c in DESIGN]].max(axis=1).values
t3 = out[(out.defined_by == "design").values & (design_max >= NMIN)]
cols = ["no", "lat", "lon", "nearest_hazzi", "hazzi_km"] + [f"n_{h}" for _, h in HS] +        ["lowest_limit_nmin", "kind", "continental_inside", "pass_id", "pass_lat", "pass_lon", "pass_z_m",
        "pass_depth_m", "pass_rank", "pass_unit", "control_pass", "first_limit_above_pass", "gap_z_m"]
t3[cols].round(3).to_csv("outputs/tables/Table3_barrier_passes.csv", index=False, encoding="utf-8")
pd.set_option("display.width", 250)
print(t3[["no", "lat", "nearest_hazzi", "hazzi_km", "n_2000", "lowest_limit_nmin", "kind", "continental_inside",
          "pass_z_m", "pass_rank", "control_pass", "first_limit_above_pass", "gap_z_m"]].round(1).to_string(index=False))
v = t3[t3.kind == "divide"]
print("divide barriers whose support first reaches", NMIN, "at the first limit above their pass:",
      int((v.lowest_limit_nmin == v.first_limit_above_pass).sum()), "of", len(v))
