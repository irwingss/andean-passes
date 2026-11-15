"""57_barreras_pasos.py - The bridge of each barrier of 56_barreras_sintesis.py: the pass of the census through which
montane habitat comes closest to joining across the barrier.
For every barrier link, 58_puentes_links.py finds the highest route over the relief across the barrier between its two
patches (maximin path on the 1 km elevation) and the saddle at its lowest point. The census pass of a link is the pass
with a depth of 100 m or more within SNAP km of that saddle whose height is closest to it (the nearest pass if none is
that deep). The bridge of a
barrier at a montane limit is the pass that most of its links at that limit point to (ties go to the higher pass): the
bridge changes with the limit, because the patches that the barrier separates change.
Table 3 gives, for every barrier of the design that seven or more variables (of 19) support in some configuration of the
design, its support at each montane limit, its kind (divide, it contains continental passes, with the lowest of them;
valley, it contains none) and its bridge at 2000 m, or at the lowest limit with links if it has none at 2000 m, or in
the configuration with most links if it has links only with other margins or patch sizes.
Writes outputs/supp/tables/S26_barrier_passes.csv (one row per barrier of Table S24 and montane limit with links) and
outputs/tables/Table3_barrier_passes.csv.
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
SNAP = 5                                                    # km from the 1 km saddle to its census pass
CTL_KM = 3                                                  # a control pass within this distance names the pass
with rasterio.open(f"{B}/barreras_id.tif") as s:
    U, TR, CRS = s.read(1), s.transform, s.crs
to_m = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)


def xy(lon, lat):
    return np.column_stack(to_m.transform(np.asarray(lon), np.asarray(lat)))


apd = pd.read_csv("outputs/database/APD_passes.csv.gz",
                  usecols=["apd_id", "lat", "lon", "z_m", "depth_m", "rank", "gmba", "tropics"])
apd = apd[apd.tropics == 1].reset_index(drop=True)
r, c = rasterio.transform.rowcol(TR, *to_m.transform(apd.lon.values, apd.lat.values))
apd["barrier_id"] = U[np.clip(r, 0, U.shape[0] - 1), np.clip(c, 0, U.shape[1] - 1)]
tree = cKDTree(xy(apd.lon, apd.lat))
ctl = pd.read_csv("data/refs/pasos_control.csv", encoding="utf-8")
ctree = cKDTree(xy(ctl.lon, ctl.lat))

pl = pd.read_csv(f"{B}/puentes_links.csv")
P = xy(pl.bridge_lon, pl.bridge_lat)
snap = []
for k, near in enumerate(tree.query_ball_point(P, SNAP * 1000)):
    near = np.asarray(near, int)
    deep = near[apd.depth_m.values[near] >= 100]
    if deep.size:
        snap.append(deep[np.abs(apd.z_m.values[deep] - pl.bridge_z1km.iat[k]).argmin()])
    else:
        snap.append(tree.query(P[k])[1])
pl["pass"] = snap

s24 = pd.read_csv(f"{SUP}/S24_barriers.csv")
HZ = dict(HS)
pl["h"] = pl.cfg.map(HZ)


def control(p):
    cd, ci = ctree.query(to_m.transform(p.lon, p.lat))
    return ctl.name.iat[ci] if cd / 1000 <= CTL_KM else ""


def bridge(e):
    cnt = e["pass"].value_counts()
    k = max(cnt[cnt == cnt.iloc[0]].index, key=lambda i: apd.z_m.iat[i])
    p = apd.loc[k]
    return {"bridge_id": p.apd_id, "bridge_lat": p.lat, "bridge_lon": p.lon, "bridge_z_m": p.z_m,
            "bridge_z1km": e.bridge_z1km[e["pass"] == k].median(), "bridge_depth_m": p.depth_m,
            "bridge_rank": p["rank"], "bridge_unit": p.gmba, "bridge_name": control(p),
            "bridge_links": int(cnt.iloc[0]), "links": len(e)}


long, rows = [], []
for _, b in s24.iterrows():
    base = {"no": b.no, "barrier_id": b.barrier_id}
    e = pl[(pl.barrier_id == b.barrier_id) & pl.h.notna()]
    br = {}
    for h, g in e.groupby("h"):
        br[int(h)] = bridge(g)
        long.append(base | {"limit_m": int(h), "n_vars": b[f"n_{[c for c, x in HS if x == h][0]}"]} | br[int(h)])
    r = base | {"defined_by": b.defined_by, "lat": b.lat, "lon": b.lon, "nearest_hazzi": b.nearest_hazzi,
                "hazzi_km": b.hazzi_km}
    for c, h in HS:
        r[f"n_{h}"] = b[f"n_{c}"]
    on = [h for c, h in HS if b[f"n_{c}"] >= NMIN]
    r["lowest_limit_nmin"] = min(on) if on else np.nan
    inside = apd[apd.barrier_id == b.barrier_id]
    cont = inside[inside["rank"] == "continental"]
    r["continental_inside"] = len(cont)
    r["kind"] = "divide" if len(cont) else "valley"
    if len(cont):
        d = cont.loc[cont.z_m.idxmin()]
        r.update({"divide_pass_id": d.apd_id, "divide_pass_z_m": d.z_m, "divide_pass_name": control(d)})
    if br:
        h = 2000 if 2000 in br else min(br)
        r.update({"bridge_cfg": [c for c, x in HS if x == h][0]} | br[h])
    else:                                       # links only with other margins or patch sizes, at 2000 m
        e = pl[pl.barrier_id == b.barrier_id]
        if len(e):
            c = e.cfg.value_counts().index[0]
            r.update({"bridge_cfg": c} | bridge(e[e.cfg == c]))
    rows.append(r)
pd.DataFrame(long).round(3).to_csv(f"{SUP}/S26_barrier_passes.csv", index=False, encoding="utf-8")
out = pd.DataFrame(rows)
design_max = s24.set_index("no").loc[out.no, [f"n_{c}" for c in DESIGN]].max(axis=1).values
t3 = out[(out.defined_by == "design").values & (design_max >= NMIN)].drop(columns=["defined_by", "barrier_id"])
t3.round(3).to_csv("outputs/tables/Table3_barrier_passes.csv", index=False, encoding="utf-8")
pd.set_option("display.width", 250)
print(t3[["no", "lat", "nearest_hazzi", "hazzi_km", "n_2000", "kind", "divide_pass_z_m", "divide_pass_name",
          "bridge_cfg", "bridge_id", "bridge_z_m", "bridge_z1km", "bridge_rank", "bridge_name", "bridge_links",
          "links"]].round(1).to_string(index=False))
lg = pd.DataFrame(long)
print("bridges below their montane limit:", int((lg.bridge_z_m < lg.limit_m).sum()), "of", len(lg),
      "| within 500 m below it:", int(((lg.bridge_z_m < lg.limit_m) & (lg.bridge_z_m >= lg.limit_m - 500)).sum()))
