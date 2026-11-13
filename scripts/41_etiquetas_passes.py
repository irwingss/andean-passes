"""41_etiquetas_passes.py - Labels for the passes of one block (PRD O1, O3, O4): GMBA level-4 unit, sea and terminal basin
of each side from HydroBASINS, rank, separated area, isolation and gateway. HydroBASINS only labels: it never adds,
moves or removes a pass.
python -P scripts/41_etiquetas_passes.py <tag>      e.g. S05p00, after 40_censo_bloque.py
Isolation here only sees the passes of this block's window, so it is an upper bound until the blocks are joined (R3).
"""
import os
import sys

import numpy as np
import pyogrio
import shapely
import yaml
from pyproj import Geod

cfg = yaml.safe_load(open("config.yml", encoding="utf-8"))
RAW, INT = cfg["paths"]["data_dir"], cfg["paths"]["interim"]
OUT = os.path.join(INT, "censo")
tag = sys.argv[1]
GEOD = Geod(ellps="WGS84")



passes = np.genfromtxt(os.path.join(OUT, f"passes_{tag}.csv"), delimiter=",", names=True, encoding="utf-8")
lat, lon = passes["lat"], passes["lon"]

# ---- GMBA level-4 unit (point in polygon, else the nearest unit: passes in the 20 km margin)
m, _, g, F = pyogrio.raw.read(os.path.join(RAW, "GMBA_Inventory_v2.0_standard_300/GMBA_Inventory_v2.0_standard_300.shp"),
                              columns=["Path"], encoding="utf-8")
path = np.array(F[0], dtype=object)
sel = np.array(["> Andes" in str(p) for p in path])
units = shapely.from_wkb(g[sel])
uname = np.array([str(p).split(" > ")[-1] for p in path[sel]])
tree = shapely.STRtree(units)
pts = shapely.points(lon, lat)
gi = tree.nearest(pts)            # a point inside a polygon is at distance 0 from it
gmba = uname[gi]

# ---- HydroBASINS level 12: basin of each side's 5 km2 point, path to the sea by NEXT_DOWN
f12 = os.path.join(RAW, "hybas_sa_lev01-12_v1c/hybas_sa_lev12_v1c.shp")
m, _, _, F = pyogrio.raw.read(f12, columns=["HYBAS_ID", "NEXT_DOWN", "MAIN_BAS", "UP_AREA", "ENDO"],
                              read_geometry=False)
H = dict(zip(m["fields"], F))
hid = H["HYBAS_ID"].astype(np.int64)
pos = {h: k for k, h in enumerate(hid)}
nxt, mainb, uparea, endo = H["NEXT_DOWN"].astype(np.int64), H["MAIN_BAS"].astype(np.int64), H["UP_AREA"], H["ENDO"]
pad = 0.5
bb = (np.nanmin([passes["lon_a5"].min(), passes["lon_b5"].min()]) - pad, np.nanmin([passes["lat_a5"].min(), passes["lat_b5"].min()]) - pad,
      np.nanmax([passes["lon_a5"].max(), passes["lon_b5"].max()]) + pad, np.nanmax([passes["lat_a5"].max(), passes["lat_b5"].max()]) + pad)
m, _, gb, Fb = pyogrio.raw.read(f12, columns=["HYBAS_ID"], bbox=bb)
bpoly, bid = shapely.from_wkb(gb), Fb[0].astype(np.int64)
btree = shapely.STRtree(bpoly)


def basin_of(lo, la):
    p = shapely.points(lo, la)
    qi = btree.nearest(p)
    return bid[qi]


# sea of each terminal basin, from its outlet position; endorheic basins are interior only if their endorheic region
# touches the Atlantic side (M1_domains: W, E, I dissolved from the same HydroBASINS); coastal sinks keep their coast.
m, _, gt, Ft = pyogrio.raw.read(f12, columns=["HYBAS_ID", "NEXT_DOWN"], where="NEXT_DOWN = 0")  # the filtered field must be read
tpt = dict(zip(Ft[0].astype(np.int64), shapely.point_on_surface(shapely.from_wkb(gt))))
m, _, gd, Fd = pyogrio.raw.read(os.path.join(INT, "M1_domains.gpkg"), columns=["dom"])
dom = dict(zip(Fd[0], shapely.from_wkb(gd)))
iparts = [p for p in shapely.get_parts(dom["I"]) if shapely.intersects(p, dom["E"])]
inter = shapely.union_all(iparts)


def sea(term):
    p = tpt[term]
    x, y = p.x, p.y
    if endo[pos[term]] == 2:
        return "interior" if shapely.intersects(inter, p) else "Pacific"
    if y > 7.5:
        return "Caribbean" if x < -62 else "Atlantic"
    if y >= -41:
        return "Atlantic" if x > -66 else "Pacific"
    return "Atlantic" if x > -69.8 else "Pacific"


def route(h):
    r = [h]
    while nxt[pos[r[-1]]] != 0:
        r.append(nxt[pos[r[-1]]])
    return r


ba, bbb = basin_of(passes["lon_a5"], passes["lat_a5"]), basin_of(passes["lon_b5"], passes["lat_b5"])
cache = {}
rank, sea_a, sea_b, area_sep, term_a, term_b = [], [], [], [], [], []
for a, b in zip(ba, bbb):
    ra, rb = cache.setdefault(a, route(a)), cache.setdefault(b, route(b))
    ta, tb = ra[-1], rb[-1]
    sa, sb = sea(ta), sea(tb)
    sea_a.append(sa); sea_b.append(sb); term_a.append(ta); term_b.append(tb)
    if sa != sb:                         # different seas, or one side interior
        rank.append("continental")
    elif ta != tb:
        rank.append("terminal")
    else:
        rank.append("tributary")
    if ta != tb:
        area_sep.append(max(uparea[pos[ta]], uparea[pos[tb]]))
    elif a == b:
        area_sep.append(np.nan)          # both sides in one level-12 basin: below HydroBASINS resolution
    else:
        sa_set = set(ra)
        j = next(x for x in rb if x in sa_set)               # first basin shared by both routes
        pre = [r[r.index(j) - 1] for r in (ra, rb) if r.index(j) > 0]
        area_sep.append(min(uparea[pos[x]] for x in pre))

# ---- isolation: distance to the nearest lower pass with depth >= 100 m (window passes only, upper bound)
z, depth = passes["z_m"], passes["depth_m"]
deep = np.flatnonzero(depth >= 100)
iso = np.full(len(passes), np.inf)
for k in range(len(passes)):
    cand = deep[z[deep] < z[k]]
    if cand.size:
        _, _, d = GEOD.inv(np.full(cand.size, lon[k]), np.full(cand.size, lat[k]), lon[cand], lat[cand])
        iso[k] = d.min() / 1000
gate = (depth >= 100) & (iso >= 100)

# ---- write
names = list(passes.dtype.names)
out = os.path.join(OUT, f"passes_{tag}_etiquetas.csv")
with open(out, "w", encoding="utf-8", newline="") as f:
    f.write(",".join(names + ["gmba", "sea_a", "sea_b", "rank", "hybas_term_a", "hybas_term_b", "area_sep_km2",
                              "isolation_km_window", "gateway_window"]) + "\n")
    for k in range(len(passes)):
        row = ["" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in passes[k].tolist()]
        row += [gmba[k], sea_a[k], sea_b[k], rank[k], str(term_a[k]), str(term_b[k]),
                "" if np.isnan(area_sep[k]) else f"{area_sep[k]:.1f}",
                "inf" if np.isinf(iso[k]) else f"{iso[k]:.2f}", str(int(gate[k]))]
        f.write(",".join(f'"{x}"' if "," in x else x for x in row) + "\n")
fd = {"z_m": z, "depth_m": depth, "width_m": passes["width_m"], "ring": passes["ring"].astype(int), "flat": passes["flat"].astype(int),
      "core": passes["core"].astype(int), "gmba": gmba.astype(str), "rank": np.array(rank, dtype=str),
      "area_sep_km2": np.array(area_sep, float), "isolation_km": iso, "gateway": gate.astype(int)}
pyogrio.raw.write(os.path.join(OUT, f"passes_{tag}.gpkg"), shapely.to_wkb(pts), list(fd.values()), list(fd.keys()),
                  geometry_type="Point", crs="EPSG:4326", driver="GPKG")
core = passes["core"] == 1
print(f"{tag}: {core.sum()} core passes; rank", {r: int(np.sum(np.array(rank)[core] == r)) for r in set(rank)},
      f"; gateways (window) {int(gate[core].sum())} -> {out}")
