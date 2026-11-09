"""46_validacion_O5.py - Validation of the joined census (PRD O5): OSM pass nodes, identity against random points of the
divide network, heights against official references, control passes and the propagation of the height error.
python -P scripts/46_validacion_O5.py      after 45_union_aislamiento.py
Inputs: data/interim/censo/cols_andes.csv, divisorias_<tag>.tif, data/refs/. The OSM nodes are downloaded once by
latitude bands and kept with their date (ODbL). Output: data/interim/censo/validacion/.
"""
import datetime
import glob
import importlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
import pyogrio
import rasterio
import shapely
import yaml
from scipy.spatial import cKDTree

cfg = yaml.safe_load(open("config.yml", encoding="utf-8"))
RAW, INT = cfg["paths"]["data_dir"], cfg["paths"]["interim"]
CEN = os.path.join(INT, "censo")
OUT = os.path.join(CEN, "validacion")
os.makedirs(OUT, exist_ok=True)
R = 6371.0088
TROP_N, TROP_S = cfg["params"]["tropics_lat"]
ON_NET_CELLS = 3          # a node is on the network if a divide cell lies within 3 cells (about 90 m)
NREP = 1000
rng = np.random.default_rng(cfg["seed"])
sys.path.insert(0, "scripts")
AREA = importlib.import_module("40_censo_bloque").area_mask()

cols = pd.read_csv(os.path.join(CEN, "cols_andes.csv"))


def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


TREE = cKDTree(xyz(cols["lat"].to_numpy(), cols["lon"].to_numpy()))
DTREE_IDX = np.flatnonzero(cols["depth_m"].to_numpy() >= 100)
DTREE = cKDTree(xyz(cols["lat"].to_numpy()[DTREE_IDX], cols["lon"].to_numpy()[DTREE_IDX]))


def nearest(lat, lon, deep=False):
    """Index of the nearest col (or nearest col with depth >= 100 m) and its distance in km."""
    t = DTREE if deep else TREE
    d, i = t.query(xyz(np.atleast_1d(lat), np.atleast_1d(lon)))
    return (DTREE_IDX[i] if deep else i), 2 * R * np.arcsin(np.minimum(d / 2, 1))


def tag_of(lat):
    north = int(np.ceil(lat))                        # a core is (north - 1, north]
    return f"{'N' if north > 0 else 'S'}{abs(north):02d}p00"


# ---------------------------------------------------------------- 1. OSM nodes
osm_files = sorted(glob.glob("data/refs/osm_pasos_*.json"))
if osm_files:
    osm_file = osm_files[-1]
else:
    osm_file = f"data/refs/osm_pasos_{datetime.date.today():%Y%m%d}.json"
    els = {}
    for n in range(13, -57, -2):
        bb = f"{n - 2},-82,{n},-55"
        q = f'[out:json][timeout:180];(node["mountain_pass"="yes"]({bb});node["natural"="saddle"]({bb}););out;'
        for attempt in range(6):
            try:
                url = "https://overpass-api.de/api/interpreter?data=" + urllib.parse.quote(q)
                req = urllib.request.Request(url, headers={"User-Agent": "andean-passes-census/1.0"})
                with urllib.request.urlopen(req, timeout=240) as r:
                    for e in json.load(r)["elements"]:
                        els[e["id"]] = e
                break
            except Exception as ex:                     # Overpass often answers 502 or 504: wait and retry
                print(f"  band {n}: {ex}; retry", flush=True)
                time.sleep(30 * (attempt + 1))
        else:
            raise RuntimeError(f"band {n} failed")
        print(f"  band {n - 2} to {n}: {len(els)} nodes so far", flush=True)
    json.dump({"source": "OpenStreetMap contributors, ODbL, via overpass-api.de",
               "date": f"{datetime.date.today()}", "elements": list(els.values())},
              open(osm_file, "w", encoding="utf-8"), ensure_ascii=False)
E = json.load(open(osm_file, encoding="utf-8"))["elements"]
osm = pd.DataFrame({"osm_id": [e["id"] for e in E], "lat": [e["lat"] for e in E], "lon": [e["lon"] for e in E],
                    "name": [e.get("tags", {}).get("name", "") for e in E],
                    "pass_tag": [e.get("tags", {}).get("mountain_pass", "") == "yes" for e in E],
                    "saddle_tag": [e.get("tags", {}).get("natural", "") == "saddle" for e in E],
                    "ele_osm": pd.to_numeric([re.sub(r"[^0-9.\-]", "", str(e.get("tags", {}).get("ele", ""))) or None
                                              for e in E], errors="coerce")})
osm = osm[shapely.contains_xy(AREA, osm["lon"], osm["lat"])].reset_index(drop=True)
i, d = nearest(osm["lat"].to_numpy(), osm["lon"].to_numpy())
osm["col_km"], osm["col_z"], osm["col_depth"] = d, cols["z_m"].to_numpy()[i], cols["depth_m"].to_numpy()[i]
osm["gmba"], osm["rank"] = cols["gmba"].to_numpy()[i], cols["rank"].to_numpy()[i]
osm["col_lat"], osm["col_lon"] = cols["lat"].to_numpy()[i], cols["lon"].to_numpy()[i]
osm["tag"] = [tag_of(v) for v in osm["lat"]]

# on the network: any divide cell within ON_NET_CELLS of the node
osm["on_net"] = False
for tg, g in osm.groupby("tag"):
    with rasterio.open(os.path.join(CEN, f"divisorias_{tg}.tif")) as src:
        for k, r in g.iterrows():
            row, col = src.index(r["lon"], r["lat"])
            win = rasterio.windows.Window(col - ON_NET_CELLS, row - ON_NET_CELLS, 2 * ON_NET_CELLS + 1, 2 * ON_NET_CELLS + 1)
            osm.at[k, "on_net"] = bool(src.read(1, window=win, boundless=True, fill_value=0).any())
osm.to_csv(os.path.join(OUT, "O5_1_nodos_osm.csv"), index=False)


def shares(g):
    return pd.Series({"nodes": len(g), "within_500m": (g["col_km"] <= 0.5).mean(), "within_1km": (g["col_km"] <= 1).mean()})


t1 = pd.concat([shares(osm).to_frame("all").T,
                osm.groupby("gmba").apply(shares, include_groups=False),
                osm.groupby("rank").apply(shares, include_groups=False)])
t1.to_csv(os.path.join(OUT, "O5_1_proporciones.csv"))
osm[osm["col_km"] > 1].sort_values("col_km", ascending=False).to_csv(os.path.join(OUT, "O5_1_sin_col_1km.csv"), index=False)
print(f"1. {len(osm)} OSM nodes in the area ({osm_file}); col within 500 m {t1.loc['all', 'within_500m']:.3f}, "
      f"within 1 km {t1.loc['all', 'within_1km']:.3f}; on the network {int(osm['on_net'].sum())}", flush=True)

# ---------------------------------------------------------------- 2. identity against random points of the network
# Null: for each on-network node, NREP divide cells of the same block core and the same GMBA unit (else the same block).
gm_all = pyogrio.raw.read(os.path.join(RAW, "GMBA_Inventory_v2.0_standard_300/GMBA_Inventory_v2.0_standard_300.shp"),
                          columns=["Path"], encoding="utf-8")
path = np.array(gm_all[3][0], dtype=object)
sel = np.array(["> Andes" in str(p) for p in path])
gtree = shapely.STRtree(shapely.from_wkb(gm_all[2][sel]))
uname = np.array([str(p).split(" > ")[-1] for p in path[sel]])
net = osm[osm["on_net"]].reset_index(drop=True)
POOL = 20000
sim500 = np.zeros((NREP, len(net)), bool)
sim1k = np.zeros((NREP, len(net)), bool)
for tg, g in net.groupby("tag"):
    north = int(tg[1:3]) * (1 if tg[0] == "N" else -1)
    with rasterio.open(os.path.join(CEN, f"divisorias_{tg}.tif")) as src:
        top = int(round((src.bounds.top - north) / src.res[1]))
        nrow = int(round(1 / src.res[1]))
        a = src.read(1, window=rasterio.windows.Window(0, top, src.width, nrow))
        flat = np.flatnonzero(a)
        pick = flat[rng.choice(flat.size, min(flat.size, 4 * POOL), replace=False)]
        pr, pc = pick // a.shape[1] + top, pick % a.shape[1]
        plon, plat = rasterio.transform.xy(src.transform, pr, pc)
    plat, plon = np.asarray(plat), np.asarray(plon)
    ok = shapely.contains_xy(AREA, plon, plat)
    plat, plon = plat[ok][:POOL], plon[ok][:POOL]
    pg = uname[gtree.nearest(shapely.points(plon, plat))]
    _, pd_ = nearest(plat, plon)
    for k, r in g.iterrows():
        same = np.flatnonzero(pg == r["gmba"])
        src_i = same if same.size >= 50 else np.arange(plat.size)
        draw = rng.choice(src_i, NREP)
        sim500[:, k], sim1k[:, k] = pd_[draw] <= 0.5, pd_[draw] <= 1


def identity(mask, label):
    o500, o1k = (net.loc[mask, "col_km"] <= 0.5).mean(), (net.loc[mask, "col_km"] <= 1).mean()
    s500, s1k = sim500[:, mask].mean(1), sim1k[:, mask].mean(1)
    return {"group": label, "nodes": int(mask.sum()), "obs_500m": o500, "null_500m_mean": s500.mean(),
            "null_500m_q975": np.quantile(s500, 0.975), "p_500m": (np.sum(s500 >= o500) + 1) / (NREP + 1),
            "obs_1km": o1k, "null_1km_mean": s1k.mean(), "null_1km_q975": np.quantile(s1k, 0.975),
            "p_1km": (np.sum(s1k >= o1k) + 1) / (NREP + 1)}


rows = [identity(np.ones(len(net), bool), "all")]
rows += [identity((net["gmba"] == u).to_numpy(), u) for u in sorted(net["gmba"].unique())]
t2 = pd.DataFrame(rows)
t2.to_csv(os.path.join(OUT, "O5_2_identidad.csv"), index=False)
print("2. identity\n", t2.round(3).to_string(index=False), flush=True)

# ---------------------------------------------------------------- 3. heights against official references
ref = pd.read_csv("data/refs/M1_reference_passes.csv")


def dms(s, neg=True):
    d, m, sec = [float(x) for x in re.findall(r"[\d.]+", s)]
    return -(d + m / 60 + sec / 3600) if neg else d + m / 60 + sec / 3600


ref["lat"] = [r.lat if pd.notna(r.lat) else dms(r.lat_dms) for r in ref.itertuples()]    # every dms reference is in Peru: S, W
ref["lon"] = [r.lon if pd.notna(r.lon) else dms(r.lon_dms) for r in ref.itertuples()]
i, d = nearest(ref["lat"].to_numpy(), ref["lon"].to_numpy())
ref["col_km"], ref["col_z"], ref["col_depth"] = d, cols["z_m"].to_numpy()[i], cols["depth_m"].to_numpy()[i]
ref["error_m"] = ref["col_z"] - ref["elev_ref_m"]
ref["zone"] = np.where((ref["lat"] <= TROP_N) & (ref["lat"] >= TROP_S), "tropical", "temperate")
on = osm[osm["col_km"] <= 0.5].dropna(subset=["ele_osm"])
on = on.assign(error_m=on["col_z"] - on["ele_osm"],
               zone=np.where((on["lat"] <= TROP_N) & (on["lat"] >= TROP_S), "tropical", "temperate"))


def errs(g):
    e = g["error_m"]
    return pd.Series({"n": len(e), "median_error_m": e.median(), "median_abs_error_m": e.abs().median(),
                      "share_within_50m": (e.abs() <= 50).mean()})


t3 = pd.concat({"official": ref.groupby("zone").apply(errs, include_groups=False),
                "osm_ele (may come from a DEM)": on.groupby("zone").apply(errs, include_groups=False)})
ref.to_csv(os.path.join(OUT, "O5_3_alturas_referencia.csv"), index=False)
t3.to_csv(os.path.join(OUT, "O5_3_error_por_zona.csv"))
print("3. heights\n", t3.round(2).to_string(), flush=True)

# ---------------------------------------------------------------- 4. control passes, relief only
ctl = pd.read_csv("data/refs/pasos_control.csv")
ctl = ctl[ctl["in_area"] == 1].reset_index(drop=True)
i, d = nearest(ctl["lat"].to_numpy(), ctl["lon"].to_numpy())
j, dd = nearest(ctl["lat"].to_numpy(), ctl["lon"].to_numpy(), deep=True)
c = cols
t4 = pd.DataFrame({"control": ctl["name"], "system": ctl["system"], "lat": ctl["lat"], "lon": ctl["lon"],
                   "gmba": c["gmba"].to_numpy()[i], "near_km": d.round(2), "near_z": c["z_m"].to_numpy()[i],
                   "near_depth": c["depth_m"].to_numpy()[i], "near_rank": c["rank"].to_numpy()[i],
                   "near_area_sep_km2": c["area_sep_km2"].to_numpy()[i], "near_isolation_km": c["isolation_km"].to_numpy()[i].round(2),
                   "deep_km": dd.round(2), "deep_z": c["z_m"].to_numpy()[j], "deep_depth": c["depth_m"].to_numpy()[j],
                   "deep_rank": c["rank"].to_numpy()[j], "deep_area_sep_km2": c["area_sep_km2"].to_numpy()[j],
                   "deep_width_m": c["width_m"].to_numpy()[j], "deep_isolation_km": c["isolation_km"].to_numpy()[j].round(2)})
t4.to_csv(os.path.join(OUT, "O5_4_pasos_control.csv"), index=False)
print(f"4. control passes: {len(t4)}, nearest col at most {t4['near_km'].max():.2f} km; "
      f"{int((t4['near_km'] < 1).sum())} under 1 km", flush=True)

# ---------------------------------------------------------------- 5. propagation of the tropical median error
e = float(t3.loc[("official", "tropical"), "median_error_m"])
g = cols[cols["gateway"] == 1]
trop = (g["lat"] <= TROP_N) & (g["lat"] >= TROP_S)
t5 = pd.DataFrame([{"shift_m": s, "gateways": len(g), "below_2000m": int((g["z_m"] - s < 2000).sum()),
                    "tropical_below_2000m": int((g.loc[trop, "z_m"] - s < 2000).sum())} for s in (-e, 0.0, e)])
# a uniform shift moves every col alike: isolation, gateways and their order by height do not change
t5.to_csv(os.path.join(OUT, "O5_5_propagacion.csv"), index=False)
print(f"5. tropical median error {e:.1f} m\n", t5.to_string(index=False), flush=True)

# ---------------------------------------------------------------- 6. cases
p = ref[ref["name"] == "Abra de Porculla"].iloc[0]
cu = ref[ref["name"] == "Cumbre"].iloc[0]
with open(os.path.join(OUT, "O5_6_casos.txt"), "w", encoding="utf-8") as f:
    f.write(f"Porculla: col {p.col_z:.1f} m at {p.col_km:.2f} km of the reference point, depth {p.col_depth:.0f} m; "
            f"official {p.elev_ref_m:.0f} m, alternative {p.elev_alt_m:.0f} m.\n")
    f.write(f"Cumbre: col {cu.col_z:.1f} m at {cu.col_km:.2f} km, depth {cu.col_depth:.0f} m; reference {cu.elev_ref_m:.0f} m. "
            f"Needs a check on official cartography.\n")
print(open(os.path.join(OUT, "O5_6_casos.txt"), encoding="utf-8").read())
