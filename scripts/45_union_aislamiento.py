"""45_union_aislamiento.py - Joins the core passes of every block and recomputes isolation and gateways over the whole
area (PRD O3). Isolation: ellipsoidal distance to the nearest lower pass with depth >= 100 m; gateway: depth >= 100 m and
isolation >= 100 km. The window values from 41_etiquetas_passes.py are kept only as a check (they are upper bounds).
python -P scripts/45_union_aislamiento.py      after 40 and 41 for all blocks
Output: data/interim/censo/passes_andes.csv
"""
import glob
import os
import re

import numpy as np
import pandas as pd
import yaml
from pyproj import Geod
from scipy.spatial import cKDTree

cfg = yaml.safe_load(open("config.yml", encoding="utf-8"))
OUT = os.path.join(cfg["paths"]["interim"], "censo")
GEOD = Geod(ellps="WGS84")
DEEP, GATE_KM = 100.0, 100.0

files = sorted(f for f in glob.glob(os.path.join(OUT, "passes_*_etiquetas.csv"))
               if re.search(r"passes_[NS]\d\dp00_etiquetas", f))          # integer blocks only; S05p50 was the overlap test
df = pd.concat([pd.read_csv(f, dtype={"hybas_term_a": str, "hybas_term_b": str}).assign(block=os.path.basename(f)[5:11])
                for f in files], ignore_index=True)
df = df[df["core"] == 1].reset_index(drop=True)            # cores tile the latitude without overlap: no duplicates
print(f"{len(files)} blocks, {len(df)} core passes")


def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


lat, lon, z = df["lat"].to_numpy(), df["lon"].to_numpy(), df["z_m"].to_numpy()
deep = np.flatnonzero(df["depth_m"].to_numpy() >= DEEP)
tree = cKDTree(xyz(lat[deep], lon[deep]))
P = xyz(lat, lon)
iso = np.full(len(df), np.nan)
todo = np.arange(len(df))
k = 16
# candidates by chord on the sphere, the winner by ellipsoidal distance among the lower ones. A pass is settled only when
# the k-th candidate lies beyond 1.01 x the best lower distance, so no closer lower pass can sit outside the k measured
# (sphere and ellipsoid distances differ by less than 0.5 %).
R = 6371.0088
while todo.size:
    kk = min(k, deep.size)
    dch, nn = tree.query(P[todo], k=kk)
    dch, nn = np.atleast_2d(dch), np.atleast_2d(nn)
    low = z[deep[nn]] < z[todo, None]
    for r in range(todo.size):
        q = todo[r]
        c = deep[nn[r][low[r]]]
        if c.size:
            _, _, d = GEOD.inv(np.full(c.size, lon[q]), np.full(c.size, lat[q]), lon[c], lat[c])
            best = d.min() / 1000
            if kk == deep.size or 2 * R * np.arcsin(min(dch[r, -1] / 2, 1)) > 1.01 * best:
                iso[q] = best
        elif kk == deep.size:
            iso[q] = np.inf                                 # no lower deep pass in the area
    todo = todo[np.isnan(iso[todo])]
    k *= 4
    print(f"  k={kk}: {todo.size} passes not settled", flush=True)

df["isolation_km"] = iso
df["gateway"] = ((df["depth_m"] >= DEEP) & (iso >= GATE_KM)).astype(int)
for d_, w_ in [(50, 100), (200, 100), (300, 100), (100, 50), (100, 250)]:
    df[f"gateway_d{d_}_w{w_}"] = ((df["depth_m"] >= d_) & (iso >= w_)).astype(int)   # ponytail: targets stay at 100 m depth

# a window can be shorter: its margin holds passes near the window edge that the neighbouring core does not confirm
bad = df["isolation_km"] > df["isolation_km_window"] + 0.006       # window values are rounded to 0.01 km
print(f"joined isolation longer than the window value in {int(bad.sum())} passes")
out = os.path.join(OUT, "passes_andes.csv")
df.to_csv(out, index=False)
g = df[df["gateway"] == 1]
print(f"gateways {len(g)} (window values gave {int(df['gateway_window'].sum())}); "
      f"by rank {g['rank'].value_counts().to_dict()}; below 2000 m {int((g['z_m'] < 2000).sum())} -> {out}")
