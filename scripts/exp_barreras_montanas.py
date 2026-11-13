"""exp_barreras_montanas.py - EXPERIMENT, not part of the article: the montane barrier as a property of the zone.
PRD C7 reads a pass with two numbers: the relative conductance of the tiles that cover it (strength of the bridge or
barrier of its zone) and its normalized current (its role inside the zone, channel or pinch point). The census holds
only the second. This script computes the first from the tile tables that already exist:
  zone conductance of a 1 km cell = mean relative conductance of the tiles that cover it, per reading and configuration
and writes, in outputs/experimento_barreras_montanas:
  zona_{alta,baja}_bio01.tif       zone conductance with BIO1, main configuration (and the minimum and maximum over the
                                   five configurations as bands 2 and 3)
  passes_zona.csv.gz                 every tropical pass with the zone conductance of both readings
python -P scripts/exp_barreras_montanas.py principal
"""
import os
import sys

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer

RUN = sys.argv[1]
D = f"data/interim/conectividad/{RUN}"
OUT = "outputs/experimento_barreras_montanas"
os.makedirs(OUT, exist_ok=True)
CFG = {"main": ("s200", ""), "q005": ("s200", "_q0.05"), "q02": ("s200", "_q0.2"), "100km": ("s100", ""),
       "400km": ("s400", "")}
with rasterio.open(f"{D}/insumos/elev.tif") as s:
    prof, H, W, TR, CRS = s.profile, s.height, s.width, s.transform, s.crs


def zone(var, lect, cfg):
    """Mean relative conductance of the tiles of one reading that cover each cell; nan where no tile does."""
    side, q = CFG[cfg]
    S = int(side[1:])
    t = pd.read_csv(f"{D}/teselas/teselas_{side}_{var}{q}.csv")
    t = t[(t.lectura == lect) & (t.resistencia.str.endswith("env") | t.resistencia.str.contains("_env_q"))]
    su, n = np.zeros((H, W)), np.zeros((H, W))
    for r, c, g in zip(t.row, t.col, t.conductancia_rel):
        a, b = max(r - S // 2, 0), max(c - S // 2, 0)
        su[a:r - S // 2 + S, b:c - S // 2 + S] += g
        n[a:r - S // 2 + S, b:c - S // 2 + S] += 1
    return np.where(n > 0, su / np.maximum(n, 1), np.nan)


db = pd.read_csv("outputs/database/APD_passes.csv.gz", dtype={"trans_robust": str, "longi_robust": str})
tr = db[db.tropics == 1].copy()
x, y = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(tr.lon.values, tr.lat.values)
rr, cc = rasterio.transform.rowcol(TR, x, y)
rr, cc = np.clip(rr, 0, H - 1), np.clip(cc, 0, W - 1)
prof.update(count=3, dtype="float32", nodata=np.nan, compress="deflate")
for lect, p in (("alta", "longi"), ("baja", "trans")):
    z = np.stack([zone("bio01", lect, c) for c in CFG])
    with np.errstate(all="ignore"):
        lo, hi = np.nanmin(z, axis=0), np.nanmax(z, axis=0)
    lo[np.isnan(z).any(axis=0)] = np.nan; hi[np.isnan(z).any(axis=0)] = np.nan     # a value in all five, as in C8
    with rasterio.open(f"{OUT}/zona_{lect}_bio01.tif", "w", **prof) as o:
        o.write(np.stack([z[0], lo, hi]).astype("float32"))
    tr[f"{p}_zone"], tr[f"{p}_zone_lo"], tr[f"{p}_zone_hi"] = z[0][rr, cc], lo[rr, cc], hi[rr, cc]
    print(lect, "cells with a value", np.isfinite(z[0]).sum(), "in all five", np.isfinite(lo).sum())
keep = ["apd_id", "lat", "lon", "z_m", "depth_m", "rank", "gmba", "cell_1km", "trans_med", "trans_robust", "trans_zone",
        "trans_zone_lo", "trans_zone_hi", "longi_med", "longi_robust", "longi_zone", "longi_zone_lo", "longi_zone_hi"]
tr[keep].round(4).to_csv(f"{OUT}/passes_zona.csv.gz", index=False)
