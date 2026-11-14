"""56_barreras_sintesis.py - Support of the barriers of 55_barreras_grafo.py across the 19 variables and the thirteen
configurations (nine montane limits from 500 to 4500 m, two margins q, two patch areas), and Fig. 5.
A barrier line of one variable supports a cell when it passes within R km of it (lines of different variables trace the
same barrier a few km apart). A barrier is a connected set of cells supported by at least one variable in any
configuration of the design (main, q, patch area, 1500 and 2500 m), inside the tropical GMBA units; cells supported
only at the other montane limits form further barriers, one limit at a time (column defined_by); its support in a configuration is the number of variables (of 19, of
the 11 of temperature and of the 8 of precipitation) whose lines touch it.
Writes outputs/supp/tables/S24_barriers.csv (one row per barrier), S25_barrier_links.csv (every link that is a barrier
for some variable and configuration), outputs/figures/Fig5_barriers.{png,pdf} and
data/interim/conectividad/<run>/barreras/apoyo_<cfg>.tif (number of variables per cell), lineas_<cfg>.tif (that number on
the barrier lines) and barreras_id.tif
(barrier_id of every cell).
python -P scripts/56_barreras_sintesis.py principal
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyogrio
import rasterio
import shapely
from matplotlib.colors import BoundaryNorm
from pyproj import Transformer
from rasterio import features
from rasterio.windows import from_bounds
from scipy import ndimage as ndi

RUN = sys.argv[1]
D = f"data/interim/conectividad/{RUN}"
B = f"{D}/barreras"
SUP, SUPF = "outputs/supp/tables", "outputs/supp"
LN, LS = 11.0, -23.5
R = 10                                                     # km: lines of different variables within R are one barrier
HS_CFG = ["h500", "h1000", "h1500", "main", "h2500", "h3000", "h3500", "h4000", "h4500"]   # the nine montane limits
CFG = ["main", "q0.05", "q0.2", "a100", "a1000", "h500", "h1000", "h1500", "h2500", "h3000", "h3500", "h4000", "h4500"]
CORE = ["main", "q0.05", "q0.2", "a100", "a1000", "h1500", "h2500"]   # the design: main and its sensitivities (D9)
CFG_LAB = {"main": "2000 m, q 0.1, 500 km²", "q0.05": "q 0.05", "q0.2": "q 0.2", "a100": "patches of 100 km²",
           "a1000": "patches of 1000 km²", **{f"h{h}": f"montane limit {h} m" for h in (500, 1000, 1500, 2500, 3000, 3500, 4000, 4500)}}
VARS = [f"bio{v:02d}" for v in range(1, 20)]
TEMP, PREC = VARS[:11], VARS[11:]
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.5, "pdf.fonttype": 42,
                     "savefig.dpi": 600})
HALO = [matplotlib.patheffects.withStroke(linewidth=1.6, foreground="w")]

with rasterio.open(f"{D}/insumos/elev.tif") as s:
    PROF, TR, CRS, H, W = s.profile, s.transform, s.crs, s.height, s.width
    Z = s.read(1, masked=True).filled(np.nan)
to_m = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
to_g = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
gmba = shapely.from_wkb(pyogrio.raw.read("data/interim/M1_gmba_andes.gpkg")[2][0])
band = shapely.intersection(gmba, shapely.box(-90, LS, -55, LN))
band_m = shapely.transform(band, lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1])))
IN = features.rasterize([(band_m, 1)], out_shape=(H, W), transform=TR, dtype=np.uint8).astype(bool)


def lines(v, cfg):
    with rasterio.open(f"{B}/lineas_{v}_{cfg}.tif") as s:
        return (s.read(1) == 1) & IN


# per configuration: lines of every variable, and the cells within R km of each
LINE, NEAR = {}, {}
for cfg in CFG:
    for v in VARS:
        l = lines(v, cfg)
        LINE[(cfg, v)] = l
        NEAR[(cfg, v)] = ndi.distance_transform_edt(~l) <= R if l.any() else l
    cnt = sum(NEAR[(cfg, v)].astype(np.uint8) for v in VARS)
    p = PROF.copy(); p.update(dtype="uint8", nodata=0, compress="deflate", count=1)
    with rasterio.open(f"{B}/apoyo_{cfg}.tif", "w", **p) as o:
        o.write(cnt, 1)
    onl = np.zeros((H, W), bool)
    for v in VARS:
        onl |= LINE[(cfg, v)]
    with rasterio.open(f"{B}/lineas_{cfg}.tif", "w", **p) as o:          # support on the barrier lines (Fig. 5, viewer)
        o.write(np.where(onl, cnt, 0).astype(np.uint8), 1)
    print(cfg, "support written", flush=True)

# barriers: connected sets of cells near a line in any configuration of the design (CORE); the cells near lines of
# the other montane limits that stay outside them form further barriers, one montane limit at a time, so that the lines
# of many heights do not chain neighbouring barriers into one
U = np.zeros((H, W), np.int32)
ORIGIN = {}
nu = 0
for grp in [CORE] + [[c] for c in CFG if c not in CORE]:
    anyl = np.zeros((H, W), bool)
    for c in grp:
        for v in VARS:
            anyl |= NEAR[(c, v)]
    lab, n = ndi.label(anyl & (U == 0), structure=np.ones((3, 3)))
    U[lab > 0] = lab[lab > 0] + nu
    for k in range(nu + 1, nu + n + 1):
        ORIGIN[k] = "design" if grp is CORE else grp[0]
    nu += n
print(nu, "barriers")
p = PROF.copy(); p.update(dtype="int32", nodata=0, compress="deflate", count=1)
with rasterio.open(f"{B}/barreras_id.tif", "w", **p) as o:            # barrier_id of every cell, for 57_barreras_pasos.py
    o.write(U, 1)
links = pd.concat([pd.read_csv(f) for f in glob.glob(f"{B}/links_bio??_*.csv")])
links = links[links.barrier & (links.lat <= LN) & (links.lat >= LS)]
lr, lc = rasterio.transform.rowcol(TR, *to_m.transform(links.lon.values, links.lat.values))
links["barrier_id"] = U[np.clip(lr, 0, H - 1), np.clip(lc, 0, W - 1)]
links.to_csv(f"{SUP}/S25_barrier_links.csv", index=False, float_format="%.4g")
hz = pd.read_csv("data/refs/hazzi_2018_barreras.csv")
hz_geom = {n: shapely.LineString(np.column_stack(to_m.transform(g.lon.values, g.lat.values)))
           for n, g in hz.groupby("name")}
rows = []
idx = ndi.find_objects(U)
for k in range(1, nu + 1):
    sl = idx[k - 1]
    m = U[sl] == k
    rr, cc = np.nonzero(m)
    rr, cc = rr + sl[0].start, cc + sl[1].start
    x, y = TR.c + (cc + 0.5) * TR.a, TR.f + (rr + 0.5) * TR.e
    lo, la = to_g.transform(x.mean(), y.mean())
    r = {"barrier_id": k, "defined_by": ORIGIN[k], "lat": la, "lon": lo, "cells": int(m.sum())}
    pts = shapely.MultiPoint(np.column_stack([x[::25], y[::25]]))
    dist = {n: pts.distance(g) / 1000 for n, g in hz_geom.items()}
    near = min(dist, key=dist.get)
    r["nearest_hazzi"], r["hazzi_km"] = near, round(dist[near], 1)
    for cfg in CFG:
        sup = [v for v in VARS if (LINE[(cfg, v)][sl] & m).any()]
        r[f"n_{cfg}"] = len(sup)
        r[f"temp_{cfg}"] = sum(v in TEMP for v in sup)
        r[f"prec_{cfg}"] = sum(v in PREC for v in sup)
        if cfg == "main":
            r["vars_main"] = " ".join(v.replace("bio", "") for v in sup)
    e = links[(links.barrier_id == k) & (links.cfg == "main")]
    r["median_excess_main"] = e.excess.median() if len(e) else np.nan
    r["median_z_m"] = float(np.nanmedian(Z[rr, cc]))
    rows.append(r)
s24 = pd.DataFrame(rows)
s24["n_max"] = s24[[f"n_{c}" for c in CFG]].max(axis=1)
s24 = s24[s24.n_max > 0].sort_values("lat", ascending=False)
s24["no"] = np.arange(1, len(s24) + 1)
s24 = s24[["no", "barrier_id", "defined_by", "lat", "lon", "cells", "median_z_m", "nearest_hazzi", "hazzi_km", "vars_main",
           "median_excess_main", "n_max"] + [f"{p}_{c}" for c in CFG for p in ("n", "temp", "prec")]]
s24.round(3).to_csv(f"{SUP}/S24_barriers.csv", index=False)
print(s24[["no", "lat", "lon", "nearest_hazzi", "hazzi_km", "n_h500", "n_h1000", "n_h1500", "n_main", "n_h2500", "n_h3000", "n_h3500",
          "n_h4000", "n_h4500", "n_max"]].round(2)
      .to_string(index=False))

# ==== Fig. 5 =====================================================================================================
EXT0 = band_m.bounds
EXT = (EXT0[0] - 120e3, EXT0[2] + 120e3, EXT0[1] - 60e3, EXT0[3] + 60e3)
win = from_bounds(EXT[0], EXT[2], EXT[1], EXT[3], TR).round_offsets().round_lengths()
e = rasterio.windows.bounds(win, TR)
EXT = (e[0], e[2], e[1], e[3])
crop = lambda a: a[win.row_off:win.row_off + win.height, win.col_off:win.col_off + win.width]
gy, gx = np.gradient(np.nan_to_num(crop(Z)) * 8 / 1000.0)
slp, asp = np.arctan(np.hypot(gx, gy)), np.arctan2(-gx, gy)
HS = np.cos(np.radians(50)) * np.cos(slp) + np.sin(np.radians(50)) * np.sin(slp) * np.cos(np.radians(225) - asp)
HS = np.ma.masked_where(np.isnan(crop(Z)), HS)
PAIS = [shapely.from_wkb(w) for w in pyogrio.raw.read("data/interim/M1_paises.gpkg")[2]]
PAIS_M = [shapely.transform(p, lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1]))) for p in PAIS]
CMAP = plt.get_cmap("plasma_r", 19)
NORM = BoundaryNorm(np.arange(0.5, 20.5, 1), 19)


ZC = crop(Z)
GREEN = matplotlib.colors.ListedColormap(["#2e9e4f"])


def draw(ax, cfg, title, s=0.6):
    ax.imshow(HS, extent=EXT, cmap="Greys_r", vmin=-0.4, vmax=1.4, interpolation="bilinear", rasterized=True)
    h = 2000 if cfg == "main" else int(cfg[1:])
    high = np.ma.masked_where(~(ZC >= h), np.ones_like(ZC))        # land at or above the montane limit
    ax.imshow(high, extent=EXT, cmap=GREEN, alpha=0.4, interpolation="nearest", zorder=1.5, rasterized=True)
    for p in PAIS_M:
        for g in getattr(p, "geoms", [p]):
            ax.plot(*g.exterior.xy, lw=0.25, c="0.45", alpha=0.85, zorder=2)
    for g in getattr(band_m, "geoms", [band_m]):
        ax.plot(*g.exterior.xy, lw=0.35, c="k", alpha=0.5, zorder=2)
    cnt = sum(NEAR[(cfg, v)].astype(np.uint8) for v in VARS)
    anyline = np.zeros((H, W), bool)
    for v in VARS:
        anyline |= LINE[(cfg, v)]
    rr, cc = np.nonzero(crop(anyline))
    val = crop(cnt)[rr, cc]
    o = np.argsort(val)                                     # the best supported on top
    x = EXT[0] + (cc[o] + 0.5) * 1000; y = EXT[3] - (rr[o] + 0.5) * 1000
    sc = ax.scatter(x, y, c=val[o], cmap=CMAP, norm=NORM, s=s, lw=0, zorder=3, rasterized=True)
    for n, g in hz_geom.items():
        ax.plot(*g.xy, c="w", lw=2.2, zorder=4)
        ax.plot(*g.xy, c="k", lw=1.1, ls=(0, (1.5, 1.5)), zorder=4.1)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    ax.annotate(title, (0, 1), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points", fontsize=7,
                va="bottom", ha="left")
    return sc


Wd = 183.0
fig = plt.figure(figsize=(Wd / 25.4, 300 / 25.4))
gs = fig.add_gridspec(3, 3, hspace=0.06, wspace=0.03, bottom=0.06)
for j, cfg in enumerate(HS_CFG):
    h = 2000 if cfg == "main" else int(cfg[1:])
    t = f"{'abcdefghi'[j]}  Montane limit {h} m" + (" (main)" if cfg == "main" else "")
    sc = draw(fig.add_subplot(gs[j // 3, j % 3]), cfg, t)
cax = fig.add_axes([0.25, 0.035, 0.5, 0.008])
cb = fig.colorbar(sc, cax=cax, orientation="horizontal", ticks=[1, 4, 7, 10, 13, 16, 19])
cb.set_label(f"Variables (of 19) with a barrier line within {R} km; green, land at or above the montane limit;\ndashed, Hazzi et al. (2018)", fontsize=6.5)
fig.savefig(f"outputs/figures/Fig5_barriers.png", bbox_inches="tight", pad_inches=0.03)
fig.savefig(f"outputs/figures/Fig5_barriers.pdf", bbox_inches="tight", pad_inches=0.03)
print("->", f"outputs/figures/Fig5_barriers")
