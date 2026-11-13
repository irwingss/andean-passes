"""exp_fig3_candidata.py - EXPERIMENT, a candidate Fig. 3 for the article: the barrier as a property of the zone.
(a, b) Zone conductance with BIO1, the mean relative conductance of the tiles that cover each cell, per reading. Low
is a barrier. The six barriers of Hazzi et al. (2018) are marked here and nowhere else, because this is the only
panel that measures the same thing they do: how hard it is to cross a zone.
(c, d) The same quantity along the continental divide, by group of variables, from the tiles that cross the divide.
Reads outputs/experimento_barreras_montanas (exp_barreras_montanas.py and exp_teselas_divisoria.py).
python -P scripts/exp_fig3_candidata.py principal
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyogrio
import rasterio
import shapely
from pyproj import Transformer

RUN = sys.argv[1]
D, OUT = f"data/interim/conectividad/{RUN}", "outputs/experimento_barreras_montanas"
sys.path.insert(0, "scripts")
from exp_estilo import BARRIER, GROUP, HALO, LAT_TICKS, lat_label, panel, save, setup       # noqa: E402

setup()
LN, LS = 11.0, -23.5
# the same coordinates as HAZZI in 52_figuras.py, digitized from their Fig. 1
HAZZI = {"I": [(-79.73, -7.02), (-77.80, -4.70)], "II": [(-74.14, -11.94), (-73.50, -13.49), (-72.81, -13.63)],
         "III": [(-68.88, -15.44), (-67.74, -14.21)], "IV": [(-77.11, 1.02), (-76.12, 2.18)],
         "V": [(-75.33, 3.20), (-74.04, 2.56)], "VI": [(-72.96, 7.99), (-71.82, 6.79)]}
with rasterio.open(f"{D}/insumos/elev.tif") as s:
    CRS, TR, H, W = s.crs, s.transform, s.height, s.width
to_m = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
gmba = shapely.from_wkb(pyogrio.raw.read("data/interim/M1_gmba_andes.gpkg")[2][0])
band_m = shapely.transform(shapely.intersection(gmba, shapely.box(-90, LS, -55, LN)),
                           lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1])))
bx0, by0, bx1, by1 = band_m.bounds
EXT = (bx0 - 120e3, bx1 + 120e3, by0 - 60e3, by1 + 60e3)
RAS = (TR.c, TR.c + W * TR.a, TR.f + H * TR.e, TR.f)
PAIS = [shapely.transform(shapely.from_wkb(w), lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1])))
        for w in pyogrio.raw.read("data/interim/M1_paises.gpkg")[2]]
db = pd.read_csv("outputs/database/APD_passes.csv.gz", dtype={"trans_robust": str, "longi_robust": str})
cont = db[(db.tropics == 1) & (db["rank"] == "continental")]
cx, cy = to_m.transform(cont.lon.values, cont.lat.values)
pr = pd.read_csv(f"{OUT}/perfil_divisoria.csv")


def frame(ax, labels=True):
    for p in PAIS:
        for g in getattr(p, "geoms", [p]):
            ax.plot(*g.exterior.xy, lw=0.25, c="0.45", alpha=0.85, zorder=4.5)
    for g in getattr(band_m, "geoms", [band_m]):
        ax.plot(*g.exterior.xy, lw=0.35, c="k", alpha=0.6, zorder=3)
    for la in LAT_TICKS:
        lo = np.linspace(-85, -55, 200)
        x, y = to_m.transform(lo, np.full_like(lo, la))
        ax.plot(x, y, lw=0.25, c="0.6", ls=":", zorder=0.5)
        if labels:
            ax.text(EXT[0] + 15e3, np.interp(EXT[0] + 15e3, x, y), lat_label(la), fontsize=5.5, color="0.35",
                    va="bottom", path_effects=HALO)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3]); ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")


fig = plt.figure(figsize=(183 / 25.4, 128 / 25.4))
mw = 0.265
for i, (lect, t) in enumerate((("alta", "Montane climates, longitudinal"), ("baja", "Lowland climates, transverse"))):
    ax = fig.add_axes([i * (mw + 0.015), 0.12, mw, 0.82])
    with rasterio.open(f"{OUT}/zona_{lect}_bio01.tif") as s:
        z = s.read(1)
    frame(ax, labels=i == 0)
    im = ax.imshow(z[::2, ::2], extent=RAS, cmap=BARRIER, vmin=0, vmax=1, interpolation="nearest", rasterized=True,
                   zorder=1)
    ax.scatter(cx, cy, s=0.1, c="w", lw=0, rasterized=True, zorder=2)        # the continental divide
    for k, pts in HAZZI.items():                                             # Hazzi et al. (2018), only here
        x, y = to_m.transform(*zip(*pts))
        ax.plot(x, y, c="w", lw=1.8, solid_capstyle="butt", zorder=6)
        ax.plot(x, y, c="k", lw=1.8, ls=(0, (1.2, 1.2)), dash_capstyle="butt", zorder=6.1)
        ax.annotate(k, (x[-1], y[-1]), xytext=(2, 1), textcoords="offset points", fontsize=5.5, fontweight="bold",
                    zorder=7, path_effects=HALO)
    panel(ax, "ab"[i], t)
cb = fig.colorbar(im, cax=fig.add_axes([0.10, 0.055, 0.32, 0.018]), orientation="horizontal")
cb.set_label("Zone conductance with BIO1: mean over the tiles that cover the cell")
cb.ax.annotate("barrier", (0, 0.5), xycoords="axes fraction", xytext=(-4, 0), textcoords="offset points", ha="right",
               va="center", fontsize=5.5)
cb.ax.annotate("open", (1, 0.5), xycoords="axes fraction", xytext=(4, 0), textcoords="offset points", ha="left",
               va="center", fontsize=5.5)
for i, rd in enumerate(("transverse", "longitudinal")):
    a = fig.add_axes([0.64 + i * 0.19, 0.12, 0.155, 0.82])
    q = pr[(pr.reading == rd) & (pr.cfg == "main")].set_index("band").sort_index()
    for g in ("temperature", "precipitation"):
        a.fill_betweenx(q.index, q[f"{g}_lo"], q[f"{g}_hi"], color=GROUP[g], lw=0, alpha=0.22, zorder=1)
        a.plot(q[g], q.index, c=GROUP[g], lw=1.2, label=g, zorder=3)
    a.axvline(1, c="0.5", lw=0.5, ls=":")
    a.set_xlim(0, 1.05); a.set_ylim(LS, LN); a.set_yticks(LAT_TICKS)
    a.set_yticklabels([lat_label(v) for v in LAT_TICKS] if i == 0 else [])
    a.set_xlabel("Relative conductance")
    a.spines[["top", "right"]].set_visible(False)
    panel(a, "cd"[i], f"{rd.capitalize()}")
    if i == 0:
        a.set_ylabel("Latitude of the tile centre")
        a.legend(fontsize=5.5, frameon=False, loc="lower left", borderaxespad=0.2, title="Tiles that cross the\n"
                 "divide; the band spans\nthe group's variables", title_fontsize=5, alignment="left")
save(fig, f"{OUT}/Fig3_candidata_conductancia")
