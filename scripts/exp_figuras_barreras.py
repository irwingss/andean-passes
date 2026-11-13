"""exp_figuras_barreras.py - EXPERIMENT, not part of the article: figures from exp_barreras_montanas.py.
  FigA_puente_y_barrera   passes that concentrate lowland current and divert montane current in all five configurations
  FigB_conductancia_zona  zone conductance with BIO1 in both readings, and its profile along the continental divide
python -P scripts/exp_figuras_barreras.py principal
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import pyogrio
import rasterio
import shapely
from pyproj import Transformer

RUN = sys.argv[1]
D, OUT = f"data/interim/conectividad/{RUN}", "outputs/experimento_barreras_montanas"
sys.path.insert(0, "scripts")
from exp_estilo import BARRIER, BOTH, CONC, DIV, HALO, LAT_TICKS, lat_label, panel, save, setup  # noqa: E402

setup()
LN, LS = 11.0, -23.5
with rasterio.open(f"{D}/insumos/elev.tif") as s:
    CRS, TR, H, W = s.crs, s.transform, s.height, s.width
to_m = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
gmba = shapely.from_wkb(pyogrio.raw.read("data/interim/M1_gmba_andes.gpkg")[2][0])
band = shapely.intersection(gmba, shapely.box(-90, LS, -55, LN))
band_m = shapely.transform(band, lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1])))
bx0, by0, bx1, by1 = band_m.bounds
EXT = (bx0 - 120e3, bx1 + 120e3, by0 - 60e3, by1 + 60e3)
RAS = (TR.c, TR.c + W * TR.a, TR.f + H * TR.e, TR.f)                 # extent of the 1 km grid
# country outlines, for orientation only: GADM 4.1 dissolved to country (scripts/exp_paises.R). Only the six Andean
# countries are in data/raw, so Argentina and Brazil show as the far side of a border, without an outline of their own.
PAIS = [shapely.transform(shapely.from_wkb(w), lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1])))
        for w in pyogrio.raw.read("data/interim/M1_paises.gpkg")[2]]


def frame(ax, labels=True, fill=True):
    for p in PAIS:
        for g in getattr(p, "geoms", [p]):
            ax.plot(*g.exterior.xy, lw=0.25, c="0.45", alpha=0.85, zorder=4.5)
    for g in getattr(band_m, "geoms", [band_m]):
        if fill:
            ax.fill(*g.exterior.xy, fc="0.95", ec="none", zorder=0)
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


t = pd.read_csv(f"{OUT}/passes_zona.csv.gz")
u = t.sort_values("depth_m", ascending=False).drop_duplicates("cell_1km")       # one pass per 1 km cell, the deepest
u = u[u.gmba.notna()]
ux, uy = to_m.transform(u.lon.values, u.lat.values)

# ==== Fig. A: what each pass does for the two guilds, in all five configurations ====================================
# A pass is a bridge when it concentrates lowland current and a barrier when it diverts montane current. Red and blue
# are those of the current, as in Fig. 2; amber is a pass that is a bridge for one guild and a barrier for the other.
cont = (u["rank"] == "continental").to_numpy()
br = (u.trans_robust == "above 1").to_numpy()
ba = (u.longi_robust == "below 1").to_numpy()
CAT = {"bridge for lowland": (br & ~ba, CONC), "barrier for montane": (~br & ba, DIV), "both": (br & ba, BOTH)}
fig = plt.figure(figsize=(183 / 25.4, 150 / 25.4))
ax = fig.add_axes([0.0, 0.02, 0.47, 0.94])
frame(ax)
for rk, ms in (("tributary", 0.28), ("terminal", 1.6), ("continental", 3.6)):
    rm = (u["rank"] == rk).to_numpy()
    for j, (m, col) in enumerate(CAT.values()):
        ax.plot(ux[m & rm], uy[m & rm], ls="", marker="o", ms=ms, mfc=col, mec="none", rasterized=True,
                zorder=2 + 3 * ms + j)                          # "both" on top, and the larger ranks over the smaller
n = {k: int(m.sum()) for k, (m, _) in CAT.items()}
ax.legend([Line2D([], [], ls="", marker="o", ms=3.6, mfc=c, mec="none") for _, c in CAT.values()],
          [f"{k} ({v:,})" for k, v in n.items()], loc="center right", fontsize=5.5,
          frameon=False, alignment="left", title_fontsize=5.5,
          title="What the pass does, in all five\nconfigurations (size is its rank)")
panel(ax, "a", f"Passes by what they do for each guild ({sum(n.values()):,} of {len(u):,} cells)")
a = fig.add_axes([0.58, 0.60, 0.40, 0.34])
bins = np.arange(0, 6301, 100)
a.hist(u.z_m, bins=bins, histtype="stepfilled", color="0.88", label=f"every pass ({len(u):,})")
a.hist([u.z_m[m] for m, _ in CAT.values()], bins=bins, histtype="barstacked", stacked=True,
       color=[c for _, c in CAT.values()], lw=0, label=list(CAT))
a.axvline(2000, c="0.5", lw=0.5, ls=":")
a.set_xlim(0, 6300); a.set_xlabel("Height of the pass (m)"); a.set_ylabel("Cells with passes")
a.spines[["top", "right"]].set_visible(False); a.legend(fontsize=5.5, frameon=False)
panel(a, "b", "Height")
a = fig.add_axes([0.58, 0.10, 0.40, 0.36])
d = pd.DataFrame({"band": np.floor(u.lat / 0.5) * 0.5 + 0.25, **{k: m for k, (m, _) in CAT.items()}})
sh = d.groupby("band").mean()
left = np.zeros(len(sh))
for k, (_, col) in CAT.items():
    a.barh(sh.index, sh[k], left=left, height=0.5, color=col, lw=0)
    left = left + sh[k].values
a.set_ylim(LS, LN); a.set_yticks(LAT_TICKS); a.set_yticklabels([lat_label(v) for v in LAT_TICKS])
a.set_xlim(0, 1); a.set_xlabel("Share of the cells with passes of the band")
a.spines[["top", "right"]].set_visible(False)
panel(a, "c", "Share by latitude")
save(fig, f"{OUT}/FigA_puente_y_barrera")

# ==== Fig. B: zone conductance, both readings, and its profile along the continental divide =======================
fig = plt.figure(figsize=(183 / 25.4, 125 / 25.4))
mw = 0.30
for i, (lect, title) in enumerate((("alta", "Montane climates (longitudinal)"),
                                   ("baja", "Lowland climates (transverse)"))):
    ax = fig.add_axes([i * (mw + 0.015), 0.10, mw, 0.84])
    with rasterio.open(f"{OUT}/zona_{lect}_bio01.tif") as s:
        z = s.read(1)
    frame(ax, labels=i == 0)
    im = ax.imshow(z[::2, ::2], extent=RAS, cmap=BARRIER, vmin=0, vmax=1, interpolation="nearest", rasterized=True,
                   zorder=1)
    ax.scatter(ux[cont], uy[cont], s=0.12, c="w", lw=0, rasterized=True, zorder=2)          # the continental divide
    panel(ax, "ab"[i], title)
cb = fig.colorbar(im, cax=fig.add_axes([0.13, 0.045, 0.34, 0.016]), orientation="horizontal")
cb.set_label("Mean relative conductance of the tiles that cover the cell, BIO1 (low is a barrier)")
a = fig.add_axes([0.70, 0.10, 0.27, 0.84])
c = u[cont].copy()
c["band"] = np.floor(c.lat / 0.5) * 0.5 + 0.25
for p, ls, lab in (("longi", "-", "longitudinal, montane"), ("trans", "--", "transverse, lowland")):
    g = c.groupby("band")[[f"{p}_zone", f"{p}_zone_lo", f"{p}_zone_hi"]].median().dropna()
    a.fill_betweenx(g.index, g[f"{p}_zone_lo"], g[f"{p}_zone_hi"], color="0.82", lw=0, alpha=0.6, zorder=1)
    a.plot(g[f"{p}_zone"], g.index, c="k", ls=ls, lw=1.1, label=lab, zorder=2)
a.set_xlim(0, 1); a.set_ylim(LS, LN); a.set_yticks(LAT_TICKS); a.set_yticklabels([lat_label(v) for v in LAT_TICKS])
a.set_xlabel("Zone conductance at the continental passes"); a.spines[["top", "right"]].set_visible(False)
a.legend(fontsize=5.5, frameon=False, loc="lower left", title="Median per 0.5° band; grey, range\nover the five configurations",
         title_fontsize=5, alignment="left")
panel(a, "c", "Along the continental divide")
save(fig, f"{OUT}/FigB_conductancia_zona")
