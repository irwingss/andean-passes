"""54_contraste_gremios.py - Contrast between the two guilds: where the same piece of landscape concentrates the current
of one guild and diverts that of the other (exploratory, decided after seeing the maps [post hoc]).
A warm inter-Andean valley carries the current of lowland climates across the cordillera and cuts the current of
montane climates along it; a cold puna does the opposite. Per 1 km cell and group of variables (temperature, BIO1 to
BIO11; precipitation, BIO12 to BIO19), with the group medians of normalized current of 51_sintesis_conectividad.py:
  contrast = log10(transverse) - log10(longitudinal)
Positive: the cell favours lowland climates over montane ones (corridor for the first, barrier for the second);
negative: the reverse. Zero when both guilds see the cell alike. A cell is a robust reciprocal filter when, in all
five configurations, its transverse current is above 1 and its longitudinal below 1 (lowland corridor, montane
barrier) or the reverse (montane corridor, lowland barrier), the same rule of robustness as for passes (C8).
Reads  data/interim/conectividad/<run>/grupos/<group>_<reading>_<tag>.tif and the census.
Writes data/interim/conectividad/<run>/contraste/contraste_<group>_<tag>.tif (float, the contrast)
       data/interim/conectividad/<run>/contraste/filtro_<group>.tif (int8: 1 lowland corridor and montane barrier,
                                     -1 the reverse, 0 not robust, -128 undefined in some configuration)
       outputs/supp/tables/S21_guild_contrast_by_latitude.csv, S22_guild_contrast_by_unit.csv,
       S23_guild_contrast_control_passes.csv
       outputs/supp/FigS8_guild_contrast.{png,pdf}
Statistics are over the cells of the tropical GMBA units (11°N to 23.5°S); the maps show the whole domain.
python -P scripts/54_contraste_gremios.py principal
"""
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
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from pyproj import Transformer
from rasterio import features
from rasterio.windows import from_bounds

RUN = sys.argv[1]
D = f"data/interim/conectividad/{RUN}"
OUT = f"{D}/contraste"
SUP, SUPF = "outputs/supp/tables", "outputs/supp"
os.makedirs(OUT, exist_ok=True)
LN, LS = 11.0, -23.5
CFG = [("s200_env", "main"), ("s200_env_q0.05", "q 0.05"), ("s200_env_q0.2", "q 0.2"), ("s100_env", "100 km"),
       ("s400_env", "400 km")]
GRP = {"temp": "temperature", "prec": "precipitation"}
# colours: red and blue mean concentrated and diverted current in the other figures, so the contrast between guilds
# gets its own pair, orange for the side that favours lowland climates and purple for the side that favours montane
LOW, MON = "#E66101", "#5E3C99"
CMAP = LinearSegmentedColormap.from_list("con", [MON, "#B2ABD2", "#F7F7F7", "#FDB863", LOW])
LIM = 1.0                                                 # colour scale from -1 to 1 (a factor of 10 either way)
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.5, "axes.titlesize": 7.5,
                     "axes.titleweight": "bold", "pdf.fonttype": 42, "savefig.dpi": 600})
HALO = [matplotlib.patheffects.withStroke(linewidth=1.6, foreground="w")]


def read(path):
    with rasterio.open(path) as s:
        return s.read(1).astype(np.float64), s.profile


# ---- contrast per configuration, and the robust reciprocal filter ---------------------------------------------
C, ABOVE, BELOW, OK = {}, {}, {}, {}
for g in GRP:
    up = down = None
    ok = None
    for tag, lab in CFG:
        t, prof = read(f"{D}/grupos/{g}_baja_{tag}.tif")
        l, _ = read(f"{D}/grupos/{g}_alta_{tag}.tif")
        with np.errstate(all="ignore"):
            c = np.log10(t) - np.log10(l)
        c[~(np.isfinite(t) & np.isfinite(l) & (t > 0) & (l > 0))] = np.nan
        C[(g, lab)] = c
        p = prof.copy(); p.update(dtype="float32", nodata=np.nan, compress="deflate", count=1)
        with rasterio.open(f"{OUT}/contraste_{g}_{tag}.tif", "w", **p) as o:
            o.write(c.astype(np.float32), 1)
        f = np.isfinite(c)
        a, b = f & (t > 1) & (l < 1), f & (t < 1) & (l > 1)
        up = a if up is None else up & a
        down = b if down is None else down & b
        ok = f if ok is None else ok & f
    k = np.where(up, 1, np.where(down, -1, 0)).astype(np.int8)
    k[~ok] = -128
    p = prof.copy(); p.update(dtype="int8", nodata=-128, compress="deflate", count=1)
    with rasterio.open(f"{OUT}/filtro_{g}.tif", "w", **p) as o:
        o.write(k, 1)
    ABOVE[g], BELOW[g], OK[g] = up, down, ok
PROF = prof
print("contrast written", OUT)

# ---- grid, GMBA units of the tropics, latitude of every cell ---------------------------------------------------
with rasterio.open(f"{D}/insumos/elev.tif") as s:
    CRS, TR, H, W = s.crs, s.transform, s.height, s.width
    ELEV = s.read(1, masked=True)
to_m = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
to_g = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
_raw = pyogrio.raw.read("data/interim/M1_gmba_andes.gpkg")
gmba = shapely.from_wkb(_raw[2][0])
band = shapely.intersection(gmba, shapely.box(-90, LS, -55, LN))
band_m = shapely.transform(band, lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1])))
IN = features.rasterize([(band_m, 1)], out_shape=(H, W), transform=TR, dtype=np.uint8).astype(bool)
cols, rows = np.meshgrid(np.arange(W), np.arange(H))
xx, yy = rasterio.transform.xy(TR, rows.ravel(), cols.ravel())
LON, LAT = (np.asarray(v).reshape(H, W) for v in to_g.transform(np.asarray(xx), np.asarray(yy)))
del xx, yy, cols, rows
Z = ELEV.filled(np.nan)

# ---- S21: by 1 degree of latitude, within the tropical GMBA units ----------------------------------------------
rows = []
for g, gl in GRP.items():
    c0 = C[(g, "main")]
    for b in np.arange(LN, LS, -1.0):
        m = IN & (LAT <= b) & (LAT > b - 1)
        v = c0[m]
        f = np.isfinite(v)
        if not f.any():
            continue
        r = {"group": gl, "lat_north": b, "lat_south": b - 1, "cells": int(f.sum()),
             "median": np.median(v[f]), "p10": np.quantile(v[f], 0.1), "p90": np.quantile(v[f], 0.9)}
        n = OK[g][m].sum()
        r["cells_all_configs"] = int(n)
        r["share_lowland_corridor_montane_barrier"] = ABOVE[g][m].sum() / n if n else np.nan
        r["share_montane_corridor_lowland_barrier"] = BELOW[g][m].sum() / n if n else np.nan
        rows.append(r)
s21 = pd.DataFrame(rows)
s21.round(4).to_csv(f"{SUP}/S21_guild_contrast_by_latitude.csv", index=False)

# ---- S22: by GMBA unit (the unit of each cell from the census passes is not a polygon, so the cells take the unit of
# the nearest tropical pass in the same 1 km cell or, where there is none, are left out) ---------------------------
passes = pd.read_csv("outputs/database/APD_passes.csv.gz", usecols=["apd_id", "lat", "lon", "z_m", "depth_m", "gmba",
                                                                    "tropics", "cell_1km"])
pt = passes[passes.tropics == 1].copy()
r_, c_ = rasterio.transform.rowcol(TR, *to_m.transform(pt.lon.values, pt.lat.values))
pt["row"], pt["col"] = np.asarray(r_), np.asarray(c_)
for g, gl in GRP.items():
    pt[f"{g}_contrast"] = C[(g, "main")][pt.row, pt.col]
    pt[f"{g}_contrast_min"] = np.nanmin([C[(g, l)][pt.row, pt.col] for _, l in CFG], axis=0)
    pt[f"{g}_contrast_max"] = np.nanmax([C[(g, l)][pt.row, pt.col] for _, l in CFG], axis=0)
    k = np.where(ABOVE[g], 1, np.where(BELOW[g], -1, 0))[pt.row, pt.col]
    k = np.where(OK[g][pt.row, pt.col], k, -9)
    pt[f"{g}_filter"] = pd.Series(k).map({1: "lowland corridor, montane barrier", -1: "montane corridor, lowland "
                                          "barrier", 0: "not robust", -9: ""}).values
one = pt.sort_values("depth_m", ascending=False).drop_duplicates("cell_1km")     # one pass per cell, the deepest
rows = []
for (u, g), x in ((ug, one) for ug in [(u, g) for u in sorted(one.gmba.dropna().unique()) for g in GRP]):
    x = x[x.gmba == u]
    v = x[f"{g}_contrast"].dropna()
    if len(v) == 0:
        continue
    f = x[f"{g}_filter"]
    rows.append({"gmba": u, "group": GRP[g], "cells_with_passes": len(v), "median_contrast": v.median(),
                 "p10": v.quantile(0.1), "p90": v.quantile(0.9),
                 "share_lowland_corridor_montane_barrier": (f == "lowland corridor, montane barrier").mean(),
                 "share_montane_corridor_lowland_barrier": (f == "montane corridor, lowland barrier").mean()})
pd.DataFrame(rows).round(4).to_csv(f"{SUP}/S22_guild_contrast_by_unit.csv", index=False)

# ---- S23: the control passes -------------------------------------------------------------------------------------
ctl = pd.read_csv("data/refs/pasos_control.csv")
ctl = ctl[(ctl.in_area == 1) & (ctl.lat <= LN) & (ctl.lat >= LS)].reset_index(drop=True)
cr, cc = rasterio.transform.rowcol(TR, *to_m.transform(ctl.lon.values, ctl.lat.values))
cr, cc = np.asarray(cr), np.asarray(cc)
s23 = ctl[["name", "system", "lat", "lon"]].copy()
s23.insert(0, "no", np.arange(1, len(ctl) + 1))
for g, gl in GRP.items():
    vals = np.array([C[(g, l)][cr, cc] for _, l in CFG])
    s23[f"{g}_contrast"] = vals[0]
    s23[f"{g}_min"], s23[f"{g}_max"] = np.nanmin(vals, 0), np.nanmax(vals, 0)
    k = np.where(ABOVE[g][cr, cc], "lowland corridor, montane barrier",
                 np.where(BELOW[g][cr, cc], "montane corridor, lowland barrier", "not robust"))
    s23[f"{g}_filter"] = np.where(OK[g][cr, cc], k, "")
s23.round(3).to_csv(f"{SUP}/S23_guild_contrast_control_passes.csv", index=False)

# ---- numbers for the text --------------------------------------------------------------------------------------
for g, gl in GRP.items():
    m = IN & OK[g]
    za, zb, zr = Z[IN & ABOVE[g]], Z[IN & BELOW[g]], Z[m]
    rho = pd.Series(C[(g, "main")][IN]).corr(pd.Series(Z[IN]), method="spearman")
    print(f"{gl}: {m.sum():,} GMBA cells defined in all five; lowland corridor and montane barrier "
          f"{ABOVE[g][m].mean():.3f} (median height {np.nanmedian(za):.0f} m), montane corridor and lowland barrier "
          f"{BELOW[g][m].mean():.3f} (median height {np.nanmedian(zb):.0f} m), all cells {np.nanmedian(zr):.0f} m; "
          f"Spearman contrast vs height {rho:.3f}")

# ==== Fig. S8 =====================================================================================================
EXT0 = band_m.bounds
EXT = (EXT0[0] - 120e3, EXT0[2] + 120e3, EXT0[1] - 60e3, EXT0[3] + 60e3)
win = from_bounds(EXT[0], EXT[2], EXT[1], EXT[3], TR).round_offsets().round_lengths()
e = rasterio.windows.bounds(win, TR)
EXT = (e[0], e[2], e[1], e[3])


def crop(a):
    return a[win.row_off:win.row_off + win.height, win.col_off:win.col_off + win.width]


def hillshade(e, az=315, alt=40, zf=8):
    gy, gx = np.gradient(np.asarray(e, float) * zf / 1000.0)
    sl, asp = np.arctan(np.hypot(gx, gy)), np.arctan2(-gx, gy)
    a, z = np.radians(360 - az + 90), np.radians(90 - alt)
    return np.cos(z) * np.cos(sl) + np.sin(z) * np.sin(sl) * np.cos(a - asp)


HS = np.ma.masked_invalid(hillshade(crop(Z)))
PAIS = [shapely.from_wkb(w) for w in pyogrio.raw.read("data/interim/M1_paises.gpkg")[2]]
PAIS_M = [shapely.transform(p, lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1]))) for p in PAIS]
cx, cy = to_m.transform(ctl.lon.values, ctl.lat.values)
ZOOM = (("Colombia and Táchira", (-77.8, -70.5, 0.6, 9.2)), ("Huancabamba Depression", (-80.6, -77.4, -7.5, -4.0)),
        ("Southern Peru", (-75.0, -66.5, -16.5, -11.3)))


def relief(ax):
    ax.imshow(HS, extent=EXT, cmap="Greys_r", vmin=-0.6, vmax=1.6, interpolation="bilinear", rasterized=True,
              zorder=0)


def frame(ax, labels=False, lw=0.25):
    for p in PAIS_M:
        for gg in getattr(p, "geoms", [p]):
            ax.plot(*gg.exterior.xy, lw=lw, c="0.45", alpha=0.85, zorder=4.5)
    for gg in getattr(band_m, "geoms", [band_m]):
        ax.plot(*gg.exterior.xy, lw=0.35, c="k", alpha=0.6, zorder=4.6)
    for la in range(-20, 11, 5):
        lo = np.linspace(-85, -55, 200)
        x, y = to_m.transform(lo, np.full_like(lo, la))
        ax.plot(x, y, lw=0.25, c="0.6", ls=":", zorder=0.5)
        if labels:
            ax.text(EXT[0] + 15e3, np.interp(EXT[0] + 15e3, x, y),
                    f"{abs(la)}°{'S' if la < 0 else 'N' if la else ''}", fontsize=5.5, color="0.35", va="bottom",
                    path_effects=HALO)
    ax.scatter(cx, cy, s=7, marker="^", fc="w", ec="k", lw=0.45, zorder=5)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")


def box_m(b):
    w, e_, s, n = b
    return to_m.transform([w, e_, e_, w, w], [s, s, n, n, s])


def panel(ax, letter, text=""):
    ax.annotate(letter, (0, 1), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points", fontsize=8,
                fontweight="bold", va="bottom", ha="left")
    if text:
        ax.annotate(text, (0, 1), xycoords="axes fraction", xytext=(9, 3), textcoords="offset points", fontsize=6.5,
                    va="bottom", ha="left")


CM = {g: np.ma.masked_invalid(crop(C[(g, "main")])) for g in GRP}
FK = {}
for g in GRP:
    k = np.where(ABOVE[g], 1.0, np.where(BELOW[g], -1.0, 0.0))
    k[~OK[g]] = np.nan
    FK[g] = np.ma.masked_invalid(crop(k))
KMAP = ListedColormap([MON, "#E8E8E8", LOW])
KNORM = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], 3)

Wd, gap, gt, bot = 183.0, 2.5, 6.0, 16.0
mw = (Wd - 3 * gap) / 4
mh = mw * (EXT[3] - EXT[2]) / (EXT[1] - EXT[0])
ZB = [shapely.Polygon(np.column_stack(box_m(b))).bounds for _, b in ZOOM]
zw = (Wd - 3 * gap) * 0.75 / 3
zh = max(zw * (y1 - y0) / (x1 - x0) for x0, y0, x1, y1 in ZB)
MID = 26.0                                                # colour bar and key between the rows
Ht = bot + zh + MID + mh + gt
fig = plt.figure(figsize=(Wd / 25.4, Ht / 25.4))
F = lambda x, y, w, h: fig.add_axes([x / Wd, y / Ht, w / Wd, h / Ht])
ytop = Ht - gt - mh
lab = iter("abcdefgh")
for i, (kind, g) in enumerate([("c", "temp"), ("c", "prec"), ("k", "temp"), ("k", "prec")]):
    a = F(i * (mw + gap), ytop, mw, mh)
    relief(a)
    if kind == "c":
        im = a.imshow(CM[g], extent=EXT, cmap=CMAP, vmin=-LIM, vmax=LIM, interpolation="nearest", rasterized=True,
                      alpha=0.92)
        txt = f"Contrast, {GRP[g]}"
    else:
        ik = a.imshow(FK[g], extent=EXT, cmap=KMAP, norm=KNORM, interpolation="nearest", rasterized=True, alpha=0.9)
        txt = f"Robust filters, {GRP[g]}"
    frame(a, labels=i == 0)
    if i == 0:
        for _, b in ZOOM:
            a.plot(*box_m(b), c="k", lw=0.6, zorder=6)
    a.set_xlim(EXT[0], EXT[1]); a.set_ylim(EXT[2], EXT[3])
    panel(a, next(lab), txt)
cax = F(0.04 * Wd, ytop - 12, 0.42 * Wd, 2.4)
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_ticks([-1, -0.5, 0, 0.5, 1]); cb.set_ticklabels(["0.1", "0.3", "1", "3", "10"])
cb.set_label("Transverse over longitudinal normalized current (median of the group, log scale)", fontsize=6)
cax.annotate("favours montane climates", (0, 0.5), xycoords="axes fraction", xytext=(0, 5), textcoords="offset points",
             ha="left", va="bottom", fontsize=5.5, color=MON)
cax.annotate("favours lowland climates", (1, 0.5), xycoords="axes fraction", xytext=(0, 5), textcoords="offset points",
             ha="right", va="bottom", fontsize=5.5, color=LOW)
kx = 0.54 * Wd
for j, (col, t) in enumerate(((LOW, "lowland corridor and montane barrier in all five configurations"),
                              (MON, "montane corridor and lowland barrier in all five"),
                              ("#E8E8E8", "not robust"))):
    a = F(kx, ytop - 9 - j * 3.4, 3, 2.4)
    a.add_patch(plt.Rectangle((0, 0), 1, 1, fc=col, ec="0.5", lw=0.3)); a.axis("off")
    fig.text((kx + 4) / Wd, (ytop - 7.8 - j * 3.4) / Ht, t, fontsize=5.5, va="center")
# zooms of the temperature contrast and the latitudinal profile
for j, ((name, _), (x0, y0, x1, y1)) in enumerate(zip(ZOOM, ZB)):
    h = zw * (y1 - y0) / (x1 - x0)
    a = F(j * (zw + gap), bot + zh - h, zw, h)
    relief(a)
    a.imshow(CM["temp"], extent=EXT, cmap=CMAP, vmin=-LIM, vmax=LIM, interpolation="nearest", rasterized=True,
             alpha=0.92)
    frame(a, lw=0.3)
    for n, (x, y) in enumerate(zip(cx, cy), 1):
        if x0 < x < x1 and y0 < y < y1:
            a.annotate(str(n), (x, y), xytext=(3, 1), textcoords="offset points", fontsize=4.5, zorder=6,
                       path_effects=HALO)
    a.set_xlim(x0, x1); a.set_ylim(y0, y1)
    panel(a, next(lab), f"{name}, temperature")
a = F(3 * (zw + gap) + 6, bot, Wd - 3 * (zw + gap) - 6, zh)
for g, ls in (("temp", "-"), ("prec", "--")):
    x = s21[s21.group == GRP[g]]
    yc = (x.lat_north + x.lat_south) / 2
    a.plot(x.share_lowland_corridor_montane_barrier, yc, c=LOW, ls=ls, lw=1)
    a.plot(x.share_montane_corridor_lowland_barrier, yc, c=MON, ls=ls, lw=1)
a.plot([], [], c="0.3", ls="-", lw=1, label="temperature"); a.plot([], [], c="0.3", ls="--", lw=1,
                                                                    label="precipitation")
a.legend(frameon=False, fontsize=5.5, loc="upper right", handlelength=2.2)
a.set_ylim(LS, LN); a.set_xlim(0, None)
a.set_yticks(range(-20, 11, 5)); a.set_yticklabels([f"{abs(v)}°{'S' if v < 0 else 'N' if v else ''}"
                                                    for v in range(-20, 11, 5)])
a.set_xlabel("Share of GMBA cells, per 1° band", fontsize=6)
a.spines[["top", "right"]].set_visible(False)
panel(a, next(lab), "Robust filters")
fig.savefig(f"{SUPF}/FigS8_guild_contrast.png", bbox_inches="tight", pad_inches=0.03)
fig.savefig(f"{SUPF}/FigS8_guild_contrast.pdf", bbox_inches="tight", pad_inches=0.03)
print("->", f"{SUPF}/FigS8_guild_contrast")
