"""52_figuras.py - Figures of the manuscript and the supplement for the census and the present connectivity.
Reads the census (data/interim/censo), the connectivity run (data/interim/conectividad/<run>) and the tables of
51_sintesis_conectividad.py. Writes PNG (600 dpi) and PDF to outputs/figures (Fig1 to Fig3) and outputs/supp (FigS1 to
FigS6).
python -P scripts/52_figuras.py principal
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
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from pyproj import Transformer
from rasterio.windows import from_bounds

RUN = sys.argv[1]
D = f"data/interim/conectividad/{RUN}"
TS = f"{D}/teselas"
FIG, SUPF = "outputs/figures", "outputs/supp"
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                     "ytick.major.width": 0.5, "axes.titlesize": 7.5, "axes.titleweight": "bold",
                     "pdf.fonttype": 42, "savefig.dpi": 600})
LN, LS = 11.0, -23.5
VARS = [f"bio{v:02d}" for v in range(1, 20)]
VLAB = {"bio01": "BIO1 annual mean temperature", "bio02": "BIO2 mean diurnal range", "bio03": "BIO3 isothermality",
        "bio04": "BIO4 temperature seasonality", "bio05": "BIO5 max. temperature of warmest month",
        "bio06": "BIO6 min. temperature of coldest month", "bio07": "BIO7 temperature annual range",
        "bio08": "BIO8 mean temperature of wettest quarter", "bio09": "BIO9 mean temperature of driest quarter",
        "bio10": "BIO10 mean temperature of warmest quarter", "bio11": "BIO11 mean temperature of coldest quarter",
        "bio12": "BIO12 annual precipitation", "bio13": "BIO13 precipitation of wettest month",
        "bio14": "BIO14 precipitation of driest month", "bio15": "BIO15 precipitation seasonality",
        "bio16": "BIO16 precipitation of wettest quarter", "bio17": "BIO17 precipitation of driest quarter",
        "bio18": "BIO18 precipitation of warmest quarter", "bio19": "BIO19 precipitation of coldest quarter"}
RANK_COL = {"continental": "#D55E00", "terminal": "#0072B2", "tributary": "#8C8C8C"}      # Okabe-Ito, plus size
CUR = LinearSegmentedColormap.from_list("cur", ["#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"])
NORM = LogNorm(0.2, 5)


def save(fig, path):
    fig.savefig(path + ".png", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(path + ".pdf", bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("->", path)


# ---- shared map frame: the 1 km equal-area grid, cropped to the tropics ------------------------------------------
with rasterio.open(f"{D}/insumos/elev.tif") as s:
    CRS, TR, H, W = s.crs, s.transform, s.height, s.width
    ELEV = s.read(1, masked=True)
to_m = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
gmba = shapely.from_wkb(pyogrio.raw.read("data/interim/M1_gmba_andes.gpkg")[2][0])
band = shapely.intersection(gmba, shapely.box(-90, LS, -55, LN))
band_m = shapely.transform(band, lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1])))
bx0, by0, bx1, by1 = band_m.bounds
EXT = (bx0 - 120e3, bx1 + 120e3, by0 - 60e3, by1 + 60e3)
win = from_bounds(EXT[0], EXT[2], EXT[1], EXT[3], TR).round_offsets().round_lengths()
EXT = rasterio.windows.bounds(win, TR)
EXT = (EXT[0], EXT[2], EXT[1], EXT[3])
ctl = pd.read_csv("data/refs/pasos_control.csv")
ctl = ctl[(ctl.in_area == 1) & (ctl.lat <= LN) & (ctl.lat >= LS)]
ctl_x, ctl_y = to_m.transform(ctl.lon.values, ctl.lat.values)


def crop(a):
    return a[win.row_off:win.row_off + win.height, win.col_off:win.col_off + win.width]


def frame(ax, labels=True, names=False):
    """GMBA outline of the tropics, graticule every 5 degrees and the control passes."""
    for g in getattr(band_m, "geoms", [band_m]):
        ax.plot(*g.exterior.xy, lw=0.35, c="k", alpha=0.6)
    for la in range(-20, 11, 5):
        lo = np.linspace(-85, -55, 200)
        x, y = to_m.transform(lo, np.full_like(lo, la))
        ax.plot(x, y, lw=0.25, c="0.6", ls=":", zorder=0.5)
        if labels:
            ax.text(EXT[0] + 15e3, np.interp(EXT[0] + 15e3, x, y), f"{abs(la)}°{'S' if la < 0 else 'N' if la else ''}",
                    fontsize=5.5, color="0.35", va="bottom")
    for lo in range(-80, -59, 5):
        la = np.linspace(-30, 15, 200)
        x, y = to_m.transform(np.full_like(la, lo), la)
        ax.plot(x, y, lw=0.25, c="0.6", ls=":", zorder=0.5)
        if labels and EXT[0] < np.interp(EXT[2] + 15e3, y, x) < EXT[1] - 60e3:
            ax.text(np.interp(EXT[2] + 15e3, y, x), EXT[2] + 15e3, f"{abs(lo)}°W", fontsize=5.5, color="0.35",
                    ha="center")
    ax.scatter(ctl_x, ctl_y, s=9, marker="^", fc="w", ec="k", lw=0.5, zorder=5)
    if names:
        for x, y, n in zip(ctl_x, ctl_y, range(1, len(ctl) + 1)):     # numbers of Table 2
            ax.annotate(str(n), (x, y), fontsize=4.5, xytext=(3, 1), textcoords="offset points", zorder=6)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")


def hillshade(e, az=315, alt=40, zf=8):
    """Shaded relief of the 1 km FABDEM grid (display only); nan stays nan."""
    gy, gx = np.gradient(np.asarray(e.filled(np.nan), float) * zf / 1000.0)
    sl, asp = np.arctan(np.hypot(gx, gy)), np.arctan2(-gx, gy)
    a, z = np.radians(360 - az + 90), np.radians(90 - alt)
    return np.cos(z) * np.cos(sl) + np.sin(z) * np.sin(sl) * np.cos(a - asp)


HS = hillshade(crop(ELEV))


def relief(ax, light=False):
    """Shaded relief under the data: cells without a value show relief, never a colour of the scale."""
    ax.imshow(np.ma.masked_invalid(HS), extent=EXT, cmap="Greys_r", vmin=-0.6 if light else -0.2,
              vmax=1.6 if light else 1.15, interpolation="bilinear", rasterized=True, zorder=0)


def current(f, step=1):
    """Normalized current, undefined where the current under uniform resistance is below 1e-6 A (as in 51)."""
    with rasterio.open(f) as s:
        a = s.read(1, masked=True, window=win)
    with rasterio.open(f.replace("_normalizada", "_corriente")) as s:
        c = s.read(1, masked=True, window=win)
    a = np.ma.masked_where(np.ma.getdata(c) / np.where(np.ma.getdata(a) > 0, np.ma.getdata(a), np.inf) < 1e-6, a)
    return a[::step, ::step]


# ---- data ----------------------------------------------------------------------------------------------------------
db = pd.read_csv("outputs/database/APD_cols.csv.gz", dtype={"trans_robust": str, "longi_robust": str})
tr = db[db.tropics == 1]
tr_x, tr_y = to_m.transform(tr.lon.values, tr.lat.values)

# ==== Fig. 1: the census ============================================================================================
fig = plt.figure(figsize=(183 / 25.4, 168 / 25.4))
ax = fig.add_axes([0.0, 0.02, 0.46, 0.96])
relief(ax, light=True)
d = tr.depth_m >= 100
for rk, sz in (("tributary", 0.25), ("terminal", 0.8), ("continental", 1.6)):
    m = (d & (tr["rank"] == rk)).to_numpy()
    ax.scatter(tr_x[m], tr_y[m], s=sz, c=RANK_COL[rk], lw=0, rasterized=True, zorder=2,
               label=f"{rk} ({m.sum():,})")
frame(ax, names=True)
ax.legend(loc="lower left", fontsize=5.5, markerscale=4, frameon=False, title="Cols with depth of 100 m or more",
          title_fontsize=5.5)
ax.set_title("a", loc="left")
z0, z1, w0, w1 = -7.5, -4.0, -80.6, -77.4
bx = to_m.transform([w0, w1, w1, w0, w0], [z0, z0, z1, z1, z0])
ax.plot(*bx, c="k", lw=0.7, zorder=6)
# b: divide network and every col around the Huancabamba Depression, from the 1 arc-second census
axb = fig.add_axes([0.50, 0.30, 0.49, 0.68])
NETC = LinearSegmentedColormap.from_list("n", ["#CFCFCF", "#CFCFCF"])
for blk, (top, bot) in (("S04p00", (-4.0, -5.0)), ("S05p00", (-5.0, -6.0)), ("S06p00", (-6.0, -7.0)),
                        ("S07p00", (-7.0, -7.5))):
    with rasterio.open(f"data/interim/censo/divisorias_{blk}.tif") as s:
        wn = from_bounds(w0, bot, w1, top, s.transform).intersection(
            rasterio.windows.Window(0, 0, s.width, s.height)).round_offsets().round_lengths()
        a = s.read(1, window=wn) == 1
        l, b, r, t = rasterio.windows.bounds(wn, s.transform)
    k = 4                                                       # 4 by 4 cells, about 120 m, for display only
    a = a[:a.shape[0] // k * k, :a.shape[1] // k * k]
    a = a.reshape(a.shape[0] // k, k, a.shape[1] // k, k).any(axis=(1, 3))
    axb.imshow(np.ma.masked_where(~a, a), extent=(l, l + a.shape[1] * k / 3600, t - a.shape[0] * k / 3600, t),
               cmap=NETC, interpolation="nearest", rasterized=True, aspect="auto")
z = tr[(tr.lat >= z0) & (tr.lat <= z1) & (tr.lon >= w0) & (tr.lon <= w1)]
for rk in ("tributary", "terminal", "continental"):
    m = z[(z["rank"] == rk) & ((rk != "tributary") | (z.depth_m >= 100))]
    axb.scatter(m.lon, m.lat, s=np.clip(m.depth_m / 60, 0.5, 14), c=RANK_COL[rk], lw=0, alpha=0.85,
                rasterized=True, label=f"{rk} ({len(m):,})")
k = ctl[(ctl.lat >= z0) & (ctl.lat <= z1) & (ctl.lon >= w0) & (ctl.lon <= w1)]
axb.scatter(k.lon, k.lat, s=22, marker="^", fc="w", ec="k", lw=0.7, zorder=5)
for _, r in k.iterrows():
    axb.annotate(r["name"], (r.lon, r.lat), fontsize=6, xytext=(4, 2), textcoords="offset points")
axb.set_xlim(w0, w1); axb.set_ylim(z0, z1)
axb.set_aspect(1 / np.cos(np.radians(5.75)))
axb.set_xlabel("Longitude (°)"); axb.set_ylabel("Latitude (°)")
axb.legend(loc="lower left", fontsize=5.5, markerscale=1.5, frameon=True, framealpha=0.8, edgecolor="none",
           title="Cols sized by depth", title_fontsize=5.5)
axb.set_title("b", loc="left")
# c: depth against height of the tropical cols with depth of 100 m or more
axc = fig.add_axes([0.58, 0.05, 0.38, 0.19])
for rk in ("tributary", "terminal", "continental"):
    m = tr[d & (tr["rank"] == rk)]
    axc.scatter(m.z_m, m.depth_m, s=0.4 if rk == "tributary" else 1.5, c=RANK_COL[rk], lw=0, rasterized=True)
axc.set_yscale("log"); axc.set_xlabel("Col height (m)"); axc.set_ylabel("Depth (m)")
axc.set_title("c", loc="left")
save(fig, f"{FIG}/Fig1_census")

# ==== Fig. 2: normalized current with BIO1, both readings, and the Huancabamba Depression ========================
fig = plt.figure(figsize=(183 / 25.4, 160 / 25.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.92], height_ratios=[1, 1], wspace=0.04, hspace=0.06)
zx0, zx1 = to_m.transform([-81.3, -76.4], [-3.0, -3.0])[0]
zy0, zy1 = to_m.transform([-78.8, -78.8], [-8.6, -2.6])[1]
for i, (l, t) in enumerate((("baja", "Transverse reading\nlowland guild, west to east"),
                            ("alta", "Longitudinal reading\nmontane guild, north to south"))):
    a = fig.add_subplot(gs[:, i])
    z = current(f"{TS}/bio01_{l}_s200_env_normalizada.tif")
    relief(a, light=True)
    im = a.imshow(z, extent=EXT, cmap=CUR, norm=NORM, interpolation="nearest", rasterized=True, alpha=0.92)
    frame(a, labels=i == 0, names=True)
    a.plot([zx0, zx1, zx1, zx0, zx0], [zy0, zy0, zy1, zy1, zy0], c="k", lw=0.6, zorder=6)
    a.set_title(f"{'ab'[i]}", loc="left"); a.set_title(t, fontsize=6.5, fontweight="normal")
    b = fig.add_subplot(gs[i, 2])
    relief(b, light=True)
    b.imshow(z, extent=EXT, cmap=CUR, norm=NORM, interpolation="nearest", rasterized=True, alpha=0.92)
    m = (tr.depth_m >= 300).to_numpy()
    b.scatter(tr_x[m], tr_y[m], s=0.7, c="k", lw=0, alpha=0.7, zorder=4, rasterized=True)
    b.scatter(ctl_x, ctl_y, s=16, marker="^", fc="w", ec="k", lw=0.6, zorder=5)
    for x, y, n in zip(ctl_x, ctl_y, ctl.name):
        if zx0 < x < zx1 and zy0 < y < zy1:
            b.annotate(n, (x, y), fontsize=5.5, xytext=(4, 2), textcoords="offset points", zorder=6,
                       path_effects=[matplotlib.patheffects.withStroke(linewidth=1.5, foreground="w")])
    b.set_xlim(zx0, zx1); b.set_ylim(zy0, zy1); b.set_xticks([]); b.set_yticks([]); b.set_aspect("equal")
    b.set_title(f"{'cd'[i]}", loc="left")
    b.set_title(f"Huancabamba, {'transverse' if i == 0 else 'longitudinal'}", fontsize=6.5, fontweight="normal",
                loc="right")
cb = fig.colorbar(im, ax=fig.axes, orientation="horizontal", shrink=0.42, pad=0.015, aspect=38)
cb.set_label("Normalized current with BIO1 (log scale; above 1 the climate concentrates the current, below 1 it diverts it)")
cb.minorticks_off(); cb.set_ticks([0.2, 0.5, 1, 2, 5]); cb.set_ticklabels(["0.2", "0.5", "1", "2", "5"])
save(fig, f"{FIG}/Fig2_current_bio1")

# ==== Fig. 3: latitudinal profiles and agreement between variables ==================================================
pr = pd.read_csv("outputs/supp/tables/S8_tile_conductance_profiles.csv")
sp = pd.read_csv("outputs/supp/tables/S10_spearman_variables.csv")
fig = plt.figure(figsize=(183 / 25.4, 132 / 25.4))
for i, (rd, t) in enumerate((("transverse", "a  Transverse"), ("longitudinal", "b  Longitudinal"))):
    a = fig.add_axes([0.07 + i * 0.21, 0.08, 0.18, 0.86])
    p = pr[(pr.cfg == "main") & (pr.reading == rd)]
    w = p.pivot(index="band", columns="var", values="median")
    a.fill_betweenx(w.index, w.min(1), w.max(1), color="0.85", lw=0, label="range of the 19 variables")
    a.plot(w.median(1), w.index, c="0.3", lw=0.8, label="median of the 19 variables")
    a.plot(w["bio01"], w.index, c="#B2182B", lw=0.9, label="BIO1")
    for lo, hi in ((-7.5, -4.0), (7.5, 8.0)):
        a.axhspan(lo, hi, color="#2166AC", alpha=0.12, lw=0)
    a.set_ylim(LS, LN); a.set_xlim(0, 1.05)
    a.set_xlabel("Relative conductance of the tiles")
    if i == 0:
        a.set_ylabel("Latitude (°)")
        a.legend(loc="lower left", fontsize=5.5, frameon=False)
    else:
        a.set_yticklabels([])
        a.text(0.03, -5.75, "Huancabamba", fontsize=5.5, va="center", color="#2166AC")
        a.text(0.03, 7.75, "Táchira", fontsize=5.5, va="center", color="#2166AC")
    a.set_title(t, loc="left")
for i, (rd, t) in enumerate((("transverse", "c  Transverse, cols"), ("longitudinal", "d  Longitudinal, cols"))):
    a = fig.add_axes([0.52, 0.55 - i * 0.47, 0.40, 0.40])
    m = sp[sp.reading == rd].pivot(index="var_a", columns="var_b", values="rho_cols").loc[VARS, VARS]
    im = a.imshow(m.values, cmap="RdBu_r", vmin=-1, vmax=1)
    a.set_xticks(range(19)); a.set_yticks(range(19))
    lab = [v.replace("bio0", "").replace("bio", "") for v in VARS]
    a.set_xticklabels(lab, fontsize=5); a.set_yticklabels(lab, fontsize=5)
    a.set_title(t, loc="left")
cb = fig.colorbar(im, cax=fig.add_axes([0.89, 0.25, 0.012, 0.5]))
cb.set_label("Spearman correlation between variables")
save(fig, f"{FIG}/Fig3_profiles_variables")

# ==== Fig. S1: validation of the census =============================================================================
V = "data/interim/censo/validacion"
osm = pd.read_csv(f"{V}/O5_1_nodos_osm.csv")
idn = pd.read_csv(f"{V}/O5_2_identidad.csv").set_index("group").loc["all"]
ref = pd.read_csv(f"{V}/O5_3_alturas_referencia.csv")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw={"wspace": 0.3})
x = np.sort(osm.col_km.values)
ax[0].plot(x, np.arange(1, len(x) + 1) / len(x), c="k", lw=0.9, label=f"OSM pass nodes ({len(x):,})")
ax[0].scatter([0.5, 1.0], [idn.null_500m_mean, idn.null_1km_mean], c="#2166AC", s=12, zorder=3,
              label="random divide cells, same core and unit (mean of 1000)")
ax[0].set_xlim(0, 3); ax[0].set_ylim(0, 1)
ax[0].set_xlabel("Distance to the nearest col (km)"); ax[0].set_ylabel("Cumulative share of nodes")
ax[0].legend(fontsize=5.5, frameon=False, loc="lower right"); ax[0].set_title("a", loc="left")
for zn, c in (("tropical", "#B2182B"), ("temperate", "#2166AC")):
    r = ref[ref.zone == zn]
    ax[1].scatter(r.elev_ref_m, r.error_m, s=10, c=c, label=f"{zn} ({len(r)})")
ax[1].axhline(0, c="k", lw=0.5); ax[1].axhspan(-50, 50, color="0.9", lw=0, zorder=0)
ax[1].set_xlabel("Official height of the pass (m)"); ax[1].set_ylabel("Nearest col minus official height (m)")
ax[1].legend(fontsize=5.5, frameon=False); ax[1].set_title("b", loc="left")
save(fig, f"{SUPF}/FigS1_validation")

# ==== Fig. S2 and S3: the 19 variables, each reading ================================================================
for k, (l, t) in enumerate((("baja", "transverse"), ("alta", "longitudinal"))):
    fig, ax = plt.subplots(4, 5, figsize=(7.2, 9.4), gridspec_kw={"wspace": 0.02, "hspace": 0.12})
    for a, v in zip(ax.flat, VARS):
        im = a.imshow(current(f"{TS}/{v}_{l}_s200_env_normalizada.tif", 2), extent=EXT, cmap=CUR, norm=NORM,
                      interpolation="nearest", rasterized=True)
        for g in getattr(band_m, "geoms", [band_m]):
            a.plot(*g.exterior.xy, lw=0.25, c="k", alpha=0.6)
        a.set_xlim(EXT[0], EXT[1]); a.set_ylim(EXT[2], EXT[3]); a.set_xticks([]); a.set_yticks([])
        a.set_title(v.upper().replace("BIO0", "BIO"), fontsize=6)
    ax.flat[-1].axis("off")
    cb = fig.colorbar(im, ax=ax.flat[-1], orientation="horizontal", fraction=0.5, aspect=12)
    cb.set_label(f"Normalized current,\n{t} reading", fontsize=6)
    cb.minorticks_off(); cb.set_ticks([0.2, 0.5, 1, 2, 5]); cb.set_ticklabels(["0.2", "0.5", "1", "2", "5"])
    save(fig, f"{SUPF}/FigS{2 + k}_variables_{t}")

# ==== Fig. S4: BIO1 under the four sensitivity configurations =======================================================
conf = (("s200_env_q0.05", "q 0.05, 200 km"), ("s200_env_q0.2", "q 0.2, 200 km"), ("s100_env", "q 0.1, 100 km"),
        ("s400_env", "q 0.1, 400 km"))
fig, ax = plt.subplots(2, 4, figsize=(7.2, 6.4), gridspec_kw={"wspace": 0.02, "hspace": 0.08})
for i, (l, t) in enumerate((("baja", "transverse"), ("alta", "longitudinal"))):
    for j, (f, lab) in enumerate(conf):
        a = ax[i, j]
        im = a.imshow(current(f"{TS}/bio01_{l}_{f}_normalizada.tif", 2), extent=EXT, cmap=CUR, norm=NORM,
                      interpolation="nearest", rasterized=True)
        frame(a, labels=False)
        a.set_title(f"{'abcdefgh'[i * 4 + j]}  {t}, {lab}", loc="left", fontsize=5.5)
cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.4, pad=0.02, aspect=30)
cb.set_label("Normalized current with BIO1"); cb.minorticks_off(); cb.set_ticks([0.2, 0.5, 1, 2, 5])
cb.set_ticklabels(["0.2", "0.5", "1", "2", "5"])
save(fig, f"{SUPF}/FigS4_sensitivity_maps")

# ==== Fig. S5: profiles with BIO1 under the five configurations ======================================================
fig, ax = plt.subplots(1, 2, figsize=(5.0, 5.2), sharey=True, gridspec_kw={"wspace": 0.08})
cc = {"main": "k", "q 0.05": "#2166AC", "q 0.2": "#92C5DE", "100 km": "#EF8A3C", "400 km": "#B2182B"}
for a, (rd, t) in zip(ax, (("transverse", "a  Transverse"), ("longitudinal", "b  Longitudinal"))):
    for c, col in cc.items():
        p = pr[(pr.cfg == c) & (pr.reading == rd) & (pr["var"] == "bio01")]
        a.plot(p["median"], p.band, c=col, lw=1.0 if c == "main" else 0.7, label=c)
    for lo, hi in ((-7.5, -4.0), (7.5, 8.0)):
        a.axhspan(lo, hi, color="#2166AC", alpha=0.12, lw=0)
    a.set_xlim(0, 1.05); a.set_ylim(LS, LN); a.set_xlabel("Relative conductance of the tiles, BIO1")
    a.set_title(t, loc="left")
ax[0].set_ylabel("Latitude (°)"); ax[0].legend(fontsize=5.5, frameon=False, loc="lower left")
save(fig, f"{SUPF}/FigS5_sensitivity_profiles")

# ==== Fig. S6: cols by the median over the 19 variables and its robustness ==========================================
fig, ax = plt.subplots(1, 2, figsize=(7.2, 6.4), gridspec_kw={"wspace": 0.03})
u = tr.sort_values("depth_m", ascending=False).drop_duplicates("cell_1km")      # deepest col per cell, as in 51
ux, uy = to_m.transform(u.lon.values, u.lat.values)
RC = {"depends": "#D9D9D9", "below 1": "#2166AC", "above 1": "#B2182B"}
for a, p, t in zip(ax, ("trans", "longi"), ("a  Transverse reading", "b  Longitudinal reading")):
    relief(a, light=True)
    for k in ("depends", "below 1", "above 1"):
        m = (u[f"{p}_robust"] == k).to_numpy()
        a.scatter(ux[m], uy[m], s=0.15, c=RC[k], lw=0, rasterized=True, label=f"{k} ({m.sum():,} cells)")
    frame(a, labels=a is ax[0])
    a.legend(loc="lower left", fontsize=5.5, markerscale=8, frameon=False,
             title="Median over the 19 variables,\nin all five configurations", title_fontsize=5.5)
    a.set_title(t, loc="left")
save(fig, f"{SUPF}/FigS6_robust_cols")
