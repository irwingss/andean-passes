"""52_figuras.py - Figures of the manuscript and the supplement for the census and the present connectivity.
Reads the census (data/interim/censo), the connectivity run (data/interim/conectividad/<run>, with the group summaries
in grupos/) and the tables of 51_sintesis_conectividad.py. Writes PNG (600 dpi) and PDF to outputs/figures (Fig1 to
Fig4) and outputs/supp (FigS1 to FigS7).
python -P scripts/52_figuras.py principal
Experiment, not part of the article: any second argument writes Figs 1 to 3 to
outputs/experimento_top_continentales/figuras, without the control passes and the barriers of Hazzi et al. (2018), and
marks the bridges for lowland climates and the barriers for montane climates of exp_top_continentales.py: both lists in Fig. 1, each list in its own direction
in Figs 2 and 3.
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
# One meaning per colour across all figures: red and blue are reserved for concentrated and diverted current; the
# rank colours appear only in Fig. 1; the two directions are told apart by panel, or by line style in Fig. 3c.
RANK_COL = {"continental": "#7B3294", "terminal": "#008837", "tributary": "#8C8C8C"}      # plus size
GROUP_COL = {"temperature": "#D55E00", "precipitation": "#0072B2"}        # warm for temperature, blue for water
CONC, DIV = "#B2182B", "#2166AC"
CUR = LinearSegmentedColormap.from_list("cur", [DIV, "#92C5DE", "#F7F7F7", "#F4A582", CONC])
NORM = LogNorm(0.2, 5)
HALO = [matplotlib.patheffects.withStroke(linewidth=1.6, foreground="w")]
# offsets (points) of the numbers of Table 2 where neighbouring passes would overlap; the rest go up and right
NUM_OFF = {3: (3, -6), 4: (2, 3), 6: (-7, 3), 7: (-8, -4), 13: (4, -1), 14: (-9, -4), 15: (4, -6), 16: (4, 2),
           18: (-9, 2), 19: (4, 1)}
GROUPS = {"temperature, mean and extremes": ["bio01", "bio05", "bio06", "bio08", "bio09", "bio10", "bio11"],
          "temperature, range and seasonality": ["bio02", "bio03", "bio04", "bio07"],
          "precipitation": [f"bio{v}" for v in range(12, 20)]}
GROUP_LS = {"temperature, mean and extremes": ("temperature", "-"),      # hue is the kind of variable,
            "temperature, range and seasonality": ("temperature", (0, (3, 1.3))),   # line style the sub-group
            "precipitation": ("precipitation", "-")}
READ_LS = {"longitudinal": dict(ls="-", mfc=None), "transverse": dict(ls="--", mfc="w")}
HB = (-80.6, -77.4, -7.5, -4.0)                    # box of the Huancabamba Depression in Fig. 1 and Fig. S7 (W, E, S, N)
# The six barriers of congruent phylogenetic breaks of endemic Andean birds (Hazzi et al. 2018, Figs 1 and 2),
# numbered as there. Shown for orientation only; they take no part in the analysis. Digitized from their Fig. 1 after
# fitting a Mercator transform to its rivers against HydroRIVERS (lon, lat; about 0.3 degrees of error).
_hz = pd.read_csv("data/refs/hazzi_2018_barreras.csv")
HAZZI = {k: (g.name.iat[0], list(zip(g.lon, g.lat))) for k, g in _hz.groupby("no")}
ROMAN = {1: "I", 2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI"}         # control passes already use 1 to 25
# zones of Fig. S7, boxed in Figs 1 and 2
ZOOM = (("Colombia and Táchira", (-77.8, -70.5, 0.6, 9.2)), ("Huancabamba Depression", HB),
        ("Southern Peru", (-75.0, -66.5, -16.5, -11.3)))

def panel(ax, letter, text=""):
    """Panel letter in bold and an optional short label in normal weight, above the top left corner."""
    ax.annotate(letter, (0, 1), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points",
                fontsize=8, fontweight="bold", va="bottom", ha="left")
    if text:
        ax.annotate(text, (0, 1), xycoords="axes fraction", xytext=(9, 3), textcoords="offset points",
                    fontsize=6.5, va="bottom", ha="left")


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


# country outlines (GADM 4.1 dissolved, scripts/exp_paises.R): orientation only, they enter no analysis. Only the six
# Andean countries are in data/raw, so Argentina and Brazil show as the far side of a border, without their own outline
PAIS = [shapely.from_wkb(w) for w in pyogrio.raw.read("data/interim/M1_paises.gpkg")[2]]
PAIS_M = [shapely.transform(g, lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1]))) for g in PAIS]


def countries(ax, deg=False, lw=0.25, zorder=4.5):
    """Thin country boundaries, in degrees or in the map projection."""
    for p in (PAIS if deg else PAIS_M):
        for g in getattr(p, "geoms", [p]):
            ax.plot(*g.exterior.xy, lw=lw, c="0.45", alpha=0.85, zorder=zorder)


def frame(ax, labels=True, names=False, passes=False):
    """Country outlines, GMBA outline of the tropics, graticule every 5 degrees and the control passes."""
    countries(ax)
    for g in getattr(band_m, "geoms", [band_m]):
        ax.plot(*g.exterior.xy, lw=0.35, c="k", alpha=0.6)
    for la in range(-20, 11, 5):
        lo = np.linspace(-85, -55, 200)
        x, y = to_m.transform(lo, np.full_like(lo, la))
        ax.plot(x, y, lw=0.25, c="0.6", ls=":", zorder=0.5)
        if labels:
            ax.text(EXT[0] + 15e3, np.interp(EXT[0] + 15e3, x, y), f"{abs(la)}°{'S' if la < 0 else 'N' if la else ''}",
                    fontsize=5.5, color="0.35", va="bottom", path_effects=HALO)
    for lo in range(-80, -59, 5):
        la = np.linspace(-30, 15, 200)
        x, y = to_m.transform(np.full_like(la, lo), la)
        ax.plot(x, y, lw=0.25, c="0.6", ls=":", zorder=0.5)
        if labels and EXT[0] < np.interp(EXT[2] + 15e3, y, x) < EXT[1] - 60e3:
            ax.text(np.interp(EXT[2] + 15e3, y, x), EXT[2] + 15e3, f"{abs(lo)}°W", fontsize=5.5, color="0.35",
                    ha="center", path_effects=HALO)
    if passes:                                                    # control passes: supplement only
        ax.scatter(ctl_x, ctl_y, s=9, marker="^", fc="w", ec="k", lw=0.5, zorder=5)
    if passes and names:
        for x, y, n in zip(ctl_x, ctl_y, range(1, len(ctl) + 1)):     # numbers of Table 2
            ax.annotate(str(n), (x, y), fontsize=4.5, xytext=NUM_OFF.get(n, (3, 1)), textcoords="offset points",
                        zorder=6, path_effects=HALO)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")


def box_m(b):
    """Outline of a (W, E, S, N) box in map units, through its four corners."""
    w, e, s, n = b
    return to_m.transform([w, e, e, w, w], [s, s, n, n, s])


def zones(ax):
    """The boxes of the zones of Fig. S7."""
    for _, b in ZOOM:
        ax.plot(*box_m(b), c="k", lw=0.6, zorder=6)


def hazzi(ax, lw=1.8, fs=5.5):
    """The six barriers of Hazzi et al. (2018), only on Fig. 3: the panel that measures what they measure."""
    for k, (_, pts) in HAZZI.items():
        x, y = to_m.transform(*zip(*pts))
        ax.plot(x, y, c="w", lw=lw, solid_capstyle="butt", zorder=6)
        ax.plot(x, y, c="k", lw=lw, ls=(0, (1.2, 1.2)), dash_capstyle="butt", zorder=6.1)
        ax.annotate(ROMAN[k], (x[-1], y[-1]), xytext=(2, 1), textcoords="offset points", fontsize=fs,
                    fontweight="bold", zorder=7, path_effects=HALO)


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
db = pd.read_csv("outputs/database/APD_passes.csv.gz", dtype={"trans_robust": str, "longi_robust": str})
tr = db[db.tropics == 1]
tr_x, tr_y = to_m.transform(tr.lon.values, tr.lat.values)
CONT = (tr["rank"] == "continental").to_numpy()                   # the continental divide, as a backbone

# ==== Fig. 1: the census ============================================================================================
fig = plt.figure(figsize=(183 / 25.4, 168 / 25.4))
ax = fig.add_axes([0.0, 0.02, 0.46, 0.94])
relief(ax, light=True)
d = tr.depth_m >= 100
for rk, sz in (("tributary", 0.25), ("terminal", 0.8), ("continental", 1.6)):
    m = (d & (tr["rank"] == rk)).to_numpy()
    ax.scatter(tr_x[m], tr_y[m], s=sz, c=RANK_COL[rk], lw=0, rasterized=True, zorder=2,
               label=f"{rk} ({m.sum():,})")
frame(ax)
ax.legend(loc="center right", fontsize=5.5, markerscale=4, frameon=False, title="Passes with a depth\nof 100 m or more",
          title_fontsize=5.5, alignment="left")
panel(ax, "a", "Passes of the tropical Andes by rank")
w0, w1, z0, z1 = HB
zones(ax)
# b: divide network and every pass around the Huancabamba Depression, from the 1 arc-second census
axb = fig.add_axes([0.53, 0.33, 0.45, 0.63])
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
countries(axb, deg=True, lw=0.3)
axb.set_xlim(w0, w1); axb.set_ylim(z0, z1)
axb.set_aspect(1 / np.cos(np.radians(5.75)))
axb.set_xticks([-80, -79, -78]); axb.set_xticklabels(["80°W", "79°W", "78°W"])
axb.set_yticks([-7, -6, -5, -4]); axb.set_yticklabels(["7°S", "6°S", "5°S", "4°S"])
axb.legend(loc="lower left", fontsize=5.5, markerscale=1.5, frameon=True, framealpha=0.85, edgecolor="none",
           title="Passes, sized by depth", title_fontsize=5.5, alignment="left")
panel(axb, "b", "Divide network and passes around the Huancabamba Depression")
# c: height of the tropical passes by rank (every pass with a depth of 30 m or more)
axc = fig.add_axes([0.59, 0.06, 0.39, 0.19])
bins = np.arange(0, 6301, 100)
for rk in ("tributary", "terminal", "continental"):
    m = tr[tr["rank"] == rk].z_m
    axc.hist(m, bins=bins, density=True, histtype="step", color=RANK_COL[rk], lw=0.9,
             label=f"{rk} ({len(m):,}; median {m.median():,.0f} m)")
    axc.axvline(m.median(), color=RANK_COL[rk], lw=0.6, ls=":")
axc.set_xlim(0, 6300); axc.set_yticks([])
axc.set_xlabel("Height of the pass (m)"); axc.set_ylabel("Density")
axc.spines[["top", "right"]].set_visible(False)
axc.legend(loc="upper left", fontsize=5.5, frameon=False, handlelength=1.2)
panel(axc, "c", "Height of the tropical passes with a depth of 30 m or more, by rank")
save(fig, f"{FIG}/Fig1_census")

# ==== Fig. 2: normalized current summarized by group of variables, both directions ================================
# four maps side by side: transverse with temperature and precipitation, then longitudinal; the zones of Fig. S7 boxed
GRP_LAB = (("temp", "temperature (BIO1 to BIO11)"), ("prec", "precipitation (BIO12 to BIO19)"))
READS = (("baja", "Transverse"), ("alta", "Longitudinal"))


def group(g, l, tag="s200_env", step=1):
    """Median normalized current over the variables of a group (51_sintesis_conectividad.py), nan where undefined."""
    with rasterio.open(f"{D}/grupos/{g}_{l}_{tag}.tif") as s:
        return np.ma.masked_invalid(s.read(1, window=win))[::step, ::step]


W, g, gt, bot = 183.0, 2.5, 6.0, 20.0
mw = (W - 3 * g) / 4
mh = mw * (EXT[3] - EXT[2]) / (EXT[1] - EXT[0])
H = bot + mh + gt
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
F = lambda x, y, w, h: fig.add_axes([x / W, y / H, w / W, h / H])
GZ = {(l, gg): group(gg, l) for l, _ in READS for gg, _ in GRP_LAB}
for i, ((l, rt), (gg, gl)) in enumerate([(r, q) for r in READS for q in GRP_LAB]):
    a = F(i * (mw + g), bot, mw, mh)
    relief(a, light=True)
    im = a.imshow(GZ[(l, gg)], extent=EXT, cmap=CUR, norm=NORM, interpolation="nearest", rasterized=True, alpha=0.92)
    frame(a, labels=i == 0)
    zones(a)
    panel(a, "abcd"[i], f"{rt}, {gl.split(' (')[0]}")
cax = F(0.3 * W, 8.0, 0.4 * W, 2.6)
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("Median normalized current over the variables of the group (log scale)")
cb.minorticks_off(); cb.set_ticks([0.2, 0.5, 1, 2, 5]); cb.set_ticklabels(["0.2", "0.5", "1", "2", "5"])
cax.annotate("diverted", (0, 0.5), xycoords="axes fraction", xytext=(-4, 0), textcoords="offset points",
             ha="right", va="center", fontsize=5.5, color=DIV)
cax.annotate("concentrated", (1, 0.5), xycoords="axes fraction", xytext=(4, 0), textcoords="offset points",
             ha="left", va="center", fontsize=5.5, color=CONC)
save(fig, f"{FIG}/Fig2_current_groups")

# ==== Fig. S7: the three zones of Fig. 2, both groups and both directions ===========================================
ZB = [shapely.Polygon(np.column_stack(box_m(b))).bounds for _, b in ZOOM]
W, g, gt = 183.0, 2.5, 5.0
iw = (W - 3 * g) / 4
hs = [iw * (y1 - y0) / (x1 - x0) for x0, y0, x1, y1 in ZB]
H = sum(hs) + 3 * gt + 16.0
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
F = lambda x, y, w, h: fig.add_axes([x / W, y / H, w / W, h / H])
m300 = (tr.depth_m >= 300).to_numpy()
ytop = H - gt
for r, ((name, _), (x0, y0, x1, y1), h) in enumerate(zip(ZOOM, ZB, hs)):
    for c, ((l, rt), (gg, gl)) in enumerate([(rr, q) for rr in READS for q in GRP_LAB]):
        b = F(c * (iw + g), ytop - h, iw, h)
        relief(b, light=True)
        im = b.imshow(GZ[(l, gg)], extent=EXT, cmap=CUR, norm=NORM, interpolation="nearest", rasterized=True,
                      alpha=0.92)
        b.scatter(tr_x[m300], tr_y[m300], s=0.5, c="k", lw=0, alpha=0.7, zorder=4, rasterized=True)
        countries(b, lw=0.3)
        b.set_xlim(x0, x1); b.set_ylim(y0, y1); b.set_xticks([]); b.set_yticks([]); b.set_aspect("equal")
        b.text(0.03, 0.97, f"{rt.lower()}, {gl.split(' (')[0]}", transform=b.transAxes, fontsize=5, ha="left",
               va="top", bbox=dict(fc="w", ec="none", pad=0.6, alpha=0.85), zorder=7)
        panel(b, "abcdefghijkl"[4 * r + c], name if c == 0 else "")
    ytop -= h + gt
cax = F(0.3 * W, 5.0, 0.4 * W, 2.6)
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("Median normalized current over the variables of the group (log scale)")
cb.minorticks_off(); cb.set_ticks([0.2, 0.5, 1, 2, 5]); cb.set_ticklabels(["0.2", "0.5", "1", "2", "5"])
save(fig, f"{SUPF}/FigS7_zones")

# ==== Fig. 3: the barrier as a property of the zone, and both directions of connectivity along the divide =======================
pr = pd.read_csv("outputs/supp/tables/S8_tile_conductance_profiles.csv")
dv = pd.read_csv("outputs/supp/tables/S17_divide_tile_profiles.csv")
sp = pd.read_csv("outputs/supp/tables/S10_spearman_variables.csv")
s16 = pd.read_csv("outputs/supp/tables/S16_by_height.csv")
s19 = pd.read_csv("outputs/supp/tables/S19_profile_bands.csv")
BAND_COL = "0.55"                                       # grey bands of the profiles, numbered as in Table S19


def bands_map(ax, b):
    """A band of Table S19 drawn in black over the stretch of the continental divide it spans, with its number."""
    m = CONT & (tr.lat.to_numpy() <= b.north) & (tr.lat.to_numpy() >= b.south)
    ax.scatter(tr_x[m], tr_y[m], s=1.2, c="k", lw=0, rasterized=True, zorder=6)
    top = b.no in (10,)                         # numbers that would cover a graticule label sit at the north end
    ax.annotate(str(b.no), (tr_x[m].min(), tr_y[m].max() if top else np.median(tr_y[m])), xytext=(-3, 0),
                textcoords="offset points", ha="right", va="center", fontsize=5.5, fontweight="bold", zorder=7,
                path_effects=HALO)


fig = plt.figure(figsize=(183 / 25.4, 128 / 25.4))
mw = 0.265
for i, (lect, t) in enumerate((("alta", "Montane climates, longitudinal"),
                               ("baja", "Lowland climates, transverse"))):
    a = fig.add_axes([i * (mw + 0.015), 0.12, mw, 0.82])
    with rasterio.open(f"{D}/zona_temperatura_{lect}.tif") as s:
        z = s.read(1, window=win, masked=True)
    im = a.imshow(z, extent=EXT, cmap="magma", vmin=0, vmax=1, interpolation="nearest", rasterized=True, zorder=1)
    frame(a, labels=i == 0)
    a.scatter(tr_x[CONT], tr_y[CONT], s=0.1, c="w", lw=0, rasterized=True, zorder=2)     # the continental divide
    hazzi(a)
    for _, b in s19[s19.reading == ("longitudinal" if lect == "alta" else "transverse")].iterrows():
        bands_map(a, b)
    panel(a, "ab"[i], t)
cb = fig.colorbar(im, cax=fig.add_axes([0.10, 0.055, 0.32, 0.018]), orientation="horizontal")
cb.set_label("Zone conductance: temperature median per tile, mean over the tiles of the cell")
cb.ax.annotate("barrier", (0, 0.5), xycoords="axes fraction", xytext=(-4, 0), textcoords="offset points", ha="right",
               va="center", fontsize=5.5)
cb.ax.annotate("open", (1, 0.5), xycoords="axes fraction", xytext=(4, 0), textcoords="offset points", ha="left",
               va="center", fontsize=5.5)
for i, rd in enumerate(("transverse", "longitudinal")):
    a = fig.add_axes([0.64 + i * 0.19, 0.12, 0.155, 0.82])
    for _, b in s19[s19.reading == rd].iterrows():               # bands of Table S19, numbered at the right
        a.axhspan(b.south, b.north, color=BAND_COL, alpha=0.35, lw=0, zorder=0.5)
        a.text(1.09, (b.north + b.south) / 2, str(b.no), fontsize=5.5, fontweight="bold", ha="left", va="center",
               clip_on=False)
    for g in ("temperature", "precipitation"):
        q = dv[(dv.cfg == "main") & (dv.reading == rd) & (dv.group == g)].set_index("band").sort_index()
        a.plot(q["median"], q.index, c=GROUP_COL[g], lw=1.2, label=g, zorder=3)
    a.axvline(1, c="0.5", lw=0.5, ls=":")
    a.set_xlim(0, 1.05); a.set_ylim(LS, LN); a.set_yticks(range(-20, 11, 5))
    a.set_yticklabels([f"{abs(v)}°{'S' if v < 0 else 'N' if v else ''}" for v in range(-20, 11, 5)] if i == 0 else [])
    a.set_xlabel("Relative conductance")
    a.spines[["top", "right"]].set_visible(False)
    panel(a, "cd"[i], rd.capitalize())
    if i == 0:
        a.set_ylabel("Latitude of the tile centre")
        a.legend(fontsize=5.5, frameon=False, loc="lower left", borderaxespad=0.2,
                 title="Tiles that cross the\ndivide; median of the\ngroup's variables; grey\nbands as in Table S19",
                 title_fontsize=5, alignment="left")
save(fig, f"{FIG}/Fig3_zone_conductance")

# ==== Fig. 4: the role of a pass against its height, and agreement between variables ==============================
fig = plt.figure(figsize=(183 / 25.4, 75 / 25.4))
a = fig.add_axes([0.07, 0.15, 0.40, 0.76])
for res, col, key in (("concentrates", CONC, "robust_above_1"), ("diverts", DIV, "robust_below_1")):
    for rd in ("longitudinal", "transverse"):
        x = s16[s16.reading == rd]
        a.plot(x.z_band + 250, x[key], c=col, lw=1.0, marker="o", ms=2.2, mec=col, **READ_LS[rd],
               label=f"{rd}, {res}")
a.axvline(2000, c="0.5", lw=0.5, ls=":")
a.text(2060, 0.98, "montane limit\nnear 2000 m", fontsize=5, va="top", color="0.35")
a.set_xlim(0, 6000); a.set_ylim(0, 1)
a.set_xlabel("Height of the pass (m), 500 m bands"); a.set_ylabel("Share of cells with passes")
a.spines[["top", "right"]].set_visible(False)
a.legend(loc="center", bbox_to_anchor=(0.66, 0.6), fontsize=5.5, frameon=False, handlelength=2.2)
panel(a, "a", "Passes that concentrate or divert current in all five configurations")
a = fig.add_axes([0.60, 0.15, 0.30, 0.76])
ORD = [v for vs in GROUPS.values() for v in vs]                 # variables in the three groups of Fig. 3
mt = sp[sp.reading == "transverse"].pivot(index="var_a", columns="var_b", values="rho_passes").loc[ORD, ORD].values
ml = sp[sp.reading == "longitudinal"].pivot(index="var_a", columns="var_b", values="rho_passes").loc[ORD, ORD].values
low = np.tril(np.ones((19, 19), bool), -1)
mm = np.where(low, mt, np.where(low.T, ml, np.nan))
im = a.imshow(mm, cmap="PuOr", vmin=-1, vmax=1, interpolation="nearest")     # not the scale of the current
for e in np.cumsum([len(vs) for vs in GROUPS.values()])[:-1] - 0.5:
    a.axhline(e, c="k", lw=0.5); a.axvline(e, c="k", lw=0.5)
a.plot([-0.5, 18.5], [-0.5, 18.5], c="k", lw=0.4)
lab = [str(int(v[3:])) for v in ORD]
a.set_xticks(range(19)); a.set_yticks(range(19)); a.set_xticklabels(lab, fontsize=4.5); a.set_yticklabels(lab, fontsize=4.5)
a.tick_params(length=1.5, pad=1)
BOX = dict(fc="w", ec="none", pad=0.6, alpha=0.85)
a.text(0.02, 0.02, "transverse", transform=a.transAxes, fontsize=5.5, va="bottom", ha="left", bbox=BOX)
a.text(0.98, 0.98, "longitudinal", transform=a.transAxes, fontsize=5.5, va="top", ha="right", bbox=BOX)
a.set_xlabel("BIO variable")
panel(a, "b", "Agreement between variables at the passes")
cb = fig.colorbar(im, cax=fig.add_axes([0.92, 0.18, 0.012, 0.68]))
cb.set_label("Spearman correlation", fontsize=6); cb.set_ticks([-1, -0.5, 0, 0.5, 1])
save(fig, f"{FIG}/Fig4_height_and_variables")

# ==== Fig. S1: validation of the census =============================================================================
V = "data/interim/censo/validacion"
osm = pd.read_csv(f"{V}/O5_1_nodos_osm.csv")
idn = pd.read_csv(f"{V}/O5_2_identidad.csv").set_index("group").loc["all"]
ref = pd.read_csv(f"{V}/O5_3_alturas_referencia.csv")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw={"wspace": 0.3})
x = np.sort(osm.census_km.values)
ax[0].plot(x, np.arange(1, len(x) + 1) / len(x), c="k", lw=0.9, label=f"OSM pass nodes ({len(x):,})")
ax[0].scatter([0.5, 1.0], [idn.null_500m_mean, idn.null_1km_mean], fc="w", ec="k", lw=0.7, s=14, zorder=3,
              label="random divide cells, same core and unit (mean of 1000)")
ax[0].set_xlim(0, 3); ax[0].set_ylim(0, 1)
ax[0].set_xlabel("Distance to the nearest pass (km)"); ax[0].set_ylabel("Cumulative share of nodes")
ax[0].legend(fontsize=5.5, frameon=False, loc="center right"); panel(ax[0], "a", "Recovery of OpenStreetMap passes")
ax[0].spines[["top", "right"]].set_visible(False)
for zn, fc in (("tropical", "k"), ("temperate", "w")):
    r = ref[ref.zone == zn]
    ax[1].scatter(r.elev_ref_m, r.error_m, s=12, fc=fc, ec="k", lw=0.6, label=f"{zn} ({len(r)})")
ax[1].axhline(0, c="k", lw=0.5); ax[1].axhspan(-50, 50, color="0.9", lw=0, zorder=0)
ax[1].set_xlabel("Official height of the pass (m)"); ax[1].set_ylabel("Nearest census pass minus official height (m)")
ax[1].legend(fontsize=5.5, frameon=False); panel(ax[1], "b", "Height against official heights")
ax[1].spines[["top", "right"]].set_visible(False)
save(fig, f"{SUPF}/FigS1_validation")

# ==== Fig. S2 and S3: the 19 variables, each direction ================================================================
for k, (l, t) in enumerate((("baja", "transverse"), ("alta", "longitudinal"))):
    fig, ax = plt.subplots(4, 5, figsize=(7.2, 9.4), gridspec_kw={"wspace": 0.02, "hspace": 0.12})
    for a, v in zip(ax.flat, VARS):
        im = a.imshow(current(f"{TS}/{v}_{l}_s200_env_normalizada.tif", 2), extent=EXT, cmap=CUR, norm=NORM,
                      interpolation="nearest", rasterized=True)
        for g in getattr(band_m, "geoms", [band_m]):
            a.plot(*g.exterior.xy, lw=0.25, c="k", alpha=0.6)
        countries(a, lw=0.15)
        a.set_xlim(EXT[0], EXT[1]); a.set_ylim(EXT[2], EXT[3]); a.set_xticks([]); a.set_yticks([])
        a.set_title(v.upper().replace("BIO0", "BIO"), fontsize=6)
    ax.flat[-1].axis("off")
    cb = fig.colorbar(im, ax=ax.flat[-1], orientation="horizontal", fraction=0.5, aspect=12)
    cb.set_label(f"Normalized current,\n{t} connectivity", fontsize=6)
    cb.minorticks_off(); cb.set_ticks([0.2, 0.5, 1, 2, 5]); cb.set_ticklabels(["0.2", "0.5", "1", "2", "5"])
    save(fig, f"{SUPF}/FigS{2 + k}_variables_{t}")

# ==== Fig. S4: the group summaries under the four sensitivity configurations, one figure per group ==================
conf = (("s200_env_q0.05", "q 0.05, 200 km"), ("s200_env_q0.2", "q 0.2, 200 km"), ("s100_env", "q 0.1, 100 km"),
        ("s400_env", "q 0.1, 400 km"))
for gg, gl in GRP_LAB:
    fig, ax = plt.subplots(2, 4, figsize=(7.2, 6.4), gridspec_kw={"wspace": 0.02, "hspace": 0.08})
    for i, (l, t) in enumerate(READS):
        for j, (f, lab) in enumerate(conf):
            a = ax[i, j]
            im = a.imshow(group(gg, l, f, 2), extent=EXT, cmap=CUR, norm=NORM, interpolation="nearest", rasterized=True)
            frame(a, labels=False, names=True, passes=True)
            panel(a, "abcdefgh"[i * 4 + j], f"{t}, {lab}")
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.4, pad=0.02, aspect=30)
    cb.set_label(f"Median normalized current, {gl}"); cb.minorticks_off(); cb.set_ticks([0.2, 0.5, 1, 2, 5])
    cb.set_ticklabels(["0.2", "0.5", "1", "2", "5"])
    save(fig, f"{SUPF}/FigS4_sensitivity_maps_{gg}")

# ==== Fig. S5: profiles of the tiles that cross the divide by group under the five configurations ====================
fig, ax = plt.subplots(2, 2, figsize=(5.0, 9.0), sharey=True, gridspec_kw={"wspace": 0.08, "hspace": 0.18})
# margins solid, tile sizes dashed; dark for the smaller value, light for the larger (no colour of the current)
cc = {"main": ("k", "-"), "q 0.05": ("0.35", "-"), "q 0.2": ("0.7", "-"), "100 km": ("0.35", "--"),
      "400 km": ("0.7", "--")}
for i, gr in enumerate(("temperature", "precipitation")):
    for k, rd in enumerate(("transverse", "longitudinal")):
        a = ax[i, k]
        a.axvline(1, c="0.5", lw=0.5, ls=":", zorder=0.5)
        for c, (col, ls) in cc.items():
            p = dv[(dv.cfg == c) & (dv.reading == rd) & (dv.group == gr)].sort_values("band")
            a.plot(p["median"], p.band, c=col, ls=ls, lw=1.2 if c == "main" else 0.8, label=c,
                   zorder=3 if c == "main" else 2)
        a.set_xlim(0, 1.05); a.set_ylim(LS, LN)
        a.set_xlabel(f"Relative conductance, {gr} median")
        a.spines[["top", "right"]].set_visible(False)
        panel(a, "abcd"[2 * i + k], f"{rd.capitalize()} connectivity, {gr}")
    ax[i, 0].set_ylabel("Latitude (°)")
ax[0, 0].legend(fontsize=5.5, frameon=False, loc="lower left")
save(fig, f"{SUPF}/FigS5_sensitivity_profiles")

# ==== Fig. S6: passes by the median over the 19 variables and its robustness ==========================================
fig, ax = plt.subplots(1, 2, figsize=(7.2, 6.4), gridspec_kw={"wspace": 0.03})
u = tr.sort_values("depth_m", ascending=False).drop_duplicates("cell_1km")      # deepest pass per cell, as in 51
u = u[u.gmba.notna()]                                                           # as in Table S13
ux, uy = to_m.transform(u.lon.values, u.lat.values)
RC = {"depends": "#D9D9D9", "below 1": DIV, "above 1": CONC}
for a, p, t in zip(ax, ("trans", "longi"), ("a", "b")):
    relief(a, light=True)
    for k in ("depends", "below 1", "above 1"):
        m = (u[f"{p}_robust"] == k).to_numpy()
        a.scatter(ux[m], uy[m], s=0.15, c=RC[k], lw=0, rasterized=True, label=f"{k} ({m.sum():,} cells)")
    frame(a, labels=a is ax[0], passes=True)
    a.legend(loc="center right", fontsize=5.5, markerscale=8, frameon=False, alignment="left",
             title="Median over the 19 variables,\nin all five configurations", title_fontsize=5.5)
    panel(a, t, "Transverse connectivity" if t == "a" else "Longitudinal connectivity")
save(fig, f"{SUPF}/FigS6_robust_passes")
