"""42_revision_bloque.py - Review sheets of one census block for the author (PRD R2), in data/interim/censo/:
revision_<tag>.pdf with (1) relief, divide network, passes by rank and depth, and each control pass with its nearest pass,
measured with relief only (D14, D16); (2) the real sinks (O2.2) and how many lie on the Amazon side; (3) 12 random
passes with depth >= 100 m on relief clips, to check none sits on a valley floor; (4) height percentile within 500 m of
sink-rim and crest passes, and the comparison with the neighbour block shifted half a degree south, if it exists.
Tables: controles_<tag>.csv, sumideros_revision_<tag>.csv, muestra_passes_<tag>.csv.
python -P scripts/42_revision_bloque.py <tag>     after 40 and 41.
"""
import csv
import importlib.util
import math
import os
import sys

import matplotlib
import numpy as np
import pyogrio
import rasterio
import shapely
import yaml
from pyproj import Geod

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

cfg = yaml.safe_load(open("config.yml", encoding="utf-8"))
OUT = os.path.join(cfg["paths"]["interim"], "censo")
tag = sys.argv[1]
GEOD = Geod(ellps="WGS84")

src = rasterio.open(os.path.join(OUT, "tmp", f"raw_{tag}.tif"))
w, n = src.transform.c, src.transform.f
res = src.transform.a
z = src.read(1)
dv = rasterio.open(os.path.join(OUT, f"divisorias_{tag}.tif")).read(1)
passes = list(csv.DictReader(open(os.path.join(OUT, f"passes_{tag}_etiquetas.csv"), encoding="utf-8")))
for c in passes:
    for k in ("lat", "lon", "z_m", "depth_m"):
        c[k] = float(c[k])
lat = np.array([c["lat"] for c in passes]); lon = np.array([c["lon"] for c in passes])
depth = np.array([c["depth_m"] for c in passes]); zc = np.array([c["z_m"] for c in passes])
core = np.array([float(c["core"]) == 1 for c in passes])
rank = np.array([c["rank"] for c in passes]); gate = np.array([float(c["gateway_window"]) == 1 for c in passes])
s_lat, e_lon = n - z.shape[0] * res, w + z.shape[1] * res
ctrl = [r for r in csv.DictReader(open("data/refs/pasos_control.csv", encoding="utf-8"))
        if s_lat < float(r["lat"]) < n and w < float(r["lon"]) < e_lon]
spec = importlib.util.spec_from_file_location("censo40", "scripts/40_censo_bloque.py")
c40 = importlib.util.module_from_spec(spec); sys.modules["censo40"] = c40; spec.loader.exec_module(c40)
AREA = c40.area_mask()

st = cfg["style"]["elev_ramp"]
RAMP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "elev", [((m - st["m"][0]) / (st["m"][-1] - st["m"][0]), c) for m, c in zip(st["m"], st["col"])])
NORM = matplotlib.colors.Normalize(st["m"][0], st["m"][-1])


def shade(a, cellx, celly):
    a = np.where(a < -1000, np.nan, a).astype(float)
    gy, gx = np.gradient(a, celly, cellx)
    sl, asp = np.arctan(np.hypot(gx, gy)), np.arctan2(-gx, gy)
    az, alt = math.radians(315), math.radians(45)
    h = np.sin(alt) * np.cos(sl) + np.cos(alt) * np.sin(sl) * np.cos(az - asp)
    return np.where(np.isnan(a), np.nan, h)


def panel(ax, i0, i1, j0, j1, f, pts=True, size=1.0, contour=None):
    sub = z[i0:i1:f, j0:j1:f]
    cy = 30.9 * f
    cx = cy * math.cos(math.radians(n - (i0 + i1) / 2 * res))
    ext = (w + j0 * res, w + j1 * res, n - i1 * res, n - i0 * res)
    ax.set_facecolor("#cfe3ef")                                   # sea
    ax.imshow(np.ma.masked_less(sub, -1000), extent=ext, cmap=RAMP, norm=NORM, interpolation="nearest")
    ax.imshow(shade(sub, cx, cy), extent=ext, cmap="gray", alpha=0.35, interpolation="nearest")
    d = dv[i0:i1, j0:j1]
    if f > 1:
        hh, ww = d.shape[0] // f * f, d.shape[1] // f * f
        d = d[:hh, :ww].reshape(hh // f, f, ww // f, f).max(axis=(1, 3))
    ax.imshow(np.ma.masked_equal(d, 0), extent=ext, cmap="Greys", vmin=0, vmax=1, alpha=0.55, interpolation="nearest")
    if contour:
        zz = np.where(sub < -1000, np.nan, sub)
        lv = np.arange(np.floor(np.nanmin(zz) / contour) * contour, np.nanmax(zz) + contour, contour)
        ax.contour(np.linspace(ext[0], ext[1], sub.shape[1]), np.linspace(ext[3], ext[2], sub.shape[0]), zz, lv,
                   colors="k", linewidths=0.3, alpha=0.6)
    if pts:
        m = (lon > ext[0]) & (lon < ext[1]) & (lat > ext[2]) & (lat < ext[3])
        for r, col in (("tributary", "#bbbbbb"), ("terminal", "#3b8bd4"), ("continental", "#d62728")):
            k = m & (rank == r)
            ax.scatter(lon[k], lat[k], s=size * (4 + 8 * (depth[k] >= 100) + 16 * (depth[k] >= 300)), c=col,
                       edgecolors="k", linewidths=0.2, zorder=3, label=r)
        k = m & gate
        ax.scatter(lon[k], lat[k], s=90 * size, marker="*", c="yellow", edgecolors="k", linewidths=0.5, zorder=4,
                   label="gateway (window)")
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    return ext


def outline(ax):
    for g in shapely.get_parts(shapely.intersection(AREA, shapely.box(w, s_lat, e_lon, n))):
        for ring in [g.exterior, *g.interiors]:
            x, y = ring.xy
            ax.plot(x, y, color="#7b2d8e", lw=0.8)


def wcsv(name, hdr, rows):
    with open(os.path.join(OUT, f"{name}_{tag}.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f); wr.writerow(hdr); wr.writerows(rows)


pdf = PdfPages(os.path.join(OUT, f"revision_{tag}.pdf"))

# ---------------------------------------------------------------- page 1: block and control passes
rows = []
for r in ctrl:
    la, lo = float(r["lat"]), float(r["lon"])
    _, _, dist = GEOD.inv(np.full(lat.size, lo), np.full(lat.size, la), lon, lat)
    a = np.argmin(dist)
    b = np.argmin(np.where(depth >= 100, dist, np.inf))
    rows.append([r["name"], la, lo, f"{dist[a] / 1000:.2f}", f"{zc[a]:.0f}", f"{depth[a]:.0f}", rank[a],
                 passes[a]["area_sep_km2"], f"{dist[b] / 1000:.2f}", f"{zc[b]:.0f}", f"{depth[b]:.0f}", rank[b],
                 passes[b]["area_sep_km2"], passes[b]["width_m"]])
hdr = ["control", "lat", "lon", "near_km", "near_z", "near_depth", "near_rank", "near_area_sep_km2",
       "deep_km", "deep_z", "deep_depth", "deep_rank", "deep_area_sep_km2", "deep_width_m"]
wcsv("controles", hdr, rows)
for r in rows:
    print(dict(zip(hdr, r)))

nz = len(ctrl)
fig = plt.figure(figsize=(16, 9 + 4.5 * math.ceil(nz / 3)))
gs = fig.add_gridspec(2 + math.ceil(nz / 3), 3)
ax = fig.add_subplot(gs[:2, :])
panel(ax, 0, z.shape[0], 0, z.shape[1], 8)
outline(ax)
for v in (max(lat[core]), min(lat[core])):
    ax.axhline(v, color="k", lw=0.6, ls="--")
for r in ctrl:
    ax.plot(float(r["lon"]), float(r["lat"]), "s", mfc="none", mec="m", ms=9, mew=1.5)
    ax.annotate(r["name"], (float(r["lon"]), float(r["lat"])), xytext=(5, 5), textcoords="offset points",
                color="m", fontsize=9, weight="bold")
ax.legend(loc="lower left", fontsize=8)
ax.set_title(f"Block {tag}: {int(core.sum())} core passes >= 30 m deep inside the area (" +
             ", ".join(f"{r} {int(((rank == r) & core).sum())}" for r in ("continental", "terminal", "tributary")) +
             f"); {int((core & (depth >= 100)).sum())} >= 100 m; gateways in window {int((core & gate).sum())}.\n"
             "Dashed: core limits. Purple: area (GMBA Andes + 20 km). Grey: divide network at 0.1 km2.", fontsize=10)
for q, r in enumerate(ctrl):
    a = fig.add_subplot(gs[2 + q // 3, q % 3])
    la, lo = float(r["lat"]), float(r["lon"])
    i, j = int((n - la) / res), int((lo - w) / res)
    h = 330
    panel(a, max(i - h, 0), i + h, max(j - h, 0), j + h, 1, size=3, contour=50)
    a.plot(lo, la, "s", mfc="none", mec="m", ms=12, mew=2)
    row = rows[q]
    a.set_title(f"{r['name']} (OSM node, square): nearest pass {row[3]} km, {row[4]} m, depth {row[5]} m ({row[6]})\n"
                f"nearest pass with depth >= 100 m: {row[8]} km, {row[9]} m, depth {row[10]} m, width {row[13]} m",
                fontsize=8)
fig.tight_layout()
pdf.savefig(fig); fig.savefig(os.path.join(OUT, f"revision_{tag}.png"), dpi=110); plt.close(fig)

# ---------------------------------------------------------------- page 2: real sinks
# Real sinks (PRD O2.2): bottom in an endorheic HydroBASINS basin and 20 m or more deep. The Amazon side is the
# Atlantic domain of M1_domains (dissolved from the same HydroBASINS); with this rule it should hold almost none.
sk = list(csv.DictReader(open(os.path.join(OUT, f"sumideros_{tag}.csv"), encoding="utf-8")))


def col_(k):
    return np.array([float(x[k]) if x[k] != "" else np.nan for x in sk])


slat, slon, sarea, sdep = col_("lat"), col_("lon"), col_("contrib_km2"), col_("depth_m")
score, sin, sedge = col_("core") == 1, col_("in_area") == 1, col_("bottom_on_edge") == 1
m, _, gd, Fd = pyogrio.raw.read(os.path.join(cfg["paths"]["interim"], "M1_domains.gpkg"), columns=["dom"])
dom = dict(zip(Fd[0], shapely.from_wkb(gd)))
amaz = shapely.contains_xy(dom["E"], slon, slat)
hdr2 = ["lat", "lon", "bottom_z_m", "fill_depth_m", "contrib_km2", "rim_census_depth_m", "in_area", "amazon_side"]
rrows = [[sk[k]["lat"], sk[k]["lon"], sk[k]["bottom_z_m"], sk[k]["fill_depth_m"], f"{sarea[k]:.2f}", sk[k]["depth_m"],
          sk[k]["in_area"], int(amaz[k])] for k in np.flatnonzero(score & ~sedge)]
wcsv("sumideros_revision", hdr2, rrows)
c = score & ~sedge
txt = (f"Sinks left in the window after breaching: {len(sk)}; with the bottom on the window edge (artefact of the "
       f"window, never a pass): {int(sedge.sum())}. The rest are real sinks: bottom in an endorheic HydroBASINS "
       f"level-12 basin and 20 m or more deep.\n"
       f"Core of the block: {int(c.sum())} sinks, {int((c & sin).sum())} inside the area; on the Amazon side "
       f"{int((c & amaz).sum())}, {int((c & amaz & sin).sum())} inside the area. With a rim pass >= 30 m deep: "
       f"{int((c & (sdep >= 30)).sum())}, {int((c & sin & (sdep >= 30)).sum())} inside the area.")
print(txt)
fig = plt.figure(figsize=(16, 18))
gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1])
ax = fig.add_subplot(gs[0])
panel(ax, 0, z.shape[0], 0, z.shape[1], 8, pts=False)
outline(ax)
k = ~sedge
sc = ax.scatter(slon[k], slat[k], c=np.log10(np.maximum(sarea[k], 1e-3)), s=12 + 30 * (sdep[k] >= 30), cmap="plasma",
                vmin=-3, vmax=3, edgecolors="k", linewidths=0.4, zorder=3)
for gp in shapely.get_parts(shapely.intersection(dom["E"], shapely.box(w, s_lat, e_lon, n))):
    x, y = gp.exterior.xy
    ax.plot(x, y, color="blue", lw=0.8)
fig.colorbar(sc, ax=ax, shrink=0.6, label="contributing area, log10 km2 (large dots: rim pass >= 30 m deep)")
ax.set_title(f"Block {tag}: real sinks. Blue: limit of the Atlantic (Amazon) side. Purple: area.")
axt = fig.add_subplot(gs[1]); axt.axis("off")
axt.text(0, 0.95, txt, fontsize=9, va="top", wrap=True)
if rrows:
    tb = axt.table(cellText=rrows[:40], colLabels=hdr2, loc="lower center")
    tb.auto_set_font_size(False); tb.set_fontsize(7)
fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
# ---------------------------------------------------------------- page 3: 12 random deep passes
rng = np.random.default_rng(cfg["seed"])
pool = np.flatnonzero(core & shapely.contains_xy(AREA, lon, lat) & (depth >= 100))
pick = np.sort(rng.choice(pool, size=min(12, pool.size), replace=False))
wcsv("muestra_passes", ["lat", "lon", "z_m", "depth_m", "width_m", "rank", "ring", "flat", "depth_lower_bound"],
     [[passes[k][c] for c in ("lat", "lon", "z_m", "depth_m", "width_m", "rank", "ring", "flat", "depth_lower_bound")]
      for k in pick])
fig, axs = plt.subplots(4, 3, figsize=(16, 21))
for a, k in zip(axs.ravel(), pick):
    i, j = int((n - lat[k]) / res), int((lon[k] - w) / res)
    h = 90                                                  # 2.8 km either side
    panel(a, max(i - h, 0), i + h, max(j - h, 0), j + h, 1, pts=False, contour=20)
    a.plot(lon[k], lat[k], "o", mfc="red", mec="k", ms=8)
    lb = " (lower bound)" if float(passes[k]["depth_lower_bound"]) == 1 else ""
    a.set_title(f"{lat[k]:.4f}, {lon[k]:.4f}: {zc[k]:.0f} m, depth {depth[k]:.0f} m{lb}, width {passes[k]['width_m']} m\n"
                f"{rank[k]}{', closed-basin rim' if float(passes[k]['ring']) == 1 else ''}"
                f"{', flat' if float(passes[k]['flat']) == 1 else ''}. Contours 20 m, grey: divides.", fontsize=8)
for a in axs.ravel()[len(pick):]:
    a.axis("off")
fig.suptitle(f"Block {tag}: 12 random passes with depth >= 100 m (core, inside the area; {pool.size} in all; "
             f"seed {cfg['seed']})", fontsize=11)
fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
# ---------------------------------------------------------------- page 4: relief position and neighbour block
# Height percentile of each core pass among the FABDEM cells within 500 m: a pass on a valley floor sits low.
ina = core & shapely.contains_xy(AREA, lon, lat)
ring = np.array([float(c["ring"]) == 1 for c in passes])
ry = 500 / 30.92
pct = np.full(lat.size, np.nan)
for k in np.flatnonzero(ina):
    i, j = int((n - lat[k]) / res), int((lon[k] - w) / res)
    rx = ry / math.cos(math.radians(lat[k]))
    a, b = int(ry), int(rx)
    sub = z[max(i - a, 0):i + a + 1, max(j - b, 0):j + b + 1]
    yy, xx = np.ogrid[max(i - a, 0) - i:min(i + a + 1, z.shape[0]) - i, max(j - b, 0) - j:min(j + b + 1, z.shape[1]) - j]
    v = sub[((yy / ry) ** 2 + (xx / rx) ** 2 <= 1) & (sub > -1000)]
    pct[k] = 100 * np.mean(v < zc[k])
lines = []
for nm, kk in (("crest", ina & ~ring), ("sink rim", ina & ring)):
    q = np.percentile(pct[kk], [10, 25, 50, 75, 90]) if kk.any() else [np.nan] * 5
    lines.append(f"{nm} passes: {int(kk.sum())}; height percentile within 500 m, p10 {q[0]:.0f}, p25 {q[1]:.0f}, "
                 f"median {q[2]:.0f}, p75 {q[3]:.0f}, p90 {q[4]:.0f}")
# neighbour block shifted half a degree south: same passes and depths in the common core strip?
north = (-1 if tag[0] == "S" else 1) * float(tag[1:].replace("p", "."))
nb_tag = f"{'S' if north - 0.5 <= 0 else 'N'}{abs(north - 0.5):05.2f}".replace(".", "p")
fn = os.path.join(OUT, f"passes_{nb_tag}.csv")
if os.path.exists(fn):
    def strip(f):
        r = [x for x in csv.DictReader(open(f, encoding="utf-8"))
             if x["core"] == "1" and north - 1 < float(x["lat"]) <= north - 0.5]
        return {(x["lat"], x["lon"]): (float(x["depth_m"]), x["ring"] == "1") for x in r}
    A, B = strip(os.path.join(OUT, f"passes_{tag}.csv")), strip(fn)
    both = A.keys() & B.keys()
    same = [p for p in both if abs(A[p][0] - B[p][0]) < 0.05]
    diff = [p for p in both if p not in same]
    onlyA, onlyB = A.keys() - B.keys(), B.keys() - A.keys()
    rg = lambda ps, D: sum(D[p][1] for p in ps)
    lines.append(f"Neighbour {nb_tag}, strip {north - 0.5:.1f} to {north - 1:.1f} (passes of the area in each core): "
                 f"{len(A)} here, {len(B)} there; same position {len(both)}, of which same depth {len(same)} and "
                 f"different depth {len(diff)} ({rg(diff, A)} sink rim). Only here {len(onlyA)} ({rg(onlyA, A)} sink rim), "
                 f"only there {len(onlyB)} ({rg(onlyB, B)} sink rim).")
txt4 = "\n".join(lines)
print(txt4)
fig, (ax, axt) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={"height_ratios": [2, 1]})
ax.hist([pct[ina & ~ring], pct[ina & ring]], bins=20, range=(0, 100), density=True, label=["crest", "sink rim"])
ax.set_xlabel("height percentile of the pass among FABDEM cells within 500 m"); ax.legend()
ax.set_title(f"Block {tag}: relief position of core passes inside the area")
axt.axis("off"); axt.text(0, 0.9, txt4, fontsize=9, va="top", wrap=True)
fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
pdf.close()
print("->", os.path.join(OUT, f"revision_{tag}.pdf"))
