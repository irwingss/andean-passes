"""47_conectividad_insumos.py - Inputs of the climatic connectivity (PRD block C) on a 1 km equal-area grid.
Writes, under data/interim/conectividad/<run>/insumos: the 19 WorldClim 2.1 bioclimatic variables (bio12 to bio19 as
log(1 + P)), FABDEM averaged to 1 km, the two guild masks (cells warmer or colder than T2000(lat), C3) and, per guild,
variable and margin q, the envelope resistance (C4): r = 1 + d / (q * IQR), where d is how far the cell falls outside the
5 to 95 percentile envelope of the variable in the guild cells within 1 degree of its latitude, and IQR is that of the
variable over the domain (margenes.csv). Files res_envbaja_bioNN.tif and res_envalta_bioNN.tif for the first q; further
q values add the suffix _q<q> (e.g. res_envbaja_bio01_q0.05.tif).
python -P scripts/47_conectividad_insumos.py <earth_engine_project> <run> <lat_north> <lat_south> <vars> <q> <shapes>
  e.g. ee-irwingssaldu piloto -3 -9 "1 12" "0.1" env      or      ee-irwingssaldu principal 11 -23.5 all "0.1 0.05 0.2" env
The domain is the GMBA area (D12) inside the reported band plus a fixed 220 km buffer, so the lowlands of both sides
enter the circuit (C1).
FABDEM comes from Earth Engine only when insumos/elev.tif is missing or on another grid; otherwise elev.tif is read
locally and Earth Engine is not initialized (the project argument is then unused).
"""
import glob
import itertools
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pyogrio
import rasterio
import rasterio.features
import shapely
from pyproj import CRS, Transformer
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject

EP, RUN, LN, LS = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
VARS = list(range(1, 20)) if sys.argv[5] == "all" else [int(v) for v in sys.argv[5].split()]
QS = [float(v) for v in sys.argv[6].split()]
SHAPES = sys.argv[7].split()
assert SHAPES == ["env"], "only the envelope shape (env) is implemented"
OUT = f"data/interim/conectividad/{RUN}"
INS = f"{OUT}/insumos"
os.makedirs(INS, exist_ok=True)
CRS_ = CRS.from_proj4("+proj=laea +lat_0=-6 +lon_0=-72 +datum=WGS84 +units=m +no_defs")   # C1
RES = 1000.0
NODATA = -9999.0

# domain: GMBA area in the band plus 220 km, on the 1 km grid
gmba = shapely.from_wkb(pyogrio.raw.read("data/interim/M1_gmba_andes.gpkg")[2][0])
band = shapely.intersection(gmba, shapely.box(-90, LS, -55, LN))
to_m = Transformer.from_crs("EPSG:4326", CRS_, always_xy=True)
dom = shapely.buffer(shapely.transform(band, lambda c: np.column_stack(to_m.transform(c[:, 0], c[:, 1]))),
                     220 * 1000.0)
x0, y0, x1, y1 = [round(v / RES) * RES for v in dom.bounds]
W, H = int((x1 - x0) / RES), int((y1 - y0) / RES)
T = from_origin(x0, y1, RES, RES)
prof = dict(driver="GTiff", width=W, height=H, count=1, dtype="float32", crs=CRS_, transform=T, nodata=NODATA,
            compress="deflate", predictor=3, tiled=True)
inside = rasterio.features.geometry_mask([dom], (H, W), T, invert=True)
print(f"grid {W} x {H} km, {inside.sum():,} cells in the domain")


def save(name, a):
    with rasterio.open(f"{INS}/{name}.tif", "w", **prof) as d:
        d.write(np.where(np.isfinite(a), a, NODATA).astype(np.float32), 1)


# WorldClim 2.1, 30 seconds, bilinear to the grid (C1)
bio = {}
tiles = sorted(glob.glob("data/raw/clima_wc/climate/tiles_2.1/tile_*_wc2.1_30s_bio.tif"))
for v in sorted(set(VARS) | {1}):
    a = np.full((H, W), np.nan, np.float32)
    for f in tiles:
        with rasterio.open(f) as s:
            b = np.full((H, W), np.nan, np.float32)
            reproject(rasterio.band(s, v), b, dst_transform=T, dst_crs=CRS_, dst_nodata=np.nan,
                      resampling=Resampling.bilinear)
        a = np.where(np.isnan(a), b, a)
    bio[v] = np.log1p(np.maximum(a, 0)) if v >= 12 else a          # C2: precipitation as log(1 + P)

# FABDEM averaged to 1 km (D19): the local elev.tif if it has this grid, else Earth Engine in 200 by 200 km chunks
same = False
if os.path.exists(f"{INS}/elev.tif"):
    with rasterio.open(f"{INS}/elev.tif") as s:
        same = (s.width, s.height, s.transform, s.crs) == (W, H, T, CRS_)
        if same:
            elev = s.read(1, masked=True).filled(np.nan).astype(np.float32)
            print("elev.tif read locally")


def fabdem_ee():
    import ee
    ee.Initialize(project=EP)
    fab = ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
    fab = fab.mosaic().setDefaultProjection(fab.first().projection()).reduceResolution(ee.Reducer.mean(), maxPixels=4096)
    wkt = CRS_.to_wkt("WKT1_GDAL")         # Earth Engine reads WKT1 only
    elev = np.full((H, W), np.nan, np.float32)

    def chunk(rc):
        r, c = rc
        h, w = min(200, H - r), min(200, W - c)
        if not inside[r:r + h, c:c + w].any():
            return
        grid = {"dimensions": {"width": w, "height": h}, "crsWkt": wkt,
                "affineTransform": {"scaleX": RES, "shearX": 0, "translateX": x0 + c * RES, "shearY": 0,
                                    "scaleY": -RES, "translateY": y1 - r * RES}}
        for k in range(6):
            try:
                a = ee.data.computePixels({"expression": fab.unmask(-32768), "fileFormat": "NUMPY_NDARRAY", "grid": grid})
                break
            except Exception:
                if k == 5:
                    raise
                time.sleep(10 * (k + 1))
        a = a[a.dtype.names[0]].astype(np.float32)
        elev[r:r + h, c:c + w] = np.where(a <= -32000, np.nan, a)

    with ThreadPoolExecutor(8) as ex:
        list(ex.map(chunk, itertools.product(range(0, H, 200), range(0, W, 200))))
    return elev


if not same:
    elev = fabdem_ee()
valid = inside & np.isfinite(elev) & np.all([np.isfinite(bio[v]) for v in bio], axis=0)
save("elev", np.where(valid, elev, np.nan))
print(f"{valid.sum():,} valid cells")

# T2000(lat): median bio1 of cells at 1900 to 2100 m within 1 degree of latitude (C3, D9)
xs = x0 + (np.arange(W) + 0.5) * RES
ys = y1 - (np.arange(H) + 0.5) * RES
lat = Transformer.from_crs(CRS_, "EPSG:4326", always_xy=True).transform(*np.meshgrid(xs, ys))[1].astype(np.float32)
sel = valid & (elev >= 1900) & (elev <= 2100)
lc, tc = lat[sel], bio[1][sel]
bins = np.round(np.arange(np.floor(lat[valid].min() * 10), np.ceil(lat[valid].max() * 10) + 1) / 10, 1)
t2000 = np.full(bins.size, np.nan)
for i, b in enumerate(bins):
    m = np.abs(lc - b) <= 1
    if m.sum() >= 200:
        t2000[i] = np.median(tc[m])
ok = np.flatnonzero(np.isfinite(t2000))
t2000 = t2000[ok[np.abs(ok[:, None] - np.arange(bins.size)).argmin(0)]]    # nearest valid window
np.savetxt(f"{INS}/t2000.csv", np.column_stack([bins, t2000]), delimiter=",", header="lat,t2000", comments="", fmt="%.2f")
tl = np.interp(lat, bins, t2000)
fb, fa = valid & (bio[1] > tl), valid & (bio[1] <= tl)
save("fuente_baja", np.where(valid, fb.astype(np.float32), np.nan))      # lowland guild, transverse reading
save("fuente_alta", np.where(valid, fa.astype(np.float32), np.nan))      # montane guild, longitudinal reading

# variables, envelope resistances and margins (C2 to C4)
order = {k: np.argsort(lat[g]) for k, g in (("baja", fb), ("alta", fa))}
lat_s = {k: lat[g][order[k]] for k, g in (("baja", fb), ("alta", fa))}


def envelope(a, k, g):
    """5 and 95 percentiles of a in guild k per 0.1 degree band, window of 1 degree, at least 200 cells, else the
    nearest valid band; interpolated to every cell by its latitude."""
    x = a[g][order[k]]
    lo, hi = np.searchsorted(lat_s[k], bins - 1, "left"), np.searchsorted(lat_s[k], bins + 1, "right")
    e = np.full((bins.size, 2), np.nan)
    for i in range(bins.size):
        if hi[i] - lo[i] >= 200:
            e[i] = np.percentile(x[lo[i]:hi[i]], [5, 95])
    ok = np.flatnonzero(np.isfinite(e[:, 0]))
    e = e[ok[np.abs(ok[:, None] - np.arange(bins.size)).argmin(0)]]
    return np.interp(lat, bins, e[:, 0]), np.interp(lat, bins, e[:, 1]), e


marg, envs = [], []
for v in VARS:
    a = np.where(valid, bio[v], np.nan)
    save(f"bio{v:02d}", a)
    q75, q25 = np.percentile(a[valid], [75, 25])
    marg.append((v, q75 - q25))
    for k, g in (("baja", fb), ("alta", fa)):
        p5, p95, e = envelope(a, k, g)
        envs += [(k, v, b, l, h) for b, (l, h) in zip(bins, e)]
        d = np.maximum(np.maximum(p5 - a, a - p95), 0)
        for j, q in enumerate(QS):
            save(f"res_env{k}_bio{v:02d}" + ("" if j == 0 else f"_q{q:g}"),
                 np.where(valid, 1 + d / (q * (q75 - q25)), np.nan))
    print(f"bio{v:02d} done", flush=True)
np.savetxt(f"{INS}/margenes.csv", marg, delimiter=",", header="var,iqr", comments="", fmt="%g")
import pandas as pd
pd.DataFrame(envs, columns=["gremio", "var", "lat", "p5", "p95"]).to_csv(f"{INS}/envolventes.csv", index=False,
                                                                       float_format="%.4g")
print("->", INS)
