"""39_descarga_fabdem.py - FABDEM V1-2 tiles for every census block window (PRD O1), into data/raw/fabdem.
Source: the FABDEM V1-2 tiles in the Earth Engine community catalogue (projects/sat-io/open-datasets/FABDEM), the same
files as data.bris (https://doi.org/10.5523/bris.s5hqmjcdj8yo2ibzi9b4ew3sn). Pixels are read on their native grid,
without resampling; a tile read this way is identical, bit for bit, to the one in the data.bris zip. Earth Engine
stores the tiles half a cell off, so the request grid is shifted by half a cell and the file keeps the data.bris origin.
python -P scripts/39_descarga_fabdem.py <earth_engine_project>     needs Earth Engine access for that project.
"""
import importlib.util
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import ee
import numpy as np
import rasterio
import shapely
from rasterio.transform import from_origin

spec = importlib.util.spec_from_file_location("c40", "scripts/40_censo_bloque.py")
c40 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c40)
DIR = os.path.join(c40.RAW, "fabdem")
os.makedirs(DIR, exist_ok=True)
ee.Initialize(project=sys.argv[1])
COL = "projects/sat-io/open-datasets/FABDEM/"
R = 1 / 3600

# tiles touched by the window of every 1-degree block (core plus 0.5 degree each side)
area = c40.area_mask()
tiles = set()
for north in range(12, -57, -1):
    g = shapely.intersection(area, shapely.box(-90, north - 1, -55, north))
    if g.is_empty:
        continue
    w0, _, e0, _ = g.bounds
    for la in range(math.floor(north - 1.5), math.ceil(north + 0.5)):
        for lo in range(math.floor(w0 - 0.5), math.ceil(e0 + 0.5)):
            tiles.add((la, lo))
w, s, e, n = area.bounds
have = set(ee.ImageCollection(COL.rstrip("/")).filterBounds(ee.Geometry.Rectangle([w - 2, s - 2, e + 2, n + 2]))
           .aggregate_array("system:index").getInfo())


def get(t):
    la, lo = t
    name = f"{c40.tile_name(la, lo)}_FABDEM_V1-2"
    out = os.path.join(DIR, name + ".tif")
    if os.path.exists(out):
        return "kept"
    if name not in have:
        return "ocean"                       # no tile in FABDEM: open sea
    img = ee.Image(COL + name)
    parts = []
    for h in (0, 1):                         # two halves: one request is limited to 48 MB
        grid = {"dimensions": {"width": 3600, "height": 1800}, "crsCode": "EPSG:4326",
                "affineTransform": {"scaleX": R, "shearX": 0, "translateX": lo - R / 2, "shearY": 0, "scaleY": -R,
                                    "translateY": la + 1 + R / 2 - h * 0.5}}
        for k in range(6):
            try:
                parts.append(ee.data.computePixels({"expression": img, "fileFormat": "NUMPY_NDARRAY",
                                                    "grid": grid})["b1"])
                break
            except Exception:
                if k == 5:
                    raise
                time.sleep(10 * (k + 1))
    a = np.vstack(parts).astype(np.float32)
    prof = dict(driver="GTiff", height=3600, width=3600, count=1, dtype="float32", crs="EPSG:4326",
                transform=from_origin(lo, la + 1, R, R), nodata=-9999.0, compress="deflate", predictor=3, tiled=True)
    with rasterio.open(out + ".part", "w", **prof) as d:
        d.write(a, 1)
    os.replace(out + ".part", out)
    return "new"


with ThreadPoolExecutor(16) as ex:      # a failed tile raises here and stops the census: no silent gaps
    res = list(ex.map(get, sorted(tiles)))
print({k: res.count(k) for k in set(res)}, "tiles ->", DIR)
