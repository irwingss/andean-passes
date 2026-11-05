"""09_divisoria_glo30.py - M1: W-E divide and crossing heights on Copernicus GLO-30 (1 arc-second).

Marker-controlled flooding (priority flood, 8 neighbours) from the W and E seeds of 08_marcadores.R.
Each cell gets the side whose water reaches it first and its level = lowest possible maximum elevation
on a path from that side's seeds (minimax). The divide is the W|E label boundary (traced with
find_contours); the crossing height of a boundary edge between W cell a and E cell b is
max(level[a], level[b]) = the lowest W-E route forced through that edge. Its minima along the divide
are the cols. Endorheic basins are free cells, so a W|I stretch carries the I|E barrier automatically.

Processed in 1-degree latitude blocks with 0.15 deg overlap; only the core rows are written.
python -P scripts/09_divisoria_glo30.py   (run from the repository root; resumable, one CSV per block)
"""
import heapq
import os
import sys

import numpy as np
import rasterio
import yaml
from numba import njit
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject
from skimage.measure import find_contours

RES = 1 / 3600
OV = 0.15
CFG = yaml.safe_load(open("config.yml", encoding="utf-8"))["paths"]
INT = CFG["interim"]
GLO = os.path.join(CFG["data_dir"], "glo30")
OUT = os.path.join(INT, "M1_glo30_divide")


@njit
def flood(elev, lab, lvl):
    nr, nc = elev.shape
    h = [(np.float32(0.0), np.int64(0))]
    h.pop()
    for r in range(nr):
        for c in range(nc):
            if lab[r, c] == 1 or lab[r, c] == 2:
                lvl[r, c] = elev[r, c]
                heapq.heappush(h, (elev[r, c], np.int64(r * nc + c)))
    while len(h):
        lv, i = heapq.heappop(h)
        r = i // nc
        c = i % nc
        for dr in range(-1, 2):
            for dc in range(-1, 2):
                rr = r + dr
                cc = c + dc
                if rr < 0 or rr >= nr or cc < 0 or cc >= nc or lab[rr, cc] != 3:
                    continue
                lab[rr, cc] = lab[r, c]
                v = max(elev[rr, cc], lv)
                lvl[rr, cc] = v
                heapq.heappush(h, (v, np.int64(rr * nc + cc)))


def tile_name(la, lo):
    return "Copernicus_DSM_COG_10_%s%02d_00_%s%03d_00_DEM" % ("S" if la < 0 else "N", abs(la), "W" if lo < 0 else "E", abs(lo))


def block(lat1, mk, mk_tr):
    """Core rows lat1-1 .. lat1 (north edge lat1)."""
    top, bot = lat1 + OV, lat1 - 1 - OV
    # corridor longitude extent in this band, from the 15 s marker raster
    r0 = max(0, int((mk_tr.f - top) / -mk_tr.e)); r1 = min(mk.shape[0], int((mk_tr.f - bot) / -mk_tr.e) + 1)
    cols = np.where((mk[r0:r1] > 0).any(axis=0))[0]
    if not len(cols):
        return None
    west = np.floor((mk_tr.c + cols.min() * mk_tr.a - 0.05) * 3600) / 3600
    east = mk_tr.c + (cols.max() + 1) * mk_tr.a + 0.05
    nr, nc = int(round((top - bot) * 3600)), int(round((east - west) * 3600))
    tr = from_origin(west, top, RES, RES)

    elev = np.full((nr, nc), np.nan, np.float32)
    for la in range(int(np.floor(bot)), int(np.ceil(top))):
        for lo in range(int(np.floor(west)), int(np.ceil(east))):
            f = os.path.join(GLO, tile_name(la, lo) + ".tif")
            if not os.path.exists(f):
                continue
            # Tiles start exactly at whole degrees with 3600 rows; placed by integer offsets (reprojecting
            # tile by tile left 1-px no-data seams that acted as walls). Above 50 deg the 2400 (1.5")
            # columns are expanded to 3600 by nearest neighbour.
            with rasterio.open(f) as s:
                a = s.read(1)
            if a.shape[1] != 3600:
                a = a[:, (np.arange(3600) * a.shape[1]) // 3600]
            r0 = int(round((top - (la + 1)) * 3600)); c0 = int(round((lo - west) * 3600))
            rs, cs = max(r0, 0), max(c0, 0)
            re, ce = min(r0 + 3600, nr), min(c0 + 3600, nc)
            if rs < re and cs < ce:
                elev[rs:re, cs:ce] = a[rs - r0:re - r0, cs - c0:ce - c0]
    lab = np.zeros((nr, nc), np.uint8)
    reproject(mk, lab, src_transform=mk_tr, src_crs="EPSG:4326", dst_transform=tr, dst_crs="EPSG:4326",
              resampling=Resampling.nearest)
    lab[np.isnan(elev)] = 0
    lvl = np.full((nr, nc), np.nan, np.float32)
    flood(np.nan_to_num(elev, nan=9999).astype(np.float32), lab, lvl)
    lab[lab == 3] = 0                                   # free cells no flood reached

    img = (lab == 1).astype(np.float32)
    out = []
    for k, cn in enumerate(find_contours(img, 0.5, mask=lab > 0)):
        r, c = cn[:, 0], cn[:, 1]
        horiz = np.abs(r - np.round(r)) < 1e-6           # vertex between two cells of the same row
        ra = np.where(horiz, np.round(r), np.floor(r)).astype(int); ca = np.where(horiz, np.floor(c), np.round(c)).astype(int)
        rb = np.where(horiz, np.round(r), np.ceil(r)).astype(int); cb = np.where(horiz, np.ceil(c), np.round(c)).astype(int)
        cross = np.maximum(lvl[ra, ca], lvl[rb, cb]); crest = np.maximum(elev[ra, ca], elev[rb, cb])
        lat = top - (r + 0.5) * RES; lon = west + (c + 0.5) * RES
        core = (lat <= lat1) & (lat > lat1 - 1)
        if not core.any():
            continue
        # split into runs of consecutive core vertices
        run = np.cumsum(np.r_[True, np.diff(core.astype(int)) != 0])
        for g in np.unique(run[core]):
            j = np.where(run == g)[0]
            out.append(np.column_stack([np.full(len(j), k), j, lon[j], lat[j], cross[j], crest[j], np.full(len(j), g)]))
    return np.vstack(out) if out else None


def main():
    os.makedirs(OUT, exist_ok=True)
    # every corridor tile must be present (ocean-only tiles that do not exist on AWS are listed as absent)
    absent = set(open(os.path.join(INT, "M1_glo30_absent.txt")).read().split()) if os.path.exists(os.path.join(INT, "M1_glo30_absent.txt")) else set()
    need = [t for t in open(os.path.join(INT, "M1_glo30_tiles.txt")).read().split() if t not in absent]
    miss = [t for t in need if not os.path.exists(os.path.join(GLO, t + ".tif"))]
    if miss:
        sys.exit("missing GLO-30 tiles: %d (e.g. %s)" % (len(miss), miss[0]))
    with rasterio.open(os.path.join(INT, "M1_markers_15s.tif")) as s:
        mk = s.read(1)
        mk_tr = s.transform
    lats = range(13, -56, -1) if len(sys.argv) < 2 else [int(x) for x in sys.argv[1:]]
    for lat1 in lats:
        f = os.path.join(OUT, "blk_%+03d.csv" % lat1)
        if os.path.exists(f):
            continue
        res = block(lat1, mk, mk_tr)
        hdr = "contour,seq,lon,lat,cross,crest,run"
        np.savetxt(f, res if res is not None else np.empty((0, 7)), delimiter=",", header=hdr, comments="",
                   fmt=["%d", "%d", "%.6f", "%.6f", "%.2f", "%.2f", "%d"])
        print(lat1, 0 if res is None else len(res), flush=True)


if __name__ == "__main__":
    main()
