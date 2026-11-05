"""09b_divisoria_altiplano.py - M1: the Altiplano-Puna stretch of the divide (17-27 S) flooded in a single window.

On the endorheic plateau the W|E contact depends on passes far to the north or south, so the 1-degree blocks of
09_divisoria_glo30.py put it on the western rim in some blocks and on the eastern rim in others (seams at 18, 21, 22,
24, 25 and 26 S; a 6" test showed that no finite overlap fixes it, only a window over the whole plateau; BITACORA
2026-11-03). Same flooding and outputs as 09, over 16.85-27.15 S at once; to fit in memory the level overwrites the
elevation as the flood advances (a reached cell's elevation is never read again) and the crest is re-read per strip.
Overwrites blk_-17.csv .. blk_-26.csv in data/interim/M1_glo30_divide (the 09 versions are kept in _blocks_09/).
python -P scripts/09b_divisoria_altiplano.py   (from the repository root)
"""
import heapq
import importlib.util
import os
import shutil

import numpy as np
import rasterio
from numba import njit
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject
from skimage.measure import find_contours

spec = importlib.util.spec_from_file_location("d09", "scripts/09_divisoria_glo30.py")
d09 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d09)
RES, OV = d09.RES, d09.OV
LAT1 = range(-17, -27, -1)                        # core strips lat1-1 .. lat1
TOP, BOT = -16.85, -27.15


@njit
def flood_inplace(elev, lab):
    """09 flood with the level written into elev (seeds keep their elevation as level)."""
    nr, nc = elev.shape
    h = [(np.float32(0.0), np.int64(0))]
    h.pop()
    for r in range(nr):
        for c in range(nc):
            if lab[r, c] == 1 or lab[r, c] == 2:
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
                elev[rr, cc] = v
                heapq.heappush(h, (v, np.int64(rr * nc + cc)))


def read_glo(top, bot, west, nr, nc):
    """GLO-30 mosaic for rows top..bot from west, placed by integer offsets as in 09."""
    elev = np.full((nr, nc), np.nan, np.float32)
    for la in range(int(np.floor(bot)), int(np.ceil(top))):
        for lo in range(int(np.floor(west)), int(np.ceil(west + nc * RES))):
            f = os.path.join(d09.GLO, d09.tile_name(la, lo) + ".tif")
            if not os.path.exists(f):
                continue
            with rasterio.open(f) as s:
                a = s.read(1)
            r0 = int(round((top - (la + 1)) * 3600)); c0 = int(round((lo - west) * 3600))
            rs, cs = max(r0, 0), max(c0, 0)
            re, ce = min(r0 + 3600, nr), min(c0 + 3600, nc)
            if rs < re and cs < ce:
                elev[rs:re, cs:ce] = a[rs - r0:re - r0, cs - c0:ce - c0]
    return elev


def main():
    with rasterio.open(os.path.join(d09.INT, "M1_markers_15s.tif")) as s:
        mk = s.read(1); mk_tr = s.transform
    r0 = max(0, int((mk_tr.f - TOP) / -mk_tr.e)); r1 = min(mk.shape[0], int((mk_tr.f - BOT) / -mk_tr.e) + 1)
    cols = np.where((mk[r0:r1] > 0).any(axis=0))[0]
    west = np.floor((mk_tr.c + cols.min() * mk_tr.a - 0.05) * 3600) / 3600
    east = mk_tr.c + (cols.max() + 1) * mk_tr.a + 0.05
    nr, nc = int(round((TOP - BOT) * 3600)), int(round((east - west) * 3600))
    tr = from_origin(west, TOP, RES, RES)
    print("window", nr, "x", nc, "=", round(nr * nc / 1e9, 2), "G cells", flush=True)

    lvl = read_glo(TOP, BOT, west, nr, nc)                 # elevation now, level after the flood
    lab = np.zeros((nr, nc), np.uint8)
    reproject(mk, lab, src_transform=mk_tr, src_crs="EPSG:4326", dst_transform=tr, dst_crs="EPSG:4326",
              resampling=Resampling.nearest)
    del mk
    lab[np.isnan(lvl)] = 0
    np.nan_to_num(lvl, copy=False, nan=9999)
    flood_inplace(lvl, lab)
    lab[lab == 3] = 0
    print("flooded", flush=True)

    bak = os.path.join(d09.OUT, "_blocks_09"); os.makedirs(bak, exist_ok=True)
    for lat1 in LAT1:
        # strip with the same 0.15 deg margin as a 09 block; contours of the global labels, core rows written
        t, b = lat1 + OV, lat1 - 1 - OV
        a0, a1 = int(round((TOP - t) * 3600)), int(round((TOP - b) * 3600))
        lb = lab[a0:a1]; lv = lvl[a0:a1]
        crest_e = read_glo(t, b, west, a1 - a0, nc)
        out = []
        for k, cn in enumerate(find_contours((lb == 1).astype(np.float32), 0.5, mask=lb > 0)):
            r, c = cn[:, 0], cn[:, 1]
            horiz = np.abs(r - np.round(r)) < 1e-6
            ra = np.where(horiz, np.round(r), np.floor(r)).astype(int); ca = np.where(horiz, np.floor(c), np.round(c)).astype(int)
            rb = np.where(horiz, np.round(r), np.ceil(r)).astype(int); cb = np.where(horiz, np.ceil(c), np.round(c)).astype(int)
            cross = np.maximum(lv[ra, ca], lv[rb, cb]); crest = np.maximum(crest_e[ra, ca], crest_e[rb, cb])
            lat = t - (r + 0.5) * RES; lon = west + (c + 0.5) * RES
            core = (lat <= lat1) & (lat > lat1 - 1)
            if not core.any():
                continue
            run = np.cumsum(np.r_[True, np.diff(core.astype(int)) != 0])
            for g in np.unique(run[core]):
                j = np.where(run == g)[0]
                out.append(np.column_stack([np.full(len(j), k), j, lon[j], lat[j], cross[j], crest[j], np.full(len(j), g)]))
        f = os.path.join(d09.OUT, "blk_%+03d.csv" % lat1)
        if os.path.exists(f) and not os.path.exists(os.path.join(bak, os.path.basename(f))):
            shutil.copy(f, bak)
        res = np.vstack(out) if out else np.empty((0, 7))
        np.savetxt(f, res, delimiter=",", header="contour,seq,lon,lat,cross,crest,run", comments="",
                   fmt=["%d", "%d", "%.6f", "%.6f", "%.2f", "%.2f", "%d"])
        print(lat1, len(res), flush=True)


if __name__ == "__main__":
    main()
