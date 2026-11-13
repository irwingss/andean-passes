"""40_censo_bloque.py - Census of passes on one latitude block, from FABDEM relief alone (PRD O1-O3).
python -P scripts/40_censo_bloque.py <north_edge_of_core_deg> [core_deg=1] [overlap_deg=0.5]
e.g. -5 -> core 5-6 S, window 4.5-6.5 S. Output: data/interim/censo/passes_<tag>.csv (core passes only).

Steps: FABDEM mosaic -> real sinks: depressions with the bottom in an endorheic HydroBASINS basin and 20 m or more
deep, kept as in FABDEM -> breaching of every other depression (Lindsay 2016, WhiteboxTools, no depth or length limit,
pits raised, no fill)
-> D8 with flats drained to their outlets (Barnes, Lehman & Mulla 2014b) -> accumulation -> channel links at 0.1 km2
-> every cell labelled by the link it drains to -> divides = cell edges between labels (never a channel cell)
-> passes and depth by a top-down merge of the divide graph (prominence principle) -> width along the crest.
Heights of passes are always read from FABDEM without conditioning.
"""
import glob
import math
import os
import sys
import zipfile
from heapq import heappop, heappush

import numba as nb
import numpy as np
import pyogrio
import rasterio
import shapely
import yaml
from pyproj import Transformer
from rasterio.transform import from_origin

cfg = yaml.safe_load(open("config.yml", encoding="utf-8"))
RAW, INT = cfg["paths"]["data_dir"], cfg["paths"]["interim"]
OUT = os.path.join(INT, "censo")
TMP = os.path.join(OUT, "tmp")
os.makedirs(TMP, exist_ok=True)

RES = 1 / 3600
CH_AREA = 0.1e6          # m2, minimum scale of a divide (O2.4)
MIN_DEPTH = 30.0         # D15
WIDTH_DEPTH = 100.0      # width only for passes this deep (O3)
NODATA = -9999.0
R_EARTH = 6371008.8
DI = np.array([-1, -1, 0, 1, 1, 1, 0, -1], np.int64)   # D8 neighbours, clockwise from north
DJ = np.array([0, 1, 1, 1, 0, -1, -1, -1], np.int64)


# ---------------------------------------------------------------- area and window
def area_mask():
    geom = shapely.from_wkb(pyogrio.raw.read(os.path.join(INT, "M1_gmba_andes.gpkg"))[2][0])
    aea = "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +datum=WGS84 +units=m"
    fw, bw = Transformer.from_crs(4326, aea, always_xy=True), Transformer.from_crs(aea, 4326, always_xy=True)
    m = shapely.transform(geom, lambda xy: np.column_stack(fw.transform(xy[:, 0], xy[:, 1])))
    m = shapely.buffer(m, cfg["params"]["gmba_margin_km"] * 1000)
    return shapely.transform(m, lambda xy: np.column_stack(bw.transform(xy[:, 0], xy[:, 1])))


def tile_name(lat, lon):          # lat, lon = lower-left integer corner
    return f"{'N' if lat >= 0 else 'S'}{abs(lat):02d}{'E' if lon >= 0 else 'W'}{abs(lon):03d}"


def zip_name(lat, lon):           # data.bris groups FABDEM in 10 x 10 degree zips named by their corners
    la, lo = math.floor(lat / 10) * 10, math.floor(lon / 10) * 10
    return f"{tile_name(la, lo)}-{tile_name(la + 10, lo + 10)}_FABDEM_V1-2.zip"


def mosaic(n, s, w, e):
    """FABDEM tiles placed by integer offset (no resampling). Missing tiles (ocean) stay NODATA."""
    nr, nc = round((n - s) * 3600), round((e - w) * 3600)
    z = np.full((nr, nc), NODATA, np.float32)
    for la in range(math.floor(s), math.ceil(n)):
        for lo in range(math.floor(w), math.ceil(e)):
            member = f"{tile_name(la, lo)}_FABDEM_V1-2.tif"
            if os.path.isdir(os.path.join(RAW, "fabdem")):     # tiles from 39_descarga_fabdem.py
                tif = os.path.join(RAW, "fabdem", member)
                if not os.path.exists(tif):
                    continue                                   # no FABDEM tile: open sea
                fh = open(tif, "rb")
            else:                                              # the data.bris zips
                zz = zipfile.ZipFile(os.path.join(RAW, zip_name(la, lo)))
                if member not in zz.namelist():
                    continue
                fh = zz.open(member)
            with fh, rasterio.MemoryFile(fh.read()) as mf, mf.open() as src:
                t = src.read(1).astype(np.float32)
                if src.nodata is not None:
                    t[t == src.nodata] = NODATA
                assert t.shape == (3600, 3600) and abs(src.transform.a - RES) < 1e-12
                top, left = src.transform.f, src.transform.c
            r0, c0 = round((n - top) * 3600), round((left - w) * 3600)
            ra, ca = max(r0, 0), max(c0, 0)
            rb, cb = min(r0 + 3600, nr), min(c0 + 3600, nc)
            if ra < rb and ca < cb:
                z[ra:rb, ca:cb] = t[ra - r0:rb - r0, ca - c0:cb - c0]
    # FABDEM gives the sea as 0: zero regions touching the window edge become NODATA (an outlet). Land below sea
    # level (Sechura) is negative, so it stays.
    from scipy import ndimage
    lab, _ = ndimage.label(z == 0)
    edge = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
    z[np.isin(lab, edge[edge > 0])] = NODATA
    return z


# ---------------------------------------------------------------- D8 with flats (Barnes et al. 2014b)
@nb.njit(cache=True)
def d8_flats(z, dx, dy, nod):
    nr, nc = z.shape
    rcv = np.full((nr, nc), -1, np.int8)     # 0..7 direction, -1 no flow, -2 outlet to the sea or a NODATA cell
    dd = np.empty(8)
    for i in range(nr):
        for j in range(nc):
            if z[i, j] == nod:
                continue
            dd[0] = dd[4] = dy
            dd[2] = dd[6] = dx[i]
            dd[1] = dd[3] = dd[5] = dd[7] = math.sqrt(dx[i] ** 2 + dy ** 2)
            best, bk, edge = 0.0, -1, False
            for k in range(8):
                ii, jj = i + DI[k], j + DJ[k]
                if ii < 0 or jj < 0 or ii >= nr or jj >= nc:     # the window edge is not an outlet (O2.5)
                    continue
                if z[ii, jj] == nod:
                    edge = True
                    continue
                s = (z[i, j] - z[ii, jj]) / dd[k]
                if s > best:
                    best, bk = s, k
            if bk >= 0:
                rcv[i, j] = bk
            elif edge:
                rcv[i, j] = -2
    # flats: label from low edges, gradients away from higher and towards lower terrain
    lab = np.zeros((nr, nc), np.int32)
    mask = np.zeros((nr, nc), np.int32)
    q = np.empty(nr * nc // 4 + 1024, np.int64)   # ring buffer
    qn = q.size
    nlab = 0
    for i in range(nr):
        for j in range(nc):
            if rcv[i, j] == -1 or z[i, j] == nod:
                continue
            low = False
            for k in range(8):
                ii, jj = i + DI[k], j + DJ[k]
                if 0 <= ii < nr and 0 <= jj < nc and rcv[ii, jj] == -1 and z[ii, jj] == z[i, j]:
                    low = True
                    break
            if not low or lab[i, j] != 0:
                continue
            nlab += 1
            lab[i, j] = nlab
            h = t = 0
            q[t] = i * nc + j
            t += 1
            while h != t:
                c = q[h]; h = (h + 1) % qn
                ci, cj = c // nc, c % nc
                for k in range(8):
                    ii, jj = ci + DI[k], cj + DJ[k]
                    if 0 <= ii < nr and 0 <= jj < nc and lab[ii, jj] == 0 and z[ii, jj] == z[ci, cj] \
                            and (rcv[ii, jj] == -1 or rcv[ii, jj] >= 0):
                        lab[ii, jj] = nlab
                        q[t] = ii * nc + jj; t = (t + 1) % qn
    fh = np.zeros(nlab + 1, np.int32)
    # away from higher terrain
    h = t = 0
    for i in range(nr):
        for j in range(nc):
            if rcv[i, j] != -1 or lab[i, j] == 0:
                continue
            for k in range(8):
                ii, jj = i + DI[k], j + DJ[k]
                if 0 <= ii < nr and 0 <= jj < nc and z[ii, jj] != nod and z[ii, jj] > z[i, j]:
                    q[t] = i * nc + j; t = (t + 1) % qn
                    break
    loops = 1
    q[t] = -1; t = (t + 1) % qn
    while h != t:
        c = q[h]; h = (h + 1) % qn
        if c == -1:
            if h == t:
                break
            loops += 1
            q[t] = -1; t = (t + 1) % qn
            continue
        ci, cj = c // nc, c % nc
        if mask[ci, cj] > 0:
            continue
        mask[ci, cj] = loops
        fh[lab[ci, cj]] = loops
        for k in range(8):
            ii, jj = ci + DI[k], cj + DJ[k]
            if 0 <= ii < nr and 0 <= jj < nc and lab[ii, jj] == lab[ci, cj] and rcv[ii, jj] == -1 \
                    and mask[ii, jj] == 0:
                q[t] = ii * nc + jj; t = (t + 1) % qn
    # towards lower terrain
    for i in range(nr):
        for j in range(nc):
            mask[i, j] = -mask[i, j]
    h = t = 0
    for i in range(nr):
        for j in range(nc):
            if lab[i, j] == 0 or rcv[i, j] == -1:
                continue
            for k in range(8):
                ii, jj = i + DI[k], j + DJ[k]
                if 0 <= ii < nr and 0 <= jj < nc and rcv[ii, jj] == -1 and lab[ii, jj] == lab[i, j]:
                    q[t] = i * nc + j; t = (t + 1) % qn
                    break
    loops = 1
    q[t] = -1; t = (t + 1) % qn
    while h != t:
        c = q[h]; h = (h + 1) % qn
        if c == -1:
            if h == t:
                break
            loops += 1
            q[t] = -1; t = (t + 1) % qn
            continue
        ci, cj = c // nc, c % nc
        if mask[ci, cj] > 0:
            continue
        if mask[ci, cj] < 0:
            mask[ci, cj] = fh[lab[ci, cj]] + mask[ci, cj] + 2 * loops
        else:
            mask[ci, cj] = 2 * loops
        for k in range(8):
            ii, jj = ci + DI[k], cj + DJ[k]
            if 0 <= ii < nr and 0 <= jj < nc and lab[ii, jj] == lab[ci, cj] and rcv[ii, jj] == -1 \
                    and mask[ii, jj] <= 0:
                q[t] = ii * nc + jj; t = (t + 1) % qn
    # directions inside flats: towards the lowest mask of the same flat
    for i in range(nr):
        for j in range(nc):
            if rcv[i, j] != -1 or lab[i, j] == 0:
                continue
            bm, bk = mask[i, j], -1
            for k in range(8):
                ii, jj = i + DI[k], j + DJ[k]
                if 0 <= ii < nr and 0 <= jj < nc and lab[ii, jj] == lab[i, j] and mask[ii, jj] < bm:
                    bm, bk = mask[ii, jj], k
            if bk >= 0:
                rcv[i, j] = bk
    # flats without outlet never got a mask: they stay -1 (bottom of a real sink)
    return rcv, lab


# ---------------------------------------------------------------- accumulation, links, labels
@nb.njit(cache=True)
def labels_from_flow(rcv, flat, cell_area):
    nr, nc = rcv.shape
    n = nr * nc
    indeg = np.zeros(n, np.int8)
    for c in range(n):
        k = rcv.flat[c]
        if k >= 0:
            indeg[(c // nc + DI[k]) * nc + c % nc + DJ[k]] += 1
    order = np.empty(n, np.int32)
    t = 0
    for c in range(n):
        if indeg[c] == 0:
            order[t] = c; t += 1
    h = 0
    acc = np.empty(n, np.float32)
    for c in range(n):
        acc[c] = cell_area[c // nc]
    while h < t:
        c = order[h]; h += 1
        k = rcv.flat[c]
        if k >= 0:
            r = (c // nc + DI[k]) * nc + c % nc + DJ[k]
            acc[r] += acc[c]
            indeg[r] -= 1
            if indeg[r] == 0:
                order[t] = r; t += 1
    ch = acc >= CH_AREA
    ndch = indeg            # reuse, all zero now
    for c in range(n):
        k = rcv.flat[c]
        if ch[c] and k >= 0:
            r = (c // nc + DI[k]) * nc + c % nc + DJ[k]
            if ch[r]:
                ndch[r] += 1
    lab = np.zeros(n, np.int32)
    term = np.empty(n, np.int64)          # where each cell ends: its pit or outlet cell, or n + flat of a sink bottom
    nid = 0
    for h in range(t):                    # upstream first: links between confluences
        c = order[h]
        if not ch[c]:
            continue
        if ndch[c] != 1:
            nid += 1
            lab[c] = nid
        k = rcv.flat[c]
        if k >= 0:
            r = (c // nc + DI[k]) * nc + c % nc + DJ[k]
            if ch[r] and ndch[r] == 1:
                lab[r] = lab[c]
    for h in range(t - 1, -1, -1):        # downstream first: every cell takes the link it drains to
        c = order[h]
        k = rcv.flat[c]
        if k >= 0:
            term[c] = term[(c // nc + DI[k]) * nc + c % nc + DJ[k]]
        elif k == -1 and flat.flat[c] > 0:
            term[c] = n + flat.flat[c]
        else:
            term[c] = c
        if ch[c]:
            continue
        if k >= 0:
            lab[c] = lab[(c // nc + DI[k]) * nc + c % nc + DJ[k]]
        elif k == -1 and flat.flat[c] > 0:     # bottom of a real sink: one label per flat
            lab[c] = n + flat.flat[c]
        else:
            nid += 1
            lab[c] = nid
    return lab.reshape(nr, nc), ch.reshape(nr, nc), term.reshape(nr, nc), acc.reshape(nr, nc), t


# ---------------------------------------------------------------- divide graph and passes
@nb.njit(cache=True)
def is_div(lab, ch, z, nod, i, j, i2, j2):
    nr, nc = lab.shape
    if i < 0 or j < 0 or i2 >= nr or j2 >= nc or i2 < 0 or j2 < 0 or i >= nr or j >= nc:
        return False
    return lab[i, j] != lab[i2, j2] and not ch[i, j] and not ch[i2, j2] and z[i, j] != nod and z[i2, j2] != nod


@nb.njit(cache=True)
def list_cracks(lab, ch, z, nod):
    nr, nc = lab.shape
    m = 0
    for i in range(nr):
        for j in range(nc):
            if is_div(lab, ch, z, nod, i, j, i, j + 1):
                m += 1
            if is_div(lab, ch, z, nod, i, j, i + 1, j):
                m += 1
    key = np.empty(m, np.int64)
    w = np.empty(m, np.float32)
    m = 0
    for i in range(nr):
        for j in range(nc):
            if is_div(lab, ch, z, nod, i, j, i, j + 1):
                key[m] = 2 * (i * nc + j); w[m] = max(z[i, j], z[i, j + 1]); m += 1
            if is_div(lab, ch, z, nod, i, j, i + 1, j):
                key[m] = 2 * (i * nc + j) + 1; w[m] = max(z[i, j], z[i + 1, j]); m += 1
    return key, w


@nb.njit(cache=True)
def group(rcv, ci, cj, e):
    """Which part of corner (ci, cj) the crack on side e (0 up, 1 right, 2 down, 3 left) belongs to. A corner crossed
    by a diagonal D8 flow is split in two, so a divide never crosses a flow path; crossed both ways, in four."""
    nr, nc = rcv.shape
    d1 = (ci >= 1 and cj >= 1 and rcv[ci - 1, cj - 1] == 3) or (ci < nr and cj < nc and rcv[ci, cj] == 7)
    d2 = (ci >= 1 and cj < nc and rcv[ci - 1, cj] == 5) or (ci < nr and cj >= 1 and rcv[ci, cj - 1] == 1)
    if d1 and d2:
        return e
    if d1:
        return 0 if e <= 1 else 1
    if d2:
        return 0 if (e == 0 or e == 3) else 1
    return 0


@nb.njit(cache=True)
def endpoints(key, rcv):
    """Node ids (4 * corner + part) at both ends of each crack. Right crack of cell (i,j) joins corners (i,j+1)
    [its side 2] and (i+1,j+1) [side 0]; down crack joins (i+1,j) [side 1] and (i+1,j+1) [side 3]."""
    nc = rcv.shape[1]
    ea, eb = np.empty(key.size, np.int64), np.empty(key.size, np.int64)
    for n in range(key.size):
        c = key[n] // 2
        i, j = c // nc, c % nc
        if key[n] % 2 == 0:
            ea[n] = 4 * (i * (nc + 1) + j + 1) + group(rcv, i, j + 1, 2)
            eb[n] = 4 * ((i + 1) * (nc + 1) + j + 1) + group(rcv, i + 1, j + 1, 0)
        else:
            ea[n] = 4 * ((i + 1) * (nc + 1) + j) + group(rcv, i + 1, j, 1)
            eb[n] = 4 * ((i + 1) * (nc + 1) + j + 1) + group(rcv, i + 1, j + 1, 3)
    return ea, eb


@nb.njit(cache=True)
def find(par, a):
    r = a
    while par[r] != r:
        r = par[r]
    while par[a] != r:
        nx = par[a]; par[a] = r; a = nx
    return r


@nb.njit(cache=True)
def merge_passes(na, nb_, w, order, nedge, min_depth):
    """Top-down merge of the divide graph: each crack where two higher parts of the network first meet is a pass;
    depth = lower of the two maxima - crack height. Returns crack indices. The depth is only a lower bound (flag) when
    the lower part reaches the window edge, where its true maximum may lie beyond."""
    nnodes = nedge.size
    par = np.full(nnodes, -1, np.int64)
    cmax = np.zeros(nnodes, np.float32)
    cedge = nedge.copy()
    m = w.size
    oi, od = np.empty(m, np.int64), np.empty(m, np.float32)
    ohi, oring = np.empty(m, np.float32), np.empty(m, np.bool_)
    no = 0
    for idx in order:
        a, b = na[idx], nb_[idx]
        wi = w[idx]
        if par[a] == -1:
            par[a] = a; cmax[a] = wi
        if par[b] == -1:
            par[b] = b; cmax[b] = wi
        ra, rb = find(par, a), find(par, b)
        if ra == rb:      # closes a loop: the rim of a closed basin, handled by sink_spills
            continue
        lo, hi = min(cmax[ra], cmax[rb]), max(cmax[ra], cmax[rb])
        if cmax[ra] < cmax[rb]:
            trunc = cedge[ra]
        elif cmax[rb] < cmax[ra]:
            trunc = cedge[rb]
        else:
            trunc = cedge[ra] and cedge[rb]
        if lo - wi >= min_depth:
            oi[no], od[no], ohi[no], oring[no] = idx, lo - wi, hi, trunc
            no += 1
        if cmax[ra] >= cmax[rb]:
            par[rb] = ra; cedge[ra] = cedge[ra] or cedge[rb]
        else:
            par[ra] = rb; cedge[rb] = cedge[ra] or cedge[rb]
    return oi[:no].copy(), od[:no].copy(), ohi[:no].copy(), oring[:no].copy()


@nb.njit(cache=True)
def crack_w(lab, ch, z, nod, k, nc):
    c = k // 2
    i, j = c // nc, c % nc
    i2, j2 = (i, j + 1) if k % 2 == 0 else (i + 1, j)
    if not is_div(lab, ch, z, nod, i, j, i2, j2):
        return -1e30
    return max(z[i, j], z[i2, j2])


@nb.njit(cache=True)
def side_width(lab, ch, z, nod, rcv, start, banned, w0, top, dx, dy, cap):
    """Shortest length along the divide from node `start` (not crossing crack `banned`) to the first crack higher
    than `top`, using only cracks between w0 and top. Returns -1 if not reached within `cap` nodes."""
    nr, nc = lab.shape
    dist = {start: 0.0}
    hp = [(0.0, start)]
    seen = 0
    while len(hp) > 0:
        d, v = heappop(hp)
        if d > dist[v]:
            continue
        seen += 1
        if seen > cap:
            return -1.0
        corner, s = v // 4, v % 4
        ci, cj = corner // (nc + 1), corner % (nc + 1)
        if ci <= 0 or cj <= 0 or ci >= nr or cj >= nc:     # reached the window edge: unknown
            return -1.0
        for e in range(4):
            if group(rcv, ci, cj, e) != s:
                continue
            if e == 0:    # up: right crack of cell (ci-1, cj-1); other end (ci-1, cj) side 2
                okk = ci >= 1 and cj >= 1 and cj < nc + 1
                k, oi, oj, oe, ln = 2 * ((ci - 1) * nc + cj - 1), ci - 1, cj, 2, dy
            elif e == 2:  # down: right crack of cell (ci, cj-1); other end (ci+1, cj) side 0
                okk = ci < nr and cj >= 1
                k, oi, oj, oe, ln = 2 * (ci * nc + cj - 1), ci + 1, cj, 0, dy
            elif e == 3:  # left: down crack of cell (ci-1, cj-1); other end (ci, cj-1) side 1
                okk = ci >= 1 and cj >= 1
                k, oi, oj, oe, ln = 2 * ((ci - 1) * nc + cj - 1) + 1, ci, cj - 1, 1, dx[max(ci - 1, 0)]
            else:         # right: down crack of cell (ci-1, cj); other end (ci, cj+1) side 3
                okk = ci >= 1 and cj < nc
                k, oi, oj, oe, ln = 2 * ((ci - 1) * nc + cj) + 1, ci, cj + 1, 3, dx[max(ci - 1, 0)]
            if not okk or k == banned:
                continue
            ww = crack_w(lab, ch, z, nod, k, nc)
            if ww < w0:
                continue
            if ww > top:
                return d + ln / 2
            u = 4 * (oi * (nc + 1) + oj) + group(rcv, oi, oj, oe)
            nd = d + ln
            if u not in dist or nd < dist[u]:
                dist[u] = nd
                heappush(hp, (nd, u))
    return -1.0


@nb.njit(cache=True)
def widths(lab, ch, z, nod, rcv, keys, ea, eb, ws, dx, dy, cap):
    nc = lab.shape[1]
    out = np.full(keys.size, np.nan)
    for n in range(keys.size):
        i = (keys[n] // 2) // nc
        own = dy if keys[n] % 2 == 0 else dx[i]
        top = ws[n] + WIDTH_DEPTH
        wa = side_width(lab, ch, z, nod, rcv, ea[n], keys[n], ws[n], top, dx, dy, cap)
        wb = side_width(lab, ch, z, nod, rcv, eb[n], keys[n], ws[n], top, dx, dy, cap)
        if wa >= 0 and wb >= 0:
            out[n] = wa + wb + own
    return out


def sink_spills(key, w, term, rcv, flat):
    """Pass of each closed basin: the lowest crack of its rim. Along a closed rim both directions climb to the same
    highest point, so depth = highest rim crack - lowest rim crack. Returns, per sink with a rim: its terminal id,
    lowest rim crack, depth, rim maximum, whether the basin reaches the window edge (rim incomplete) and whether its
    bottom is on the edge (an artefact of the window)."""
    nr, nc = term.shape
    n = nr * nc
    bord = np.concatenate([term[0], term[-1], term[:, 0], term[:, -1]])
    bflat = np.concatenate([flat[0], flat[-1], flat[:, 0], flat[:, -1]])
    bottom_on_edge = np.concatenate([np.unique(bord[bord < n]), n + np.unique(bflat[bflat > 0])])
    touches = np.unique(bord)
    c = key // 2
    a = term.ravel()[c]
    b = term.ravel()[np.where(key % 2 == 0, c + 1, c + nc)]
    t = np.concatenate([a, b])
    other = np.concatenate([b, a])
    idx = np.concatenate([np.arange(key.size)] * 2)
    sink = (t >= n) | (rcv.ravel()[np.minimum(t, n - 1)] == -1)
    sel = sink & (t != other)
    t, idx = t[sel], idx[sel]
    if t.size == 0:
        e = np.zeros(0, np.int64)
        return e, e, np.zeros(0, np.float32), np.zeros(0, np.float32), np.zeros(0, bool), np.zeros(0, bool)
    o = np.lexsort((w[idx], t))
    t, idx = t[o], idx[o]
    first = np.r_[0, np.flatnonzero(np.diff(t)) + 1]
    last = np.r_[first[1:] - 1, t.size - 1]
    lo, hi = idx[first], idx[last]
    tt = t[first]
    return tt, lo, w[hi] - w[lo], w[hi], np.isin(tt, touches), np.isin(tt, bottom_on_edge)


def write_sinks(tag, term, rcv, flat, raw, zfill, n, w, dx, dy, area, core, north, st, si, sd, shi, str_, sedge, key, wcr):
    """Every real sink kept by the breaching rule, for review: bottom, contributing area and rim pass."""
    nr, nc = rcv.shape
    N = nr * nc
    t = term.ravel()
    bot = np.flatnonzero((rcv.ravel() == -1) & (raw.ravel() != NODATA))
    bt = t[bot]
    o = np.lexsort((raw.ravel()[bot], bt))             # lowest bottom cell of each sink
    bt, bot = bt[o], bot[o]
    f = np.r_[True, bt[1:] != bt[:-1]]
    sterm, scell = bt[f], bot[f]
    insink = np.isin(t, sterm)
    ids, inv = np.unique(t[insink], return_inverse=True)
    rows = np.flatnonzero(insink) // nc
    a_km2 = np.bincount(inv, weights=(dx[rows] * dy)) / 1e6
    a_km2 = a_km2[np.searchsorted(ids, sterm)]
    pos = {v: k for k, v in enumerate(st)}
    j = np.array([pos.get(v, -1) for v in sterm])
    bi, bj = scell // nc, scell % nc
    la, lo = n - (bi + 0.5) * RES, w + (bj + 0.5) * RES
    import csv
    with open(os.path.join(OUT, f"sumideros_{tag}.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["lat", "lon", "bottom_z_m", "fill_depth_m", "contrib_km2", "flat_bottom", "spill_lat", "spill_lon", "spill_z_m",
                     "depth_m", "rim_max_m", "rim_incomplete", "bottom_on_edge", "core", "in_area"])
        inarea = shapely.contains_xy(area, lo, la)
        for k in range(sterm.size):
            sp = ["", "", "", "", "", "", ""]
            if j[k] >= 0:
                c = key[si[j[k]]] // 2
                ci, cj = c // nc, c % nc
                right = key[si[j[k]]] % 2 == 0
                sp = [f"{n - (ci + 0.5) * RES if right else n - (ci + 1) * RES:.6f}",
                      f"{w + (cj + 1) * RES if right else w + (cj + 0.5) * RES:.6f}",
                      f"{wcr[si[j[k]]]:.1f}", f"{sd[j[k]]:.1f}", f"{shi[j[k]]:.1f}", int(str_[j[k]]), int(sedge[j[k]])]
            wr.writerow([f"{la[k]:.6f}", f"{lo[k]:.6f}", f"{raw.ravel()[scell[k]]:.1f}",
                         f"{zfill.ravel()[scell[k]] - raw.ravel()[scell[k]]:.1f}", f"{a_km2[k]:.4f}",
                         int(sterm[k] >= N)] + sp + [int(north - core < la[k] <= north), int(inarea[k])])


@nb.njit(cache=True)
def trace_down(rcv, acc, ii, jj, thr):
    """Follow D8 from each cell until the path drains `thr` m2, or ends in a sink, the sea or the window edge."""
    nr, nc = rcv.shape
    oi, oj = ii.copy(), jj.copy()
    for n in range(ii.size):
        i, j = ii[n], jj[n]
        while acc[i, j] < thr and rcv[i, j] >= 0:
            k = rcv[i, j]
            i, j = i + DI[k], j + DJ[k]
        oi[n], oj[n] = i, j
    return oi, oj


# ---------------------------------------------------------------- main
def run(north, core=1.0, ov=0.5):
    tag = f"{'S' if north <= 0 else 'N'}{abs(north):05.2f}".replace(".", "p")
    area = area_mask()
    band = shapely.box(-90, north - core, -55, north)
    w0, s0, e0, n0 = shapely.intersection(area, band).bounds
    n, s = north + ov, north - core - ov
    w, e = math.floor((w0 - ov) * 3600) / 3600, math.ceil((e0 + ov) * 3600) / 3600
    print(f"{tag}: window {n} to {s}, {w:.3f} to {e:.3f}", flush=True)

    raw = mosaic(n, s, w, e)
    nr, nc = raw.shape
    tr = from_origin(w, n, RES, RES)
    prof = dict(driver="GTiff", height=nr, width=nc, count=1, dtype="float32", crs="EPSG:4326", transform=tr,
                nodata=NODATA, compress="deflate", tiled=True, BIGTIFF="YES")
    fr = os.path.abspath(os.path.join(TMP, f"raw_{tag}.tif"))
    fi, fb = (os.path.abspath(os.path.join(TMP, f"{p}_{tag}.tif")) for p in ("breach_in", "breach_endo"))
    ff = os.path.abspath(os.path.join(TMP, f"fill_{tag}.tif"))
    with rasterio.open(fr, "w", **prof) as d:
        d.write(raw, 1)
    import whitebox
    wbt = whitebox.WhiteboxTools()
    wbt.set_verbose_mode(True)              # its messages are kept and printed only if a tool fails
    msgs = []

    def wbt_ok(ret, out):
        if ret != 0 or not os.path.exists(out):
            print("WhiteboxTools failed:", *msgs[-20:], flush=True)
            return False
        return True
    if not os.path.exists(ff):          # fill level of raw FABDEM: depth of every depression
        assert wbt_ok(wbt.fill_depressions(fr, ff, fix_flats=False, callback=msgs.append), ff)
    # O2.2: a depression is real if its bottom (lowest raw cell) lies in a HydroBASINS level-12 basin marked
    # endorheic and its spill is 20 m or more above that bottom. Its bottom cells become NODATA, which WhiteboxTools
    # treats as an outlet, so the breach drains the depression to its bottom and never cuts its rim. Everything
    # else is breached, with no depth or length limit; single-cell pits are raised first.
    from scipy import ndimage
    with rasterio.open(ff) as src:
        zfill = src.read(1)
    dep_l, ndep = ndimage.label((zfill > raw + 1e-4) & (raw != NODATA), structure=np.ones((3, 3)))
    ids = np.arange(1, ndep + 1)
    bot_z = ndimage.minimum(raw, dep_l, ids)
    depth_d = ndimage.maximum(zfill, dep_l, ids) - bot_z
    del zfill
    bpos = np.array(ndimage.minimum_position(raw, dep_l, ids)).reshape(-1, 2)
    blat, blon = n - (bpos[:, 0] + 0.5) * RES, w + (bpos[:, 1] + 0.5) * RES
    f12 = os.path.join(RAW, "hybas_sa_lev01-12_v1c/hybas_sa_lev12_v1c.shp")
    _, _, g, F = pyogrio.raw.read(f12, columns=["ENDO"], bbox=(w, s, e, n))
    endo_poly = shapely.union_all(shapely.from_wkb(g)[F[0] > 0])
    real = (depth_d >= 20) & shapely.contains_xy(endo_poly, blon, blat)
    lut = np.zeros(ndep + 1, bool)
    lut[ids[real]] = True
    lut_z = np.full(ndep + 1, np.nan)
    lut_z[1:] = bot_z
    bottom = lut[dep_l] & (raw == lut_z[dep_l])      # the lowest cell, or the whole flat bottom of a real sink
    del dep_l, lut_z
    print(f"  depressions: {ndep}; >= 20 m {int((depth_d >= 20).sum())}; real (endorheic and >= 20 m) "
          f"{int(real.sum())}", flush=True)
    if not os.path.exists(fb):
        zin = raw.copy()
        zin[bottom] = NODATA
        with rasterio.open(fi, "w", **prof) as d:
            d.write(zin, 1)
        del zin
        assert wbt_ok(wbt.breach_depressions(fi, fb, fill_pits=True, callback=msgs.append), fb)
    with rasterio.open(fb) as src:
        zc = src.read(1).astype(np.float64)
    zc[bottom] = raw[bottom]                         # real sinks go back to FABDEM without conditioning
    zc[raw == NODATA] = NODATA
    nb_low = ndimage.minimum_filter(np.where(zc == NODATA, np.inf, zc), size=3)
    print(f"  breaching: {int((zc < raw - 1e-6).sum())} cells lowered; sink bottoms with a lower neighbour after "
          f"breaching: {int((bottom & (nb_low < zc - 1e-6)).sum())}", flush=True)
    del nb_low, bottom

    lat_c = n - (np.arange(nr) + 0.5) * RES
    dy = R_EARTH * math.radians(RES)
    dx = R_EARTH * np.cos(np.radians(lat_c)) * math.radians(RES)
    rcv, flat = d8_flats(zc, dx, dy, NODATA)
    nsink = int(((rcv == -1) & (raw != NODATA)).sum())
    del zc
    lab, ch, term, acc, nord = labels_from_flow(rcv, flat, (dx * dy).astype(np.float32))
    assert nord == nr * nc, "flow graph has a cycle"
    key, wcr = list_cracks(lab, ch, raw, np.float32(NODATA))
    ea, eb = endpoints(key, rcv)
    uniq, inv = np.unique(np.concatenate([ea, eb]), return_inverse=True)
    na, nb_ = inv[:key.size], inv[key.size:]
    del inv
    cor = uniq // 4
    cri, crj = cor // (nc + 1), cor % (nc + 1)
    nedge = (cri <= 1) | (crj <= 1) | (cri >= nr - 1) | (crj >= nc - 1)
    del cor, cri, crj
    order = np.argsort(-wcr, kind="stable")
    ci_, cd, chi, ctr = merge_passes(na, nb_, wcr, order, nedge, MIN_DEPTH)
    del na, nb_, uniq, order, nedge
    st, si, sd, shi, str_, sedge = sink_spills(key, wcr, term, rcv, flat)
    with rasterio.open(ff) as src:
        zfill = src.read(1)
    write_sinks(tag, term, rcv, flat, raw, zfill, n, w, dx, dy, area, core, north, st, si, sd, shi, str_, sedge, key,
                wcr)
    del zfill
    del term
    k = (sd >= MIN_DEPTH) & ~sedge
    si, sd, shi, str_ = si[k], sd[k], shi[k], str_[k]
    ci_, cd, chi = np.r_[ci_, si], np.r_[cd, sd], np.r_[chi, shi]
    cring = np.r_[np.zeros(len(ctr), bool), np.ones(len(si), bool)]
    ctr = np.r_[ctr, str_]
    ck, cw, cea, ceb = key[ci_], wcr[ci_], ea[ci_], eb[ci_]
    print(f"  cells {nr}x{nc}, sink cells {nsink}, divide cracks {key.size}, passes >= {MIN_DEPTH} m: {ck.size}",
          flush=True)

    cell = ck // 2
    ci, cj = cell // nc, cell % nc
    right = ck % 2 == 0
    lat = np.where(right, n - (ci + 0.5) * RES, n - (ci + 1) * RES)
    lon = np.where(right, w + (cj + 1) * RES, w + (cj + 0.5) * RES)
    i2, j2 = np.where(right, ci, ci + 1), np.where(right, cj + 1, cj)
    fa, fb2 = flat[ci, cj], flat[i2, j2]
    en_plano = (fa > 0) & (fa == fb2)
    flat_cells = np.bincount(flat[flat > 0].ravel(), minlength=flat.max() + 1)
    plano_km2 = np.where(en_plano, flat_cells[fa] * (dx[ci] * dy) / 1e6, np.nan)

    inarea = shapely.contains_xy(area, lon, lat)
    core_ = (lat <= north) & (lat > north - core)
    deep = inarea & (cd >= WIDTH_DEPTH)
    width = np.full(ck.size, np.nan)
    width[deep] = widths(lab, ch, raw, np.float32(NODATA), rcv, ck[deep], cea[deep], ceb[deep], cw[deep], dx, dy,
                         200_000)
    # each side followed downstream until it drains 5 km2, so the HydroBASINS label is read inside a basin
    tai, taj = trace_down(rcv, acc, ci, cj, 5e6)
    tbi, tbj = trace_down(rcv, acc, i2, j2, 5e6)

    import csv
    fo = os.path.join(OUT, f"passes_{tag}.csv")
    keep = inarea
    with open(fo, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["lat", "lon", "z_m", "depth_m", "upper_max_m", "width_m", "ring", "flat", "flat_km2", "core", "depth_lower_bound",
                     "lat_a", "lon_a", "lat_b", "lon_b", "lat_a5", "lon_a5", "lat_b5", "lon_b5"])
        f6 = lambda v: f"{v:.6f}"
        for k in np.flatnonzero(keep):
            wr.writerow([f6(lat[k]), f6(lon[k]), f"{cw[k]:.1f}", f"{cd[k]:.1f}", f"{chi[k]:.1f}",
                         "" if np.isnan(width[k]) else f"{width[k]:.0f}", int(cring[k]), int(en_plano[k]),
                         "" if np.isnan(plano_km2[k]) else f"{plano_km2[k]:.3f}", int(core_[k]), int(ctr[k]),
                         f6(n - (ci[k] + 0.5) * RES), f6(w + (cj[k] + 0.5) * RES),
                         f6(n - (i2[k] + 0.5) * RES), f6(w + (j2[k] + 0.5) * RES),
                         f6(n - (tai[k] + 0.5) * RES), f6(w + (taj[k] + 0.5) * RES),
                         f6(n - (tbi[k] + 0.5) * RES), f6(w + (tbj[k] + 0.5) * RES)])
    keep = inarea & core_
    print(f"  core passes: {keep.sum()}, with depth >= {WIDTH_DEPTH:.0f} m: {deep.sum()} -> {fo}", flush=True)
    # divide mask for review (a cell is marked if it borders a divide crack)
    dm = np.zeros((nr, nc), np.uint8)
    c = key // 2
    a_i, a_j = c // nc, c % nc
    dm[a_i, a_j] = 1
    dm[np.where(key % 2 == 0, a_i, a_i + 1), np.where(key % 2 == 0, a_j + 1, a_j)] = 1
    p2 = dict(prof, dtype="uint8", nodata=0)
    with rasterio.open(os.path.join(OUT, f"divisorias_{tag}.tif"), "w", **p2) as d:
        d.write(dm, 1)


if __name__ == "__main__":
    a = [float(x) for x in sys.argv[1:]]
    run(*a)
