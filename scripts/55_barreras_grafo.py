"""55_barreras_grafo.py - Barriers for montane climates as links of a patch-based graph (minimum planar graph; Fall et
al. 2007; Galpern, Manseau & Fall 2011), one variable at a time, on the resistances of the connectivity (C4).
  Habitat: montane-guild cells (colder than T2000(lat), C3) inside the montane envelope of the variable (resistance 1).
  The montane limit is 2000 m, with 1500 and 2500 m as in D9 and 500, 1000, 3000, 3500, 4000 and 4500 m to look for barriers at
  other heights: for those the guild, its envelopes and the resistances are rebuilt here with the recipe of 47_conectividad_insumos.py (T_h(lat) from the cells within 100 m of h).
  Patches: connected habitat of AMIN km2 or more.
  Every cell is allocated to the patch it reaches at least cost (cost-distance Voronoi over the resistance); two
  patches that touch in that allocation are joined by a link, whose least-cost path crosses their common boundary.
  The excess cost of a link is its cost minus the length of its path: zero if the path stays inside the envelope, so
  the variable sets no climatic barrier there. A link is a barrier for the variable when its excess cost exceeds EMIN
  (1 km-equivalent), and its barrier line is the boundary between the two patches in the allocation.
  The cell of highest resistance along the least-cost path of a link is the floor of the gap it crosses (gap_lon,
  gap_lat, gap_r), which 57_barreras_pasos.py joins to the census passes.
Configurations: montane limit (500 to 4500 m every 500 m), margin q (0.1, 0.05, 0.2) and minimum patch area (500, 100,
1000 km2), each varied alone from the main one (2000 m, q 0.1, 500 km2).
Writes data/interim/conectividad/<run>/barreras/links_<var>_<cfg>.csv (every link with its cost, length and excess)
       data/interim/conectividad/<run>/barreras/lineas_<var>_<cfg>.tif (uint8, 1 on the lines of the barriers)
python -P scripts/55_barreras_grafo.py principal [vars] [cfgs]       e.g. principal "1 12" main
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from skimage.graph import MCP_Geometric
from skimage.measure import label

RUN = sys.argv[1]
D = f"data/interim/conectividad/{RUN}"
I, OUT = f"{D}/insumos", f"{D}/barreras"
os.makedirs(OUT, exist_ok=True)
CFG = {"main": (2000, 0.1, 500), "q0.05": (2000, 0.05, 500), "q0.2": (2000, 0.2, 500), "a100": (2000, 0.1, 100),
       "a1000": (2000, 0.1, 1000), "h500": (500, 0.1, 500), "h1000": (1000, 0.1, 500), "h1500": (1500, 0.1, 500), "h2500": (2500, 0.1, 500),
       "h3000": (3000, 0.1, 500), "h3500": (3500, 0.1, 500), "h4000": (4000, 0.1, 500),
       "h4500": (4500, 0.1, 500)}
EMIN = 1.0
VARS = [int(v) for v in sys.argv[2].split()] if len(sys.argv) > 2 and sys.argv[2] != "all" else list(range(1, 20))
CFGS = sys.argv[3].split() if len(sys.argv) > 3 else list(CFG)

with rasterio.open(f"{I}/fuente_alta.tif") as s:
    MON2000 = s.read(1, masked=True).filled(0) == 1
    PROF, T, CRS = s.profile, s.transform, s.crs
H, W = MON2000.shape
tr = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
IQR = pd.read_csv(f"{I}/margenes.csv").set_index("var").iqr
_cache = {}


def rd(n):
    with rasterio.open(f"{I}/{n}.tif") as s:
        return s.read(1, masked=True).filled(np.nan)


def grid_lat():
    if "lat" not in _cache:
        xs = T.c + (np.arange(W) + 0.5) * T.a
        ys = T.f + (np.arange(H) + 0.5) * T.e
        _cache["lat"] = tr.transform(*np.meshgrid(xs, ys))[1].astype(np.float32)
    return _cache["lat"]


def guild(h):
    """Montane guild for limit h (m) and what the envelopes need, as in 47_conectividad_insumos.py."""
    if h == 2000:
        return MON2000, None
    if h not in _cache:
        lat, elev, b1 = grid_lat(), rd("elev"), rd("bio01")
        valid = np.isfinite(elev) & np.isfinite(b1)
        sel = valid & (np.abs(elev - h) <= 100)
        lc, tc = lat[sel], b1[sel]
        bins = np.round(np.arange(np.floor(lat[valid].min() * 10), np.ceil(lat[valid].max() * 10) + 1) / 10, 1)
        th = np.full(bins.size, np.nan)
        for i, b in enumerate(bins):
            m = np.abs(lc - b) <= 1
            if m.sum() >= 200:
                th[i] = np.median(tc[m])
        ok = np.flatnonzero(np.isfinite(th))
        th = th[ok[np.abs(ok[:, None] - np.arange(bins.size)).argmin(0)]]
        mon = valid & (b1 <= np.interp(lat, bins, th))
        order = np.argsort(lat[mon])
        _cache[h] = (mon, (bins, order, lat[mon][order]))
    return _cache[h]


def resistance(v, h, q):
    """Montane-guild resistance of variable v for limit h and margin q (C4)."""
    if h == 2000:
        return rd(f"res_envalta_bio{v:02d}" + ("" if q == 0.1 else f"_q{q:g}")).astype(np.float64)
    mon, (bins, order, ls) = guild(h)
    a = rd(f"bio{v:02d}")
    x = a[mon][order]
    lo, hi = np.searchsorted(ls, bins - 1, "left"), np.searchsorted(ls, bins + 1, "right")
    e = np.full((bins.size, 2), np.nan)
    for i in range(bins.size):
        if hi[i] - lo[i] >= 200:
            e[i] = np.percentile(x[lo[i]:hi[i]], [5, 95])
    ok = np.flatnonzero(np.isfinite(e[:, 0]))
    e = e[ok[np.abs(ok[:, None] - np.arange(bins.size)).argmin(0)]]
    lat = grid_lat()
    p5, p95 = np.interp(lat, bins, e[:, 0]), np.interp(lat, bins, e[:, 1])
    d = np.maximum(np.maximum(p5 - a, a - p95), 0)
    return (1 + d / (q * IQR[v])).astype(np.float64)


def lonlat(r, c):
    return tr.transform(T.c + (np.asarray(c) + 0.5) * T.a, T.f + (np.asarray(r) + 0.5) * T.e)


def run(v, cfg):
    h, q, amin = CFG[cfg]
    r = resistance(v, h, q)
    hab = guild(h)[0] & (r <= 1.0 + 1e-9)
    lab = label(hab, connectivity=2)
    n = np.bincount(lab.ravel())
    big = np.flatnonzero(n >= amin); big = big[big > 0]
    pat = np.zeros(n.size, np.int32); pat[big] = np.arange(1, big.size + 1)
    P = pat[lab]                                           # patch id, 0 outside the patches
    cost_arr = np.where(np.isfinite(r), r, np.inf)
    mcp = MCP_Geometric(cost_arr)
    starts = np.argwhere(P > 0)
    cum, tb = mcp.find_costs(starts)
    off = np.asarray(mcp.offsets).astype(np.int64)
    # predecessor of every cell along its least-cost path (traceback holds the step that led into the cell)
    flat = np.arange(H * W)
    tbf = tb.ravel()
    par = flat.copy()
    m = tbf >= 0
    par[m] = flat[m] - (off[tbf[m], 0] * W + off[tbf[m], 1])
    lenstep = np.zeros(H * W); lenstep[m] = np.hypot(off[tbf[m], 0], off[tbf[m], 1])
    # allocation by pointer doubling: every cell takes the patch at the start of its least-cost path
    alloc = P.ravel().copy()
    jump = par.copy()
    for _ in range(40):
        nxt = jump[jump]
        if np.array_equal(nxt, jump):
            break
        jump = nxt
    alloc = alloc[jump]
    A = alloc.reshape(H, W)
    cum = np.where(np.isfinite(cum), cum, np.nan)
    rows = []
    lines = np.zeros((H, W), np.uint8)
    contacts = []
    Ap = np.pad(A, 1)
    for dy, dx in ((0, 1), (1, 0), (1, 1), (1, -1)):
        nb = Ap[1 + dy:1 + dy + H, 1 + dx:1 + dx + W]
        ra1, ca1 = np.nonzero((A > 0) & (nb > 0) & (A != nb))
        ra2, ca2 = ra1 + dy, ca1 + dx
        c1, c2 = cum[ra1, ca1], cum[ra2, ca2]
        e = 0.5 * (cost_arr[ra1, ca1] + cost_arr[ra2, ca2]) * np.hypot(dy, dx)
        contacts.append(pd.DataFrame({"a": A[ra1, ca1], "b": A[ra2, ca2], "r1": ra1, "c1": ca1, "r2": ra2,
                                      "c2": ca2, "cost": c1 + c2 + e, "step": np.hypot(dy, dx)}))
    ct = pd.concat(contacts, ignore_index=True)
    ct = ct[np.isfinite(ct.cost)]
    ct["A"], ct["B"] = np.minimum(ct.a, ct.b), np.maximum(ct.a, ct.b)
    best = ct.loc[ct.groupby(["A", "B"]).cost.idxmin()]
    area = np.bincount(P.ravel(), minlength=big.size + 1)
    rf = cost_arr.ravel()
    for _, x in best.iterrows():
        L = x.step
        gap = int(x.r1) * W + int(x.c1)                   # the cell of highest resistance along the path: the gap floor
        for rr, cc in ((int(x.r1), int(x.c1)), (int(x.r2), int(x.c2))):
            i = rr * W + cc
            while True:
                if rf[i] > rf[gap]:
                    gap = i
                if par[i] == i:
                    break
                L += lenstep[i]; i = par[i]
        lo, la = lonlat(int(x.r1), int(x.c1))
        glo, gla = lonlat(gap // W, gap % W)
        rows.append({"var": f"bio{v:02d}", "cfg": cfg, "A": int(x.A), "B": int(x.B), "area_A": int(area[int(x.A)]),
                     "area_B": int(area[int(x.B)]), "cost": x.cost, "length_km": L, "excess": x.cost - L,
                     "lon": float(lo), "lat": float(la), "gap_lon": float(glo), "gap_lat": float(gla),
                     "gap_r": float(rf[gap])})
    lk = pd.DataFrame(rows)
    lk["barrier"] = lk.excess > EMIN
    bar = set(zip(lk.A[lk.barrier], lk.B[lk.barrier]))
    sel = np.array([p in bar for p in zip(ct.A, ct.B)])
    lines[ct.r1.values[sel], ct.c1.values[sel]] = 1
    lines[ct.r2.values[sel], ct.c2.values[sel]] = 1
    lk.to_csv(f"{OUT}/links_bio{v:02d}_{cfg}.csv", index=False, float_format="%.4g")
    p = PROF.copy(); p.update(dtype="uint8", nodata=0, compress="deflate", count=1)
    with rasterio.open(f"{OUT}/lineas_bio{v:02d}_{cfg}.tif", "w", **p) as o:
        o.write(lines, 1)
    return big.size, len(lk), int(lk.barrier.sum())


if __name__ == "__main__":
    for cfg in CFGS:
        for v in VARS:
            t0 = time.time()
            np_, nl, nb = run(v, cfg)
            print(f"bio{v:02d} {cfg}: {np_} patches, {nl} links, {nb} barriers, {time.time() - t0:.0f} s", flush=True)
