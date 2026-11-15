"""58_puentes_links.py - The bridge of every barrier link of 55_barreras_grafo.py: the highest route over the relief
across the barrier between its two patches of montane habitat. The route may only use cells whose nearest patch
(straight-line distance) is one of the two, so it crosses the boundary between them, which is where the barrier line
lies, and cannot go around through other patches. Among those routes, the bridge route is the one whose lowest cell is
highest (maximin path over the elevation of insumos/elev.tif, 1 km); that lowest cell is the saddle where the montane
belt comes closest to continuous across the barrier, and its height is the height of the bridge. The patches are
rebuilt exactly as in 55 (same resistance, guild and labels).
The search runs inside the bounding box of the barrier (barreras_id.tif) widened by PAD km, or the whole grid if a
patch is not in it.
Writes data/interim/conectividad/<run>/barreras/puentes_links.csv (one row per barrier link of Table S25).
python -P scripts/58_puentes_links.py principal
"""
import importlib.util
import sys
import time

import numpy as np
import pandas as pd
import rasterio
from scipy import ndimage

RUN = sys.argv[1]
B = f"data/interim/conectividad/{RUN}/barreras"
PAD = 200
argv = sys.argv
sys.argv = [argv[0], RUN]
spec = importlib.util.spec_from_file_location("g55", "scripts/55_barreras_grafo.py")
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)
sys.argv = argv
from skimage.measure import label

Z = g.rd("elev")
ZF = np.where(np.isfinite(Z), Z, -1e4)
with rasterio.open(f"{B}/barreras_id.tif") as s:
    U = s.read(1)
EIGHT = np.ones((3, 3), bool)


def patches(v, cfg):
    h, q, amin = g.CFG[cfg]
    r = g.resistance(v, h, q)
    hab = g.guild(h)[0] & (r <= 1.0 + 1e-9)
    lab = label(hab, connectivity=2)
    n = np.bincount(lab.ravel())
    big = np.flatnonzero(n >= amin); big = big[big > 0]
    pat = np.zeros(n.size, np.int32); pat[big] = np.arange(1, big.size + 1)
    return pat[lab]


def window(bid):
    r, c = np.nonzero(U == bid)
    return (slice(max(r.min() - PAD, 0), r.max() + PAD + 1), slice(max(c.min() - PAD, 0), c.max() + PAD + 1))


def joined(z, a, b, t, ok):
    lab, _ = ndimage.label(((z >= t) | a | b) & ok, EIGHT)
    return np.intersect1d(np.unique(lab[a]), np.unique(lab[b])).size > 0, lab


def bridge(z, a, b, pw):
    _, (ir, ic) = ndimage.distance_transform_edt(pw == 0, return_indices=True)
    near = pw[ir, ic]                                         # nearest patch of every cell
    ok = (near == pw[a][0]) | (near == pw[b][0])
    lo, hi = -500.0, float(z[ok].max()) + 1                  # joined at lo (all land), apart above the highest cell
    if not joined(z, a, b, lo, ok)[0]:
        return None
    while hi - lo > 1:
        t = np.floor((lo + hi) / 2)
        if joined(z, a, b, t, ok)[0]:
            lo = t
        else:
            hi = t
    _, lab = joined(z, a, b, hi, ok)                              # apart at hi: the saddle cell joins the two sides
    la, lb = np.isin(lab, np.unique(lab[a])), np.isin(lab, np.unique(lab[b]))
    cand = ok & (z >= lo) & (z < hi) & ndimage.binary_dilation(la, EIGHT) & ndimage.binary_dilation(lb, EIGHT)
    if not cand.any():
        cand = ok & (z >= lo) & (z < hi) & ndimage.binary_dilation(la, EIGHT)
    r, c = np.argwhere(cand)[0]
    return lo, r, c


lk = pd.read_csv("outputs/supp/tables/S25_barrier_links.csv")
lk = lk[lk.barrier & (lk.barrier_id > 0)].reset_index()
rows = []
t0 = time.time()
for (var, cfg), e in lk.groupby(["var", "cfg"]):
    P = patches(int(var[3:]), cfg)
    for _, x in e.iterrows():
        w = window(x.barrier_id)
        a, b = P[w] == x.A, P[w] == x.B
        if not (a.any() and b.any()):
            w = (slice(None), slice(None))
            a, b = P == x.A, P == x.B
        res = bridge(ZF[w], a, b, P[w])
        if res is None:
            continue
        z, r, c = res
        r, c = r + (w[0].start or 0), c + (w[1].start or 0)
        lo, la = g.lonlat(r, c)
        rows.append({"index": x["index"], "barrier_id": x.barrier_id, "var": var, "cfg": cfg, "A": x.A, "B": x.B,
                     "bridge_z1km": z, "bridge_lon": float(lo), "bridge_lat": float(la)})
    print(var, cfg, len(e), "links", f"{time.time() - t0:.0f} s", flush=True)
pd.DataFrame(rows).to_csv(f"{B}/puentes_links.csv", index=False, float_format="%.5g")
