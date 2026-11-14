"""50_conectividad_teselas.py - Climatic connectivity by tiled wall-to-wall circuits (Pelletier et al. 2014), on the
inputs of 47_conectividad_insumos.py. Runs on a laptop: one sparse solve per tile and direction instead of one per
source block.
  Transverse reading (lowland species): current from the west edge to the east edge of each tile, injected and
  grounded only at lowland cells (warmer than T2000(lat)).
  Longitudinal reading (montane species): current from the north edge to the south edge, only at montane cells.
Longitudinal tiles are squares of side S. Transverse tiles are S km from north to south and, from west to east, S km
or the smallest of 2S, 3S, ... up to LMAX = 800 km that puts lowland cells on both walls: a tile grows only where the
cordillera is wider than S, as on the Altiplano, so it spans one cordillera from the lowland on one side to the lowland
on the other and never several. Longitudinal tiles do not grow, since a longer one would only add lowland gaps that a
montane species does not cross. Tile centres are S/4 apart; currents are summed over tiles and divided by the same sum with uniform resistance, so the
normalized current is above 1 where climate concentrates movement and below 1 where it hinders it.
Each tile enters the sums with a Hann window over its own size, highest at its centre and near zero at its edges, so a
cell takes most of its value from the tiles centred near it and no step remains where a tile begins or ends; the edges,
where current is injected and grounded, weigh least.
python -P scripts/50_conectividad_teselas.py <run> <side_km> <vars|all> <shapes>     e.g. principal 200 "1" env
The shape "env" reads res_env{baja,alta}_bioNN.tif (first margin q of 47, 0.1); "env_q0.05" reads
res_env{baja,alta}_bioNN_q0.05.tif.
Writes data/interim/conectividad/<run>/teselas/<var>_<reading>_s<side>_<shape>_{corriente,normalizada}.tif and
teselas_s<side>_<var>[_q<q>].csv (relative conductance of every tile, C7), e.g. bio01_alta_s200_env_q0.05_corriente.tif
and teselas_s200_bio01_q0.05.csv; the q suffix appears only for shapes that carry it.
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import rasterio
from rasterio.windows import Window
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

RUN = sys.argv[1]
S = int(sys.argv[2])
INS = f"data/interim/conectividad/{RUN}/insumos"
OUT = f"data/interim/conectividad/{RUN}/teselas"
NB = [(0, 1, 1.0), (1, -1, 2 ** 0.5), (1, 0, 1.0), (1, 1, 2 ** 0.5)]      # each neighbour pair once; distance in km
LMAX = 800                                                                 # longest tile along the current, km
SIZES = [S * k for k in range(1, max(LMAX // S, 1) + 1)]                   # tried in this order


def hann(h, w):
    """Tile weight: 2-D Hann window over the tile; with tiles of side S every S/4 the windows of a cell add to 16."""
    a, b = (2 * np.sin(np.pi * (np.arange(n) + 0.5) / n) ** 2 for n in (h, w))
    return np.outer(a, b)


def read(n):
    with rasterio.open(f"{INS}/{n}.tif") as s:
        return s.read(1, masked=True).filled(np.nan), s.profile


def solve(r, wall0, wall1):
    """Current per cell between two walls (Dirichlet 1 and 0) on resistance r (nan = outside). None if no path."""
    h, w = r.shape
    ok = np.isfinite(r)
    idx = np.full(r.shape, -1)
    idx[ok] = np.arange(ok.sum())
    P = np.pad(idx, 1, constant_values=-1)
    Rp = np.pad(r, 1, constant_values=np.nan)
    ia, ra = P[1:-1, 1:-1], Rp[1:-1, 1:-1]
    I, J, C = [], [], []
    for dy, dx, d in NB:
        ib, rb = P[1 + dy:h + 1 + dy, 1 + dx:w + 1 + dx], Rp[1 + dy:h + 1 + dy, 1 + dx:w + 1 + dx]
        m = (ia >= 0) & (ib >= 0)
        I.append(ia[m]); J.append(ib[m]); C.append(2 / (d * (ra[m] + rb[m])))
    I, J, C = map(np.concatenate, (I, J, C))
    n = ok.sum()
    G = sparse.coo_matrix((C, (I, J)), (n, n)).tocsr()
    G = G + G.T
    # keep only the components that touch both walls
    w0, w1 = idx[wall0 & ok], idx[wall1 & ok]
    _, lab = connected_components(G, directed=False)
    both = np.intersect1d(np.unique(lab[w0]), np.unique(lab[w1]))
    if both.size == 0:
        return None
    keep = np.isin(lab, both)
    v = np.zeros(n)
    fixed = np.zeros(n, bool)
    fixed[w0] = fixed[w1] = True
    v[w0] = 1.0
    free = keep & ~fixed
    L = sparse.diags(np.asarray(G.sum(1)).ravel()) - G
    Lf = L[free][:, free].tocsc()
    v[free] = spsolve(Lf, -L[free][:, fixed] @ v[fixed])
    v[~keep] = np.nan
    cur = np.zeros(n)
    dv = np.abs(v[I] - v[J]) * C
    dv[np.isnan(dv)] = 0
    np.add.at(cur, I, dv); np.add.at(cur, J, dv)
    total = (L @ np.nan_to_num(v))[w0].sum()          # current leaving the source wall
    out = np.full(r.shape, np.nan)
    out[ok] = cur / 2 / total                          # one ampere per tile: the map shows redistribution only
    return out, total


def walls(guild, axis):
    """Per row (axis 1, west-east) or column (axis 0, north-south): the outermost guild cell at each end, if it lies in
    the outer quarter of the tile. A coast or a domain edge inside the tile acts as its edge."""
    g = ~np.isnan(guild) & (guild == 1)
    gg = g if axis == 1 else g.T
    n = gg.shape[1]
    w0, w1 = np.zeros_like(gg), np.zeros_like(gg)
    rows = np.flatnonzero(gg.any(1))
    first = gg.argmax(1)
    last = n - 1 - gg[:, ::-1].argmax(1)
    for i in rows:
        if first[i] < n // 4:
            w0[i, first[i]] = True
        if last[i] >= n - n // 4:
            w1[i, last[i]] = True
    return (w0, w1) if axis == 1 else (w0.T, w1.T)


def tile(args):
    """One tile of one reading around a centre: the shortest size along the current that has walls on both sides, then
    one solve per resistance. Each worker reads only the window of its tile, so memory does not grow with the domain."""
    lect, cy, cx = args
    axis = 1 if lect == "baja" else 0
    for L in (SIZES if axis == 1 else SIZES[:1]):
        h, w = (S, L) if axis == 1 else (L, S)
        y0, x0 = max(cy - h // 2, 0), max(cx - w // 2, 0)
        y1, x1 = min(cy - h // 2 + h, H_), min(cx - w // 2 + w, W_)
        win = Window(x0, y0, x1 - x0, y1 - y0)
        w0, w1 = walls(readw(f"fuente_{lect}", win), axis)
        if w0.sum() >= 3 and w1.sum() >= 3:
            break
    else:
        return lect, None, {}
    mask = readw(F_["baja"][0], win)
    uni = solve(np.where(np.isnan(mask), np.nan, 1.0), w0, w1)
    if uni is None:
        return lect, None, {}
    res = {"uni": uni}
    for k, n in enumerate(F_[lect]):
        c = solve(readw(n, win), w0, w1)
        if c is not None:
            res[k] = (c[0], c[1] / uni[1])                  # conductance of the tile relative to uniform
    return lect, (y0, x0, y1 - y0, x1 - x0), res


def readw(n, win):
    with rasterio.open(f"{INS}/{n}.tif") as s:
        return s.read(1, window=win, masked=True).filled(np.nan)


def init(files, hw):
    global F_, H_, W_
    F_ = files
    H_, W_ = hw


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    vars_ = list(range(1, 20)) if sys.argv[3] == "all" else [int(v) for v in sys.argv[3].split()]
    shapes = sys.argv[4].split()
    names = [f"bio{v:02d}_{s}" for v in vars_ for s in shapes]
    # "env" resistances depend on the reading (res_envbaja_*, res_envalta_*); the others are shared
    def rfile(v, s, l):
        base, _, q = s.partition("_")                 # "env_q0.05" -> res_envbaja_bio01_q0.05
        return f"res_{base}{l if base == 'env' else ''}_bio{v:02d}" + (f"_{q}" if q else "")
    files = {l: [rfile(v, s, l) for v in vars_ for s in shapes] for l in ("baja", "alta")}
    r0, prof = read(files["baja"][0])
    H, W = r0.shape
    step = S // 4
    ys = sorted(set(list(range(0, max(H - S, 0) + 1, step)) + [max(H - S, 0)]))
    xs = sorted(set(list(range(0, max(W - S, 0) + 1, step)) + [max(W - S, 0)]))
    cen = [(y + S // 2, x + S // 2) for y in ys for x in xs
           if np.isfinite(r0[y:y + S, x:x + S]).sum() > S * S // 10]
    jobs = [(l, cy, cx) for l in ("baja", "alta") for cy, cx in cen]
    print(f"{len(cen)} tile centres of {S} km (up to {SIZES[-1]} km along the current), {len(names)} resistances, "
          "both readings")
    t0 = time.time()
    tiles = []
    acc = {(l, k): np.zeros((H, W)) for l in ("baja", "alta") for k in list(range(len(names))) + ["uni"]}
    with Pool(int(os.environ.get("NPROC", os.cpu_count())), initializer=init, initargs=(files, (H, W))) as p:
        for lect, box, res in p.imap_unordered(tile, jobs, chunksize=2):
            if box is None:
                continue
            y, x, h, w = box
            wt = hann(h, w)
            for key, (c, g) in res.items():
                acc[(lect, key)][y:y + h, x:x + w] += np.nan_to_num(c) * wt
                if key != "uni":
                    tiles.append((y + h // 2, x + w // 2, h, w, lect, names[key], g))
    print(f"solved in {time.time() - t0:.0f} s")
    import pandas as pd
    from pyproj import Transformer
    tt = pd.DataFrame(tiles, columns=["row", "col", "alto_km", "ancho_km", "lectura", "resistencia", "conductancia_rel"])
    xx, yy = rasterio.transform.xy(prof["transform"], tt.row, tt.col)
    tt["lon"], tt["lat"] = Transformer.from_crs(prof["crs"], "EPSG:4326", always_xy=True).transform(xx, yy)
    for nm, g in tt.groupby("resistencia"):        # one file per variable and q, so calls do not overwrite
        v, s = nm.split("_", 1)
        g.to_csv(f"{OUT}/teselas_s{S}_{v}{s[len(s.split('_')[0]):]}.csv", index=False, float_format="%.5g")
    prof.update(dtype="float32", nodata=-9999.0)
    valid = np.isfinite(r0)
    for lect in ("baja", "alta"):
        u = acc[(lect, "uni")]
        for k, nm in enumerate(names):
            v, s = nm.split("_", 1)
            c = acc[(lect, k)]
            for kind, a in (("corriente", c), ("normalizada", np.where(u > 0, c / np.where(u > 0, u, 1), np.nan))):
                with rasterio.open(f"{OUT}/{v}_{lect}_s{S}_{s}_{kind}.tif", "w", **prof) as d:
                    d.write(np.where(valid & np.isfinite(a), a, -9999.0).astype(np.float32), 1)
    print("->", OUT)
