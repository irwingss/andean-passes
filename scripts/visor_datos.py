"""Layers of the web viewer (docs/, GitHub Pages) of the Andean Passes Database.
The normalized current of the main connectivity run (BIO1, both readings, as in Fig. 2) becomes XYZ
tiles of indexed PNG; human modification and protected areas too; the cols with a depth of 100 m or more and the control
passes go to docs/capas/indice.json, column by column to keep it small.
python -P scripts/visor_datos.py [tiles]   (without "tiles" the tile pyramids are kept as they are)
"""
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyogrio
import rasterio
import shapely
from PIL import Image
from rasterio import features
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject

REPO = Path(__file__).resolve().parents[1]
INT, DB = REPO / 'data/interim', REPO / 'outputs/database'
RUN = INT / 'conectividad/principal'
DST = REPO / 'docs/capas'
ZMIN, ZMAX = 3, 9
U_MIN = 1e-6                                                     # as in 51_sintesis_conectividad.py

# normalized current in seven classes on a log scale; palette index 0 is always transparent
CUR = [('#2166AC', .85, 'below 0.5'), ('#67A9CF', .8, '0.5 to 0.8'), ('#D1E5F0', .7, '0.8 to 0.95'),
       ('#F7F7F7', .45, '0.95 to 1.05'), ('#FDDBC7', .7, '1.05 to 1.25'), ('#EF8A62', .8, '1.25 to 2'),
       ('#B2182B', .85, '2 or more')]
CUR_BREAKS = (0.5, 0.8, 0.95, 1.05, 1.25, 2)
PAL = {'trans_bio01': CUR, 'longi_bio01': CUR,
       'ghm': [('#fee8c8', .7, '0.1 to 0.2'), ('#fdbb84', .75, '0.2 to 0.4'), ('#d7301f', .8, '0.4 or more (high)')],
       'pa': [('#009E73', .6, 'Protected areas (WDPA, November 2026)')]}


def geoms(path, layer=None, columns=()):
    meta, _, wkb, campos = pyogrio.raw.read(path, layer=layer, columns=list(columns))
    return shapely.from_wkb(wkb), pd.DataFrame({k: v for k, v in zip(meta['fields'], campos)})


def tiles(nombre, idx, pal, aff, zmax=ZMAX):
    """XYZ tiles by nearest pixel. Plate carree to Web Mercator is separable: one lookup for rows and one for columns."""
    dst = DST / nombre
    shutil.rmtree(dst, ignore_errors=True)
    lut = [0, 0, 0] + [v for c, *_ in pal for v in (int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16))]
    alfa = bytes([0] + [round(255 * a) for _, a, _ in pal])
    alto, ancho = idx.shape
    oeste, norte = aff.c, aff.f
    este, sur = oeste + ancho * aff.a, norte + alto * aff.e
    merc = lambda lat: (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2  # noqa: E731
    n = 0
    for z in range(ZMIN, zmax + 1):
        m = 2 ** z
        for x in range(int((oeste + 180) / 360 * m), int((este + 180) / 360 * m) + 1):
            col = np.floor(((x + (np.arange(256) + .5) / 256) / m * 360 - 180 - oeste) / aff.a).astype(int)
            okc = (col >= 0) & (col < ancho)
            if not okc.any():
                continue
            for y in range(int(merc(norte) * m), int(merc(sur) * m) + 1):
                lat = np.degrees(np.arctan(np.sinh(math.pi * (1 - 2 * (y + (np.arange(256) + .5) / 256) / m))))
                fil = np.floor((lat - norte) / aff.e).astype(int)
                okf = (fil >= 0) & (fil < alto)
                if not okf.any():
                    continue
                t = idx[np.ix_(fil.clip(0, alto - 1), col.clip(0, ancho - 1))] * okf[:, None] * okc[None, :]
                if not t.any():
                    continue
                f = dst / f'{z}/{x}/{y}.png'
                f.parent.mkdir(parents=True, exist_ok=True)
                im = Image.fromarray(t.astype(np.uint8), 'P')
                im.putpalette(lut)
                im.save(f, optimize=True, transparency=alfa)
                n += 1
    print(nombre, n, 'tiles', flush=True)


def corriente(var, lect):
    """Normalized current of the main run (q 0.1, 200 km), undefined where the uniform current is below U_MIN."""
    with rasterio.open(RUN / f'teselas/{var}_{lect}_s200_env_normalizada.tif') as s:
        a, prof = s.read(1, masked=True, window=None).filled(np.nan), s.profile
    with rasterio.open(RUN / f'teselas/{var}_{lect}_s200_env_corriente.tif') as s:
        c = s.read(1, masked=True, window=None).filled(np.nan)
    with np.errstate(all='ignore'):
        return np.where(c / a >= U_MIN, a, np.nan).astype(np.float32), prof


def rasters():
    G = from_origin(-85, 14, 1 / 120, 1 / 120)                   # 30" geographic grid over the connectivity domain
    shape = (round(40 * 120), round(30 * 120))
    for lect, pref in (('baja', 'trans'), ('alta', 'longi')):
        a, prof = corriente('bio01', lect)
        for nombre, x in ((f'{pref}_bio01', a),):
            g = np.full(shape, np.nan, np.float32)
            reproject(x.astype(np.float32), g, src_transform=prof['transform'], src_crs=prof['crs'], src_nodata=np.nan,
                      dst_transform=G, dst_crs='EPSG:4326', dst_nodata=np.nan, resampling=Resampling.nearest)
            k = np.where(np.isfinite(g), np.digitize(g, CUR_BREAKS) + 1, 0).astype(np.uint8)
            tiles(nombre, k, CUR, G, 9)
    with rasterio.open(INT / 'M3_ghm_andes.tif') as r:
        g, ag = r.read(1, masked=True).filled(0), r.transform
    tiles('ghm', np.digitize(np.nan_to_num(g), (.1, .2, .4)).astype(np.uint8), PAL['ghm'], ag, 8)
    W, S, E, N, RES = -82, -56, -60, 13, 1 / 240
    pa, _ = geoms(INT / 'M3_wdpa_andes.gpkg')
    # ponytail: outlines drawn only as pixels, the WDPA itself is not redistributed
    tiles('pa', features.rasterize(((g, 1) for g in pa if g is not None), out_shape=(round((N - S) / RES),
          round((E - W) / RES)), transform=from_origin(W, N, RES, RES), dtype=np.uint8), PAL['pa'], from_origin(W, N, RES, RES))


def andes():
    a, _ = geoms(INT / 'M1_gmba_andes.gpkg')
    fc = {'type': 'FeatureCollection', 'features': [{'type': 'Feature', 'properties': {}, 'geometry':
          shapely.geometry.mapping(shapely.set_precision(shapely.simplify(g, .02), 1e-4))} for g in a]}
    (DST / 'andes.geojson').write_text(json.dumps(fc, separators=(',', ':')), encoding='utf-8')


def indice():
    c = pd.read_csv(DB / 'APD_cols.csv.gz', dtype={'trans_robust': str, 'longi_robust': str})
    c = c[c.depth_m >= 100].reset_index(drop=True)
    units = sorted(c.gmba.unique())
    ranks = ['continental', 'terminal', 'tributary']
    rob = {'above 1': 1, 'below 1': -1, 'depends': 0}

    def arr(s, d):
        return [None if not np.isfinite(v) else (round(float(v), d) if d else int(round(v))) for v in s]

    cols = {'id': c.apd_id.str[3:].astype(int).tolist(), 'lat': arr(c.lat, 4), 'lon': arr(c.lon, 4),
            'z': arr(c.z_m, 0), 'd': arr(c.depth_m, 0), 'w': arr(c.width_m, 0), 'iso': arr(c.isolation_km, 1),
            'u': [units.index(g) for g in c.gmba], 'r': [ranks.index(r) for r in c['rank']],
            't': arr(c.trans_med, 2), 'tmin': arr(c.trans_min, 2), 'tmax': arr(c.trans_max, 2),
            'l': arr(c.longi_med, 2), 'lmin': arr(c.longi_min, 2), 'lmax': arr(c.longi_max, 2),
            'tr': [rob.get(v) for v in c.trans_robust.fillna('')], 'lr': [rob.get(v) for v in c.longi_robust.fillna('')]}
    ctl = pd.read_csv(REPO / 'data/refs/pasos_control.csv')
    ctl = ctl[ctl.in_area == 1]
    out = {'units': units, 'ranks': ranks, 'pal': PAL, 'cols': cols,
           'control': [[r['name'], r.lat, r.lon] for _, r in ctl.iterrows()]}
    f = DST / 'indice.json'
    f.write_text(json.dumps(out, separators=(',', ':'), ensure_ascii=False), encoding='utf-8')
    print('indice.json', len(c), 'cols', round(f.stat().st_size / 1e6, 1), 'MB')


if __name__ == '__main__':
    DST.mkdir(parents=True, exist_ok=True)
    for old in ('atlas', 'miniaturas', 'dominios', 'umbrales', 'brechas.geojson', 'corredor.geojson',
                'divisoria.geojson', 'rutas.geojson', 'trans_med', 'longi_med'):                      # layers of the earlier census
        p = DST / old
        shutil.rmtree(p, ignore_errors=True) if p.is_dir() else p.unlink(missing_ok=True)
    andes()
    indice()
    if 'tiles' in sys.argv:
        rasters()
