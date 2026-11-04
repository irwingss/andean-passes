"""Layers of the web viewer (docs/, GitHub Pages) of the Andean Passes Database.
Rasters become XYZ tiles of indexed PNG (Leaflet tile layers, same breaks and colours as the figures), the divide, the montane
gaps and the descent routes become simplified GeoJSON, and the cols go to docs/capas/indice.json with the sheets of File S4.
Run from anywhere: python visor_datos.py [tiles]   (without "tiles" the tile pyramids are kept as they are)
"""
import json, math, re, shutil, sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyogrio
import rasterio
import shapely
from PIL import Image
from rasterio import features
from rasterio.transform import from_origin
from rasterio.windows import from_bounds

REPO = Path(__file__).resolve().parent
L1 = REPO.parents[1]
INT, DB, TAB = L1 / 'data/interim', L1 / 'outputs/database', L1 / 'outputs/tables'
DST = REPO / 'docs/capas'

W, S, E, N, RES = -82, -56, -60, 13, 1 / 240                     # the 15" grid of M1_domains_15s.tif
AFF = from_origin(W, N, RES, RES)
SHAPE = (round((N - S) / RES), round((E - W) / RES))
ZMIN, ZMAX = 3, 9                                                # z9 is the first zoom finer than 15"

# palette index 0 is always transparent; colours as in the figures
PAL = {
    'dominios': [('#E69F00', .55, 'Pacific + Caribbean (W)'), ('#56B4E9', .55, 'Atlantic (E)'), ('#CC79A7', .55, 'Endorheic (I)')],
    'umbrales': [('#08306b', .8, 'Below 2000 m'), ('#2b7bba', .75, '2000 to 2500 m'), ('#a6cee3', .7, '2500 to 3000 m')],
    'ghm': [('#fee8c8', .7, '0.1 to 0.2'), ('#fdbb84', .75, '0.2 to 0.4'), ('#d7301f', .8, '0.4 or more (high)')],
    'pa': [('#009E73', .6, 'Protected areas (WDPA, November 2026)')],
}


def geoms(path, layer=None, columns=()):
    meta, _, wkb, campos = pyogrio.raw.read(path, layer=layer, columns=list(columns))
    return shapely.from_wkb(wkb), pd.DataFrame({k: v for k, v in zip(meta['fields'], campos)})


def tiles(nombre, idx, pal, aff=AFF, zmax=ZMAX):
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


def relieve():
    """Hillshades of the three focal windows as grey PNG in Web Mercator rows; returns name, file and Leaflet bounds."""
    out = []
    y = lambda lat: math.asinh(math.tan(math.radians(lat)))  # noqa: E731
    (DST / 'relieve').mkdir(exist_ok=True)
    for k, nombre in (('huancabamba', 'Huancabamba'), ('scolombia', 'Southern Colombia'), ('tachira', 'Táchira')):
        with rasterio.open(INT / f'M1_hillshade_{k}.tif') as r:
            a, b = r.read(1), r.bounds
        n = math.ceil(a.shape[0] * (y(b.top) - y(b.bottom)) / math.radians(b.top - b.bottom))
        lat = np.degrees(np.arctan(np.sinh(y(b.top) - (np.arange(n) + .5) * (y(b.top) - y(b.bottom)) / n)))
        fil = ((b.top - lat) / (b.top - b.bottom) * a.shape[0]).astype(int).clip(0, a.shape[0] - 1)
        im = Image.fromarray((np.nan_to_num(a[fil]).clip(0, 1) * 255).astype(np.uint8), 'L')
        im.quantize(32).save(DST / f'relieve/{k}.png', optimize=True)
        out.append({'nombre': nombre, 'png': f'relieve/{k}.png', 'bounds': [[b.bottom, b.left], [b.top, b.right]]})
    return out


def rasters():
    dom, c = geoms(INT / 'M1_domains.gpkg', columns=['dom'])
    d = features.rasterize(((g, 'WEI'.index(v) + 1) for g, v in zip(dom, c['dom'])), out_shape=SHAPE, transform=AFF, dtype=np.uint8)
    tiles('dominios', d, PAL['dominios'])

    with rasterio.open(L1 / 'data/raw/hyd_sa_dem_15s/hyd_sa_dem_15s.tif') as r:
        z = r.read(1, window=from_bounds(W, S, E, N, r.transform))
    assert z.shape == SHAPE
    k = np.select([z < 2000, z < 2500, z < 3000], [1, 2, 3], 0).astype(np.uint8)
    k[d == 0] = 0                                                # footprint of the drainage domains
    tiles('umbrales', k, PAL['umbrales'])

    with rasterio.open(INT / 'M3_ghm_andes.tif') as r:
        g, ag = r.read(1, masked=True).filled(0), r.transform
    tiles('ghm', np.digitize(np.nan_to_num(g), (.1, .2, .4)).astype(np.uint8), PAL['ghm'], ag, 8)

    pa, _ = geoms(INT / 'M3_wdpa_andes.gpkg')
    # ponytail: outlines drawn only as pixels, the WDPA itself is not redistributed
    tiles('pa', features.rasterize(((g, 1) for g in pa if g is not None), out_shape=SHAPE, transform=AFF, dtype=np.uint8), PAL['pa'])


def geojson(nombre, feats):
    """feats: (shapely geometry, properties). Coordinates to 4 decimals (11 m)."""
    fc = {'type': 'FeatureCollection', 'features': [{'type': 'Feature', 'properties': p,
          'geometry': shapely.geometry.mapping(shapely.set_precision(g, 1e-4))} for g, p in feats if not g.is_empty]}
    (DST / nombre).write_text(json.dumps(fc, separators=(',', ':'), ensure_ascii=False), encoding='utf-8')
    print(nombre, len(fc['features']), round((DST / nombre).stat().st_size / 1e6, 2), 'MB', flush=True)


def tramos(p, marca):
    """Runs of consecutive bins with marca True: one line per run, split where the divide itself jumps."""
    grupo = (marca != marca.shift()).cumsum()
    for _, r in p[marca].groupby(grupo[marca]):
        partes = np.cumsum(r.gap_before_km.to_numpy() > 0)
        lineas = [r[partes == i][['lon', 'lat']].to_numpy() for i in np.unique(partes)]
        lineas = [l if len(l) > 1 else np.repeat(l, 2, axis=0) for l in lineas]
        yield shapely.MultiLineString([shapely.LineString(l) for l in lineas]), r


def vectores(gw):
    d, _ = geoms(DB / 'APD_divide.gpkg')
    ie, c = geoms(INT / 'M1_divide_raw.gpkg', columns=['type'])
    geojson('divisoria.geojson', [(shapely.simplify(g, .002), {'tipo': 'W|E'}) for g in d] +
            [(shapely.simplify(g, .002), {'tipo': 'I|E'}) for g, t in zip(ie, c['type']) if t == 'I|E'])
    a, _ = geoms(INT / 'M1_gmba_andes.gpkg')
    geojson('andes.geojson', [(shapely.simplify(g, .02), {}) for g in a])

    prof = pd.read_csv(INT / 'M1_profile.csv').sort_values('bin').reset_index(drop=True)
    tabla = pd.read_csv(INT / 'M1_montane_gaps.csv')
    fb = []
    for th in (2000, 2500, 3000):
        n = suma = 0
        for g, r in tramos(prof, prof.elev_cross < th):
            km = r.d_km.iloc[-1] - r.d_km.iloc[0] + .5
            fb.append((g, {'umbral': th, 'km': round(km, 1), 'min': round(r.elev_cross.min()), 'lat': [round(r.lat.iloc[0], 2), round(r.lat.iloc[-1], 2)]}))
            n += 1
            suma += km
        t = tabla[tabla.threshold_m == th]
        assert n == len(t) and abs(suma - t.length_km.sum()) < 1, th   # same runs as 04_collados.R
    geojson('brechas.geojson', fb)

    b = pd.read_csv(DB / 'APD_montane_corridor.csv').sort_values('bin').reset_index(drop=True)
    b['gap_before_km'] = (b.d_km.diff().fillna(0) > 5) * 1.0     # same jumps of the ruler as the profile
    clase = np.select([b.gap_now, b.loss_robust_370_2090], ['gap', 'lost'], 'kept')
    fc = []
    for k in ('kept', 'lost', 'gap'):
        for g, r in tramos(b, pd.Series(clase == k)):
            fc.append((g, {'clase': k, 'km': round(r.d_km.iloc[-1] - r.d_km.iloc[0] + .5, 1), 'pa': round(100 * r.in_PA.mean()), 'ghm': round(r.ghm.mean(), 2)}))
    geojson('corredor.geojson', fc)

    rt, c = geoms(DB / 'APD_flanks.gpkg', 'flank_routes', ['col_id', 'side_lab'])
    geojson('rutas.geojson', [(shapely.simplify(g, .001), {'col': i, 'lado': s}) for g, i, s in zip(rt, c['col_id'], c['side_lab']) if i in gw])


def num(x, d=2):
    return round(float(x), d) if np.isfinite(x) else None


def pasos():
    c = pd.read_csv(DB / 'APD_cols.csv')
    h3 = pd.read_csv(TAB / 'M2_H3_tipologia.csv').set_index('col_id')
    t3 = pd.read_csv(TAB / 'Table3_conservation_priority.csv').set_index('col_id')
    con = pd.read_csv(DB / 'APD_conservation_gateways.csv').set_index('col_id')
    fut = pd.read_csv(DB / 'APD_climate_future.csv').query("ssp == 370 and period == '2081-2100'").set_index('col_id')
    atlas = {}
    (DST / 'atlas').mkdir(exist_ok=True)
    (DST / 'miniaturas').mkdir(exist_ok=True)
    for f in sorted((L1 / 'outputs/supp/atlas').glob('S4_*.png')):
        i = re.search(r'APD\d+', f.name)[0]
        shutil.copy(f, DST / f'atlas/{i}.png')
        with Image.open(f) as im:                                # 2x the popup width, for phones
            assert im.size == (2007, 1417), f                    # the popup reserves this aspect ratio
            im.convert('RGB').resize((540, 381), Image.LANCZOS).save(DST / f'miniaturas/{i}.jpg', quality=75, optimize=True)
        atlas[i] = f'atlas/{i}.png'
    pts, cp = geoms(DB / 'APD_flanks.gpkg', 'flank_points', ['col_id', 'side_lab', 'flank_m'])
    fl = {}
    for g, i, s, m in zip(pts, cp['col_id'], cp['side_lab'], cp['flank_m']):
        if m == 1000 and g is not None and not shapely.is_empty(g):
            fl.setdefault(i, []).append([round(g.y, 4), round(g.x, 4), s])
    out = []
    for r in c.itertuples():
        e = {'id': r.col_id, 'nombre': r.name if isinstance(r.name, str) else None, 'lat': r.lat, 'lon': r.lon, 'alt': round(r.elev_cross),
             'prof': round(r.depth), 'iso': num(r.iso_km, 0), 'ancho': num(r.width_km, 1), 'grupo': 'tropical' if r.tropical else 'temperate'}
        if r.gateway:
            e['puerta'] = 1
        if r.terminal:
            e['terminal'] = 1
        if r.col_id in fl and r.gateway:
            e['flancos'] = fl[r.col_id]
        if r.col_id in h3.index:                                 # the 30 tropical gateways of Table 3
            a, b, k, u = h3.loc[r.col_id], t3.loc[r.col_id], con.loc[r.col_id], fut.loc[r.col_id]
            e['nombre'] = e['nombre'] or b['name']
            e |= {'pais': b.country, 'tipo': a.h3_type, 'to': num(a.TO_ann), 'ai': num(a.AI_gap), 'dto': f'{u["dTO.mean"]:+.2f} [{u["dTO.min"]:+.2f}, {u["dTO.max"]:+.2f}]',
                  'dpg': f'{u["dPg.mean"]:+.2f} [{u["dPg.min"]:+.2f}, {u["dPg.max"]:+.2f}]', 'frac': num(k.frac_dTstar_370_2090), 'pa5': num(k.pct_PA_5km, 0),
                  'ghm': num(k.ghm_crest), 'loss': num(k.montane_loss_km, 0), 'critica': bool(k.critical), 'crest_pa': bool(k.crest_in_PA), 'pres': bool(b.high_pressure)}
        if r.col_id in atlas:
            e['atlas'] = atlas[r.col_id]
        out.append(e)
    return out, set(c.col_id[c.gateway])


if __name__ == '__main__':
    DST.mkdir(parents=True, exist_ok=True)
    cols, gw = pasos()
    vectores(gw)
    if 'tiles' in sys.argv:
        rasters()
    (DST / 'indice.json').write_text(json.dumps({'pal': PAL, 'relieve': relieve(), 'cols': cols}, separators=(',', ':'), ensure_ascii=False, default=lambda x: x.item()), encoding='utf-8')
