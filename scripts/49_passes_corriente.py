"""49_passes_corriente.py - Reads the tiled-circuit maps of one connectivity run (50_conectividad_teselas.py) at every pass
of the census and builds the latitudinal profile of tile relative conductance (PRD C7, C8). The control C9 reads only
the tiles that cross the divide, so it is in 51_sintesis_conectividad.py.
python -P scripts/49_passes_corriente.py <run> <lat_north> <lat_south>
Writes, under data/interim/conectividad/<run>:
  passes_corriente.csv.gz   one row per pass, variable, reading, margin q and tile side (pass_id is the row of
                          data/interim/censo/passes_andes.csv): corriente and normalizada
  perfil_conductancia.csv median relative conductance of the tiles per variable, reading, q, side and 0.5 degree band
  mapa_<var>_s<side>[_q<q>].png  normalized current of both readings with the control passes: every variable for
                          q = 0.1 and 200 km, bio01 for the other configurations
q and side come from the file names of 50 (no q suffix means q = 0.1, the first margin of 47).
"""
import glob
import gzip
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer

RUN, LN, LS = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
D = f"data/interim/conectividad/{RUN}"
T = f"{D}/teselas"
Q0 = 0.1       # q of the unsuffixed env files (first margin of 47)
passes = pd.read_csv("data/interim/censo/passes_andes.csv", usecols=["lat", "lon"])
passes["pass_id"] = np.arange(len(passes))
passes = passes[(passes.lat <= LN) & (passes.lat >= LS)].copy()
ctl = pd.read_csv("data/refs/pasos_control.csv")
ctl = ctl[(ctl.lat <= LN) & (ctl.lat >= LS)]
with rasterio.open(f"{D}/insumos/elev.tif") as r:      # every map shares this grid
    to = Transformer.from_crs("EPSG:4326", r.crs, always_xy=True)
    rr, cc = rasterio.transform.rowcol(r.transform, *to.transform(passes.lon.values, passes.lat.values))
    rr, cc = np.asarray(rr), np.asarray(cc)
    inb = (rr >= 0) & (rr < r.height) & (cc >= 0) & (cc < r.width)
    rr, cc = np.where(inb, rr, 0), np.where(inb, cc, 0)
PAT = re.compile(r"(bio\d\d)_(baja|alta)_s(\d+)_env(?:_q([\d.]+))?_")


def at_passes(f):
    """Values of a raster at every pass, indexing the array (much faster than rasterio sample)."""
    with rasterio.open(f) as s:
        v = s.read(1)[rr, cc].astype(float)
        return np.where(inb & (v != s.nodata), v, np.nan)


n = 0
with gzip.open(f"{D}/passes_corriente.csv.gz", "wt", newline="") as gz:      # written per map, so RAM stays flat
    for i, f in enumerate(sorted(glob.glob(f"{T}/bio??_*_s*_env*_corriente.tif"))):
        var, lect, side, q = PAT.match(os.path.basename(f)).groups()
        d = pd.DataFrame({"pass_id": passes.pass_id.values, "var": var, "lectura": lect, "q": float(q or Q0),
                          "lado_km": int(side), "corriente": at_passes(f),
                          "normalizada": at_passes(f.replace("_corriente.tif", "_normalizada.tif"))})
        d.to_csv(gz, index=False, header=i == 0, float_format="%.5g")
        n += 1
print(f"{n * len(passes):,} rows, {len(passes):,} passes, {n} maps")

# latitudinal profile of tile relative conductance, 0.5 degree bands (C8)
TP = re.compile(r"teselas_s(\d+)_bio\d\d(?:_q([\d.]+))?\.csv$")
tt = pd.concat(pd.read_csv(f).assign(lado_km=int(m.group(1)), q=float(m.group(2) or Q0))
               for f in sorted(glob.glob(f"{T}/teselas_s*_bio??*.csv")) if (m := TP.search(f)))
tt = tt[tt.resistencia.str.split("_").str[1] == "env"].assign(var=lambda d: d.resistencia.str[:5])
tt["banda"] = np.floor(tt.lat / 0.5) * 0.5 + 0.25
K = ["var", "lectura", "q", "lado_km", "banda"]
prof = (tt.groupby(K).conductancia_rel
        .agg(n_teselas="size", mediana="median", p25=lambda x: x.quantile(0.25), p75=lambda x: x.quantile(0.75))
        .reset_index().sort_values(K, ascending=[True, True, True, True, False]))
prof.to_csv(f"{D}/perfil_conductancia.csv", index=False, float_format="%.4g")

# maps of normalized current, both readings: every variable for q 0.1 and 200 km, bio01 for the other configurations
for f in sorted(glob.glob(f"{T}/bio??_baja_s*_env*_normalizada.tif")):
    var, _, side, q = PAT.match(os.path.basename(f)).groups()
    if var != "bio01" and (q or side != "200"):
        continue
    fig, ax = plt.subplots(1, 2, figsize=(12, 11), constrained_layout=True)
    for a, lect in zip(ax, ("baja", "alta")):
        with rasterio.open(f.replace("_baja_", f"_{lect}_")) as s:
            z = s.read(1, masked=True)
            ext = [s.bounds.left / 1e3, s.bounds.right / 1e3, s.bounds.bottom / 1e3, s.bounds.top / 1e3]
        im = a.imshow(np.log10(np.clip(z, 0.05, 20)), extent=ext, cmap="RdBu_r", vmin=-1.3, vmax=1.3,
                      interpolation="nearest")
        kx, ky = to.transform(ctl.lon.values, ctl.lat.values)
        a.scatter(np.array(kx) / 1e3, np.array(ky) / 1e3, s=14, c="k", marker="^")
        for x, y, n in zip(kx, ky, ctl.name):
            a.annotate(n, (x / 1e3, y / 1e3), fontsize=6, xytext=(3, 3), textcoords="offset points")
        a.set_title(f"{var}, {'tierra baja, oeste a este' if lect == 'baja' else 'montaña, norte a sur'}, "
                    f"q {float(q or Q0):g}, {side} km")
        a.set_xlabel("km"); a.set_ylabel("km")
    fig.colorbar(im, ax=ax, shrink=0.5, label="log10 corriente normalizada (mayor que 0 concentra, menor que 0 desvía)")
    fig.savefig(f"{D}/mapa_{var}_s{side}" + (f"_q{q}" if q else "") + ".png", dpi=150)
    plt.close(fig)
print("->", D)
