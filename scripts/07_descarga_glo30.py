"""07_descarga_glo30.py - Copernicus GLO-30 tiles of the divide corridor (list from 08_marcadores.R) into data_dir/glo30.
Tiles that do not exist on AWS (ocean only) go to M1_glo30_absent.txt. Idempotent. The tiles take about 8 GB.
python -P scripts/07_descarga_glo30.py   (from the repository root, after 08_marcadores.R)
"""
import csv
import os
import shutil
import urllib.error
import urllib.request

import yaml

cfg = yaml.safe_load(open("config.yml", encoding="utf-8"))
INT, dst = cfg["paths"]["interim"], os.path.join(cfg["paths"]["data_dir"], "glo30")
os.makedirs(dst, exist_ok=True)
tiles = [r["name"] for r in csv.DictReader(open(os.path.join(INT, "M1_glo30_tiles.csv")))]
open(os.path.join(INT, "M1_glo30_tiles.txt"), "w").write("\n".join(tiles) + "\n")
absent = []
for t in tiles:
    f = os.path.join(dst, t + ".tif")
    if os.path.exists(f):
        continue
    try:
        with urllib.request.urlopen(f"{cfg['remote']['glo30_bucket']}/{t}/{t}.tif") as r, open(f + ".part", "wb") as o:
            shutil.copyfileobj(r, o)
        os.replace(f + ".part", f)
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
        absent.append(t)
if absent:
    open(os.path.join(INT, "M1_glo30_absent.txt"), "w").write("\n".join(absent) + "\n")
print(len(tiles), "tiles,", len(absent), "absent")
