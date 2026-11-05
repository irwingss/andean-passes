"""00_verificar_insumos.py - compare data_dir with data/inputs_manifest.csv (path, bytes and SHA-256 of every input used for the article).
python -P scripts/00_verificar_insumos.py              checks the inputs (from the repository root)
python -P scripts/00_verificar_insumos.py --escribir   rewrites the manifest
"""
import csv
import hashlib
import os
import sys

import yaml

raw = yaml.safe_load(open("config.yml", encoding="utf-8"))["paths"]["data_dir"]
MAN = "data/inputs_manifest.csv"
# Unzipped folders are covered by their zip; the rest were never used for this article
SKIP = ("hybas_sa_lev01-12_v1c/", "hyd_sa_dem_15s/", "hyd_sa_dir_15s/", "GMBA_Inventory_v2.0_standard_300/", "ghm/gHM/",
        "pisco/", "HydroRIVERS", "S10W080", "refs/Almazroui", "PROVENANCE.csv")


def sha(f):
    h = hashlib.sha256()
    with open(f, "rb") as s:
        for b in iter(lambda: s.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


if "--escribir" in sys.argv:
    rows = []
    for d, _, fs in os.walk(raw, followlinks=True):
        for x in fs:
            p = os.path.join(d, x)
            r = os.path.relpath(p, raw).replace(os.sep, "/")
            if not r.startswith(SKIP) and "HydroRIVERS" not in r:
                rows.append((r, os.path.getsize(p), sha(p)))
    with open(MAN, "w", newline="", encoding="utf-8") as o:
        w = csv.writer(o, lineterminator="\n")
        w.writerow(["path", "bytes", "sha256"])
        w.writerows(sorted(rows))
    print(len(rows), "files,", round(sum(r[1] for r in rows) / 1e9, 1), "GB")
else:
    bad = 0
    for r in csv.DictReader(open(MAN, encoding="utf-8")):
        p = os.path.join(raw, r["path"])
        if not os.path.exists(p) or str(os.path.getsize(p)) != r["bytes"] or sha(p) != r["sha256"]:
            bad += 1
            print("DIFFERS OR MISSING:", r["path"])
    print("all inputs match" if not bad else f"{bad} inputs differ")
