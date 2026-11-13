"""exp_top_continentales.py - EXPERIMENT, not part of the article: continental passes of the tropical Andes ranked by
the median normalized current over the 19 variables (main configuration), one table per reading. WE lists the
bridges for lowland climates: the highest transverse current, passes that concentrate it (above 1). NS lists the
barriers for montane climates: the lowest longitudinal current, passes that divert it (below 1). Neighbouring passes of
the same crossing are grouped into places: going down the ranking, a pass opens a new place when it stays on that side of 1
in all five configurations and lies 25 km or more from every earlier place.
Writes outputs/experimento_top_continentales/continental_passes_{WE,NS}.csv (every pass) and places_{WE,NS}.csv (top 20).
python -P scripts/exp_top_continentales.py
"""
import os

import numpy as np
import pandas as pd

OUT = "outputs/experimento_top_continentales"
SEP_KM, N = 25.0, 20
os.makedirs(OUT, exist_ok=True)
db = pd.read_csv("outputs/database/APD_passes.csv.gz", dtype={"trans_robust": str, "longi_robust": str})
c = db[(db.tropics == 1) & (db["rank"] == "continental")]
c = c.sort_values("depth_m", ascending=False).drop_duplicates("cell_1km")          # one pass per 1 km cell, the deepest


def km(lat, lon, la, lo):
    """Great-circle distance (km) from one point to arrays of points."""
    p, q = np.radians(lat), np.radians(la)
    h = np.sin((q - p) / 2) ** 2 + np.cos(p) * np.cos(q) * np.sin(np.radians(lo - lon) / 2) ** 2
    return 12742 * np.arcsin(np.sqrt(h))


for p, tag, side in (("trans", "WE", "above 1"), ("longi", "NS", "below 1")):
    cfg = [f"{p}_med"] + [f"{p}_med_{s}" for s in ("q005", "q02", "100km", "400km")]
    t = c.dropna(subset=[f"{p}_med"]).sort_values(f"{p}_med", ascending=side == "below 1").reset_index(drop=True)
    t["cfg_lo"], t["cfg_hi"] = t[cfg].min(axis=1, skipna=False), t[cfg].max(axis=1, skipna=False)
    t.insert(0, "rank_pass", t.index + 1)
    place, pla, plo = np.full(len(t), np.nan), [], []
    for i, r in t[t[f"{p}_robust"] == side].iterrows():
        if not pla or km(r.lat, r.lon, np.array(pla), np.array(plo)).min() >= SEP_KM:
            pla.append(r.lat); plo.append(r.lon); place[i] = len(pla)
    t["place"] = place
    cols = ["rank_pass", "place", "apd_id", "lat", "lon", "z_m", "depth_m", "gmba", f"{p}_med", f"{p}_min", f"{p}_max",
            "cfg_lo", "cfg_hi", f"{p}_robust"]
    t = t[cols].round({"lat": 4, "lon": 4, "z_m": 0, "depth_m": 0, f"{p}_med": 3, f"{p}_min": 3, f"{p}_max": 3,
                       "cfg_lo": 3, "cfg_hi": 3})
    t.to_csv(f"{OUT}/continental_passes_{tag}.csv", index=False)
    top = t[t.place <= N].copy()
    top["place"] = top.place.astype(int)
    top.to_csv(f"{OUT}/places_{tag}.csv", index=False)
    print(tag, len(t), "passes;", int(np.nanmax(place)), "places", side, "in all five configurations")
    print(top.drop(columns=["rank_pass"]).to_string(index=False))
