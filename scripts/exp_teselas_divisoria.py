"""exp_teselas_divisoria.py - EXPERIMENT, not part of the article. Latitudinal profile of the relative conductance of
the tiles that cross the continental divide: a tile crosses it when continental passes lie in its central half along the
direction of the reading, that is, between its two walls (PRD C5 puts the walls in the outer quarters).
The 19 variables are read as two groups, temperature (BIO1 to BIO11) and precipitation (BIO12 to BIO19), each by its
median over the group inside the tile and by the spread from its most to its least limiting variable. They are not
fused into one number: the two groups say different things, and that contrast is the result.
Also Geary's C of the normalized current of the continental passes by distance class (written to geary_passes.csv, not
plotted), to see at what distance neighbouring passes stop resembling each other.
Writes outputs/experimento_barreras_montanas/perfil_divisoria.csv, limitantes.csv, geary_passes.csv and
FigC_perfil_divisoria.
python -P scripts/exp_teselas_divisoria.py principal
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod, Transformer

RUN = sys.argv[1]
D, OUT = f"data/interim/conectividad/{RUN}", "outputs/experimento_barreras_montanas"
sys.path.insert(0, "scripts")
from exp_estilo import GROUP, LAT_TICKS, VSHORT, axis_lat, panel, save, setup       # noqa: E402

setup()
LN, LS = 11.0, -23.5
CFG = {"main": ("s200", ""), "q 0.05": ("s200", "_q0.05"), "q 0.2": ("s200", "_q0.2"), "100 km": ("s100", ""),
       "400 km": ("s400", "")}
VARS = [f"bio{v:02d}" for v in range(1, 20)]
GROUPS = {"temperature": VARS[:11], "precipitation": VARS[11:]}
with rasterio.open(f"{D}/insumos/elev.tif") as s:
    H, W, TR, CRS = s.height, s.width, s.transform, s.crs
to_ll = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
db = pd.read_csv("outputs/database/APD_passes.csv.gz", dtype={"trans_robust": str, "longi_robust": str})
c = db[(db.tropics == 1) & (db["rank"] == "continental")]
x, y = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(c.lon.values, c.lat.values)
rr, cc = rasterio.transform.rowcol(TR, x, y)
div = np.zeros((H, W), int)
div[np.clip(rr, 0, H - 1), np.clip(cc, 0, W - 1)] = 1
I = np.pad(div.cumsum(0).cumsum(1), ((1, 0), (1, 0)))                 # integral image: passes in any rectangle


def count(r0, r1, c0, c1):
    r0, r1, c0, c1 = (np.clip(v, 0, m) for v, m in ((r0, H), (r1, H), (c0, W), (c1, W)))
    return I[r1, c1] - I[r0, c1] - I[r1, c0] + I[r0, c0]


rows, lim = [], []
for cfg, (side, q) in CFG.items():
    S = int(side[1:])
    w = []
    for v in VARS:
        t = pd.read_csv(f"{D}/teselas/teselas_{side}_{v}{q}.csv")
        t = t[t.resistencia.str.endswith("env") | t.resistencia.str.contains("_env_q")]
        w.append(t.set_index(["row", "col", "lectura"]).conductancia_rel.rename(v))
    w = pd.concat(w, axis=1).reset_index()                            # one row per tile, one column per variable
    xx, yy = rasterio.transform.xy(TR, w.row, w.col)
    w["lat"] = to_ll.transform(xx, yy)[1]
    r0, c0 = w.row.values - S // 2, w.col.values - S // 2
    tr = w.lectura.values == "baja"                                   # transverse: walls west and east
    cross = np.where(tr, count(r0, r0 + S, c0 + S // 4, c0 + S - S // 4),
                     count(r0 + S // 4, r0 + S - S // 4, c0, c0 + S)) > 0
    w["band"] = np.floor(w.lat / 0.5) * 0.5 + 0.25
    w["reading"] = np.where(tr, "transverse", "longitudinal")
    w = w[cross & (w.lat <= LN) & (w.lat >= LS)]
    keep = []
    for g, vs in GROUPS.items():
        A = w[vs].values
        w[g], w[f"{g}_lo"], w[f"{g}_hi"] = np.median(A, axis=1), A.min(axis=1), A.max(axis=1)
        keep += [g, f"{g}_lo", f"{g}_hi"]
        if cfg == "main":                                             # which variable limits the group, per tile
            s = pd.DataFrame({"reading": w.reading.values, "group": g,
                              "var": np.array(vs)[A.argmin(axis=1)]})
            lim.append(s.groupby(["reading", "group"])["var"].value_counts(normalize=True).rename("share"))
    g = w.groupby(["reading", "band"])
    rows.append(g[keep].median().join(g.size().rename("tiles")).reset_index().assign(cfg=cfg))
pr = pd.concat(rows)
pr.round(4).to_csv(f"{OUT}/perfil_divisoria.csv", index=False)
lim = pd.concat(lim).reset_index()
lim.round(4).to_csv(f"{OUT}/limitantes.csv", index=False)
m = pr[pr.cfg == "main"]
for rd in ("transverse", "longitudinal"):
    q = m[m.reading == rd].set_index("band").sort_index(ascending=False)
    print(f"\n{rd}\n", q[["temperature", "temperature_lo", "precipitation", "precipitation_lo", "tiles"]].round(2)
          .to_string())

# ---- Geary's C of the normalized current of the continental passes, by distance class (table only) ------------------
u = c.sort_values("depth_m", ascending=False).drop_duplicates("cell_1km")
geo = Geod(ellps="WGS84")
bins = np.arange(0, 205, 5)
gc = {"km": bins[:-1] + 2.5}
for p in ("trans_med", "longi_med"):
    v = u.dropna(subset=[p])
    la, lo, z = v.lat.values, v.lon.values, v[p].values
    num, cnt = np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)
    for i in range(len(v) - 1):
        near = np.abs(la[i + 1:] - la[i]) < 2.0
        d = geo.inv(np.full(near.sum(), lo[i]), np.full(near.sum(), la[i]), lo[i + 1:][near], la[i + 1:][near])[2] / 1e3
        k = np.digitize(d, bins) - 1
        ok = k < len(num)
        np.add.at(num, k[ok], (z[i] - z[i + 1:][near][ok]) ** 2); np.add.at(cnt, k[ok], 1)
    gc[p] = num / cnt / (2 * z.var(ddof=1))                       # 0 identical neighbours, 1 no spatial structure
pd.DataFrame(gc).round(3).to_csv(f"{OUT}/geary_passes.csv", index=False)

# ---- figure -------------------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(183 / 25.4, 105 / 25.4))
for i, rd in enumerate(("transverse", "longitudinal")):
    a = fig.add_axes([0.06 + i * 0.26, 0.11, 0.22, 0.82])
    q = m[m.reading == rd].set_index("band").sort_index()
    for g in GROUPS:
        a.fill_betweenx(q.index, q[f"{g}_lo"], q[f"{g}_hi"], color=GROUP[g], lw=0, alpha=0.22, zorder=1)
        a.plot(q[g], q.index, c=GROUP[g], lw=1.2, label=g, zorder=3)
    a.axvline(1, c="0.5", lw=0.5, ls=":")
    a.set_xlim(0, 1.05)
    axis_lat(a, LS, LN, labels=i == 0)
    a.set_xlabel("Relative conductance of the tiles")
    panel(a, "ab"[i], f"{rd.capitalize()} reading: {'lowland' if i == 0 else 'montane'} climates")
    if i == 0:
        a.set_ylabel("Latitude of the tile centre")
        a.legend(fontsize=5.5, frameon=False, loc="lower left", borderaxespad=0.2, title="Median of the group per\n"
                 "0.5° band; the band spans its\nmost and least limiting variable", title_fontsize=5, alignment="left")
# c: which variable limits each group, as a share of the tiles that cross the divide
a = fig.add_axes([0.70, 0.11, 0.26, 0.82])
yt, ylab, h = [], [], 0.0
for g, vs in GROUPS.items():
    s = lim[lim.group == g].pivot_table(index="var", columns="reading", values="share").fillna(0)
    s = s.loc[s.max(axis=1).sort_values().index]
    s = s[s.max(axis=1) >= 0.04]
    for v in s.index:
        a.barh(h - 0.19, s.loc[v, "transverse"], height=0.34, color=GROUP[g], alpha=0.45, lw=0)
        a.barh(h + 0.19, s.loc[v, "longitudinal"], height=0.34, color=GROUP[g], lw=0)
        yt.append(h); ylab.append(VSHORT[v]); h += 1
    h += 0.7
a.set_yticks(yt); a.set_yticklabels(ylab, fontsize=5.5); a.set_ylim(-0.9, h - 0.9)
a.set_xlim(0, 0.7); a.set_xlabel("Share of the tiles that cross the divide")
a.tick_params(axis="y", length=0)
a.spines[["top", "right", "left"]].set_visible(False)
for lab, al in (("longitudinal", 1.0), ("transverse", 0.45)):
    a.barh(np.nan, 0, color="0.3", alpha=al, label=lab)
a.legend(fontsize=5.5, frameon=False, loc="lower right", title="Reading", title_fontsize=5, alignment="left")
panel(a, "c", "Variable that limits its group")
save(fig, f"{OUT}/FigC_perfil_divisoria")
