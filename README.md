# Andean Passes Database (APD v1.0)

Code, derived data, figures and supplementary tables for:

> Saldaña, I.S., Vallejos Bardales, L.M., Ugaz-Cherre, A. et al. Climatic Permeability of Andean Mountain passes: a continental census under present and future climates (submitted).

## Interactive viewer

**https://irwingss.github.io/andean-passes/**

The 4579 cols of the Pacific–Atlantic divide of the Andes (11°N to 55°S), the drainage domains, the montane gaps at 1500, 2000 and 2500 m, the descent routes of the gateways, the montane corridor under SSP3-7.0, protected areas and human modification, in a Leaflet viewer.

## About

The pipeline extracts the drainage divide of the whole cordillera from the 30 m Copernicus DEM by marker-controlled priority flooding, finds every col along it and measures its crossing height, depth, isolation and width. For each gateway it then measures how far the climate of the crest overlaps the climate of the lowland flank on each side, tests Janzen's prediction with a Gaussian GAM, explores a moisture axis in the tropical Andes, and projects the thermal and precipitation contrasts with eight CMIP6 models under three scenarios. Finally it maps where the montane corridor along the divide would break, and how much of those places is protected or modified by human use.

Code comments and some file and column names are in Spanish; the paper, figures and tables are in English. Some comments cite the project's decision log (BITACORA) or plan, which are not part of this archive. A few columns and rows named `huanca` hold extra cols for a companion study; they do not enter any result of this article.

## Contents

| Folder | What it holds |
|---|---|
| `scripts/` | The whole pipeline, numbered in the order the analysis was built. `run_all.sh` gives the order in which to run it |
| `config.yml` | Paths and every fixed parameter |
| `data/refs/` | Reference passes and their sources, gateway names, and the OpenStreetMap passes used for names |
| `data/interim/` | Intermediate results needed to rerun the statistics, tables and database without the inputs; includes the flood blocks of the divide |
| `data/inputs_manifest.csv`, `data/inputs_provenance.csv` | Path, size and SHA-256 of every input file used; URL, date and licence of each download |
| `outputs/database/` | The Andean Passes Database (`README.md` and data dictionary inside) |
| `outputs/figures/` | Figs. 1–5 (600 dpi PNG and vector PDF) and their captions |
| `outputs/tables/` | Tables 1–3 and the tables behind every number of the text |
| `outputs/supp/` | Supporting Information (Word and Excel), its figures and tables, and the S4 atlas |
| `outputs/figs_qc/`, `outputs/kml/` | Quality-control figures, including the map of drainage domains (`M1_dominios.png`, Fig. S1.1), and the cols as KML |
| `docs/` | Web viewer (GitHub Pages: branch main, folder `/docs`), written by `scripts/visor_datos.py` |

Licences differ by file: see `LICENSE` (code) and `LICENSE-DATA.md` (everything else).

## Input data

| Data | Source | How it is obtained | Terms |
|---|---|---|---|
| Elevation | Copernicus GLO-30, public release | `07_descarga_glo30.py` | Free, with attribution |
| Basins | HydroBASINS v1c, level 12 (Lehner & Grill 2013) | `01_download_M1.R` | HydroSHEDS licence |
| Elevation and flow direction, 15″ | HydroSHEDS v1 (Lehner et al. 2008) | `01_download_M1.R` | HydroSHEDS licence |
| Mountain limits | GMBA Mountain Inventory v2, standard 300 (Snethlage et al. 2022) | `01_download_M1.R` | CC BY 4.0 |
| Climate, 1970–2000 | WorldClim 2.1, 30″ (Fick & Hijmans 2017) | `18_descarga_worldclim_M2.R`, `19_descarga_M3.R` | Non-commercial, no redistribution |
| Future climate | WorldClim 2.1 CMIP6, 10′, 8 models, 3 scenarios, 3 periods | `19_descarga_M3.R` | Non-commercial, no redistribution |
| Aridity | Global Aridity Index v3.1 (Zomer et al. 2022) | `13_descarga_clima_M2.R` | CC BY 4.0 |
| Cloud frequency | EarthEnv MODIS (Wilson & Jetz 2016) | `13_descarga_clima_M2.R` | CC BY-NC 4.0 |
| Human modification | Global Human Modification (Kennedy et al. 2019) | `01_download_M1.R` | CC BY 4.0 |
| Protected areas | WDPA, November 2026 (UNEP-WCMC & IUCN 2026) | Manual, from Protected Planet (`01_download_M1.R` prints the steps) | No redistribution; non-commercial |
| Map outlines | GADM 4.1, level 1 | `01_download_M1.R` | Non-commercial, no redistribution |
| Model sensitivities | Hausfather et al. (2022) | `01_download_M1.R` | CC BY 4.0 |
| Reference passes | INEI, Unidad de Pasos Fronterizos de Chile, literature (Tables S1.2, S1.3) | In `data/refs/` | Compiled for this study |
| Pass names | OpenStreetMap, Overpass export of 2 November 2026 | In `data/refs/` | ODbL |

The WDPA changes every month, and the November 2026 release may no longer be offered; a later release gives slightly different protection results. WorldClim also gives no version archive, so `scripts/00_verificar_insumos.py` tells you whether your inputs are the files used for the article.

## Reproducing the analysis

Requirements: Python 3.13 (`requirements.txt`), R 4.6 (`Rscript scripts/00_setup.R` installs the packages; versions below) and Inkscape 1.4 for Fig. 5. Run the commands from the repository root. Inputs take about 25 GB in `data/raw/`, or wherever `data_dir` points in `config.yml`.

**Without the inputs.** With the repository alone these steps rerun: cols and validation, statistics of Janzen's prediction and of the moisture axis, tables and database.

```bash
Rscript scripts/04_collados.R --force
Rscript scripts/06_validacion_huancabamba.R
Rscript scripts/17_janzen_M2.R --force
Rscript scripts/27_tablas_manuscrito.R
Rscript scripts/28_APD.R
```

**From the inputs.**

```bash
# 0. Inputs: downloads, then check against the files used for the article
Rscript scripts/00_setup.R
Rscript scripts/01_download_M1.R
Rscript scripts/13_descarga_clima_M2.R
Rscript scripts/18_descarga_worldclim_M2.R
Rscript scripts/19_descarga_M3.R
# (the GLO-30 tiles come later, inside run_all.sh, once 08 has listed them)
# 1. Everything else, in order
rm -rf data/interim/M1_glo30_divide   # only to redo the flood; otherwise its blocks are reused
bash scripts/run_all.sh
python -P scripts/00_verificar_insumos.py
```

| Step | Scripts | What it does |
|---|---|---|
| Census | `02`, `03`, `08`, `07`, `09`, `09b`, `10`, `04`, `06` | Drainage domains, seeds, GLO-30 tiles, flood in 1° blocks and in one window over the Altiplano–Puna, profile, cols, validation and sensitivity, census layers of the APD |
| Quality control | `11`, `12` | `outputs/figs_qc/`, KML |
| Present climate | `14`, `15`, `17`, `31` | Flanks, climate metrics, test of Janzen's prediction and moisture axis, diagnosis of the monthly index |
| Future climate | `20`, `29`, `30` | Future permeability, montane corridor, protection and human modification |
| Products | `26`, `27`, `28`, `32`, `33`, `34`, `visor_datos.py` | Figs. 1–4, tables, database, Fig. 5, supplement, viewer layers |

The whole chain took about one hour on a desktop computer, most of it in the flood. `scripts/visor_datos.py` keeps the map tiles unless called with `tiles`. Fig. 5 caches its 3D renders in `data/interim/fig5_renders/`.

A full rerun from a clean copy of this repository reproduced every table of `outputs/` byte for byte. ggplot2 may swap the order of two legends between runs, which changes Fig. 4 and the atlas pages in layout only.

Versions used: R 4.6.1 with ape 5.8-1, data.table 1.18.6.1, digest 0.6.39, geodata 0.6-9, geosphere 1.6-8, ggplot2 4.0.3, ggrepel 0.9.8, jsonlite 2.0.0, maps 3.4.4, mgcv 1.9-4, openxlsx 4.2.9, patchwork 1.3.2, png 0.1-9, ragg 1.5.2, rayshader 0.41.2, rgl 1.3.36, scales 1.4.0, sf 1.1-3, terra 1.9-50, yaml 2.3.12. Python 3.13.15. The scripts contain no random step that enters a reported result.

## Citation

See `CITATION.cff`. The DOI is assigned by Zenodo at release.
