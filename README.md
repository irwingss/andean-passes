# Andean Passes Database (APD)

Code, derived data, figures and supplementary tables for:

> Saldaña, I.S., Vallejos Bardales, L.M., Ugaz-Cherre, A. et al. Census and climatic connectivity of the cols of the tropical Andes (working title; in preparation).

## About

The census finds every col of the divide network of the Andes in FABDEM V1-2 (Copernicus GLO-30 with forests and buildings removed, 1 arc-second), from 12°N to 56°S: 465 585 cols with a depth of 30 m or more, each with height, depth, width, isolation, GMBA unit, the sea each side drains to and its rank (continental, terminal or tributary). HydroBASINS only labels the sides; it never decides whether a col exists.

The climatic connectivity reads every col of the tropical Andes (11°N to 23.5°S) in two ways, with tiled wall-to-wall circuits on a 1 km equal-area grid and each of the 19 bioclimatic variables of WorldClim 2.1:

- **Transverse reading**, for lowland species: current from west to east between cells warmer than the montane limit T2000(lat).
- **Longitudinal reading**, for montane species: current from north to south between cells colder than that limit.

The resistance of a cell is how far it falls outside the climatic envelope of its guild at its latitude. The normalized current of a col says whether climate concentrates the current there (above 1) or diverts it (below 1); the relative conductance of the tiles says how hard a zone is to cross. Every result is given in five configurations (margins q = 0.1, 0.05 and 0.2; tiles of 200, 100 and 400 km). The analysis is descriptive and makes no statistical test between cols.

Code comments and some file and column names are in Spanish; the paper, figures and tables are in English. Some comments cite the project's plan (PRD) or decision log, which are not part of this archive.

## Contents

| Folder | What it holds |
|---|---|
| `scripts/` | The pipeline; `run_all.sh` gives the order |
| `config.yml` | Paths and fixed parameters |
| `data/refs/` | Control passes, official heights of passes and their sources, and the OpenStreetMap pass nodes of 8 November 2026 |
| `data/interim/` | Study area (`M1_gmba_andes.gpkg`) and sea domains (`M1_domains.gpkg`); the census and connectivity runs (about 17 GB) are not in git |
| `data/inputs_manifest.csv`, `data/inputs_provenance.csv` | Size and SHA-256 of input files; URL, date and licence of each download |
| `outputs/database/` | `APD_cols.csv.gz`: one row per col of the census, with the connectivity of the tropical cols (dictionary in the Supplementary Information) |
| `outputs/figures/` | Figs 1 to 3 (600 dpi PNG and vector PDF) |
| `outputs/tables/` | Tables 1 and 2 |
| `outputs/supp/` | Supplementary Information (Word), Supplementary Tables (Excel), Figs S1 to S6 and Tables S1 to S15 as CSV |
| `docs/` | Web viewer (GitHub Pages), written by `scripts/visor_datos.py`: the normalized current of BIO1 in both readings (Fig. 2) and the cols with a depth of 100 m or more, coloured by rank or by their median over the 19 variables |

Licences differ by file: see `LICENSE` (code) and `LICENSE-DATA.md` (everything else).

## Input data

| Data | Source | How it is obtained | Terms |
|---|---|---|---|
| Elevation | FABDEM V1-2 (Hawker et al. 2022; Neal & Hawker 2023), Earth Engine community catalogue | `39_descarga_fabdem.py`; at 1 km inside `47_conectividad_insumos.py` | CC BY-NC-SA 4.0 |
| Basins | HydroBASINS v1c, level 12 (Lehner & Grill 2013) | `01_download_M1.R` | HydroSHEDS licence |
| Mountain area | GMBA Mountain Inventory v2, standard 300 (Snethlage et al. 2022) | `01_download_M1.R` | CC BY 4.0 |
| Climate, 1970 to 2000 | WorldClim 2.1, 19 bioclimatic variables at 30 arc-seconds (Fick & Hijmans 2017) | `18_descarga_worldclim_M2.R` | Non-commercial, no redistribution |
| Pass nodes | OpenStreetMap, `mountain_pass=yes` and `natural=saddle`, 8 November 2026 | In `data/refs/` | ODbL |
| Official heights | INEI (Peru) and Unidad de Pasos Fronterizos (Chile) | In `data/refs/` | Compiled for this study |

`01_download_M1.R` also fetches layers for the later blocks of the study (human modification, protected areas, map outlines).

## Reproducing the analysis

Requirements: Python 3.13 (`requirements.txt`) and R 4.6 with sf, terra and yaml (`Rscript scripts/00_setup.R`). Run every command from the repository root. Inputs go to `data/raw/`, or wherever `data_dir` points in `config.yml`; FABDEM is read from Earth Engine, which needs a registered Earth Engine project.

```bash
Rscript scripts/00_setup.R
Rscript scripts/01_download_M1.R
Rscript scripts/18_descarga_worldclim_M2.R
bash scripts/run_all.sh
```

| Step | Scripts | What it does |
|---|---|---|
| Census | `02`, `39`, `40`, `41`, `45` | Study area and sea domains, FABDEM, divide network and cols per block of 1°, labels and rank, blocks joined and isolation |
| Census in the cloud | `44_censo_nube.sh`, `44_vm_censo.sh` | The same census on one Google Cloud machine, 16 blocks at a time (31 min, about 1.2 USD) |
| Validation | `46`, `42` | OpenStreetMap nodes, identity against random divide cells, official heights, control passes; review sheet of a test block |
| Connectivity | `47`, `50`, `49` | Grid, variables, guilds and resistances; tiled circuits, one variable per call; values at every col |
| Products | `51`, `52`, `53` | Database and tables, figures, Supplementary Information |

Without the census and connectivity runs, `51`, `52` and `53` cannot rerun; their outputs are in `outputs/`. The connectivity of the main configuration took 74 min on a 16-core desktop, and the four sensitivity configurations about 6 h more.

## Citation

See `CITATION.cff`. The DOI is assigned by Zenodo at release.
