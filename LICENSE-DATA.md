# Licences for data, figures and tables

The code (`scripts/`) is under the MIT licence (`LICENSE`). Everything else in this repository follows the rules below, which come from the terms of the inputs. Where a file falls under more than one rule, the most restrictive one applies.

## 1. Files derived from FABDEM: CC BY-NC-SA 4.0

FABDEM V1-2 is distributed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). Every height, depth, width and isolation of the census comes from it, so the census and everything built on it are licensed under CC BY-NC-SA 4.0, with the attributions of section 4:

- `outputs/database/APD_cols.csv.gz`
- `outputs/tables/`, `outputs/supp/` (tables, figures and the Supplementary Information)
- `data/interim/` and `data/refs/pasos_control.csv`, `data/refs/M1_reference_*.csv`

## 2. Files derived from WorldClim

WorldClim 2.1 is free for academic and other non-commercial use, and its redistribution or commercial use needs prior permission. The connectivity values of the database and the tables and figures computed from them contain results derived from WorldClim; they are distributed for non-commercial use and grant no right over the original rasters. Maps made with WorldClim may be published in research articles, as the provider states.

## 3. OpenStreetMap: ODbL

`data/refs/osm_pasos_20261108.json` and the validation tables that repeat its nodes (Tables S3 and S6) come from OpenStreetMap, © OpenStreetMap contributors, and remain under the [Open Database Licence](https://opendatacommons.org/licenses/odbl/).

## 4. Required attributions

- **FABDEM V1-2:** Hawker, L. et al. (2022), Environmental Research Letters 17, 024016; Neal, J. & Hawker, L. (2023), University of Bristol, https://doi.org/10.5523/bris.s5hqmjcdj8yo2ibzi9b4ew3sn. FABDEM is derived from the Copernicus DEM GLO-30, © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018, provided under COPERNICUS by the European Union and ESA.
- **HydroSHEDS and HydroBASINS:** The Andean Passes Database incorporates data from the HydroSHEDS version 1 database which is © World Wildlife Fund, Inc. (2006-2022) and has been used herein under license. WWF has not evaluated the data as altered and incorporated within the Andean Passes Database, and therefore gives no warranty regarding its accuracy, completeness, currency or suitability for any particular purpose. Portions of the HydroSHEDS v1 database incorporate data which are the intellectual property rights of © USGS (2006-2008), NASA (2000-2005), ESRI (1992-1998), CIAT (2004-2006), UNEP-WCMC (1993), WWF (2004), Commonwealth of Australia (2007), and Her Royal Majesty and the British Crown and are used under license. The HydroSHEDS v1 database and more information are available at https://www.hydrosheds.org. HydroBASINS: Lehner & Grill (2013).
- **GMBA Mountain Inventory v2:** Snethlage et al. (2022), https://doi.org/10.48601/earthenv-t9k2-1407, CC BY 4.0.
- **WorldClim 2.1:** Fick & Hijmans (2017).

## 5. Inputs not included

The third-party inputs keep their own terms and are not in this repository, either because their terms forbid it (WorldClim) or because of their size. `data/inputs_manifest.csv` gives the path, size and SHA-256 of each file used, and the README gives the script that obtains it.
