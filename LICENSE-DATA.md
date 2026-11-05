# Licences for data, figures and tables

The code (`scripts/`) is under the MIT licence (`LICENSE`). Everything else in this repository follows the rules below, which come from the terms of the inputs. Where a file falls under more than one rule, the most restrictive one applies.

## 1. Census files: CC BY 4.0

Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), with the attributions of section 5:

- `outputs/database/APD_divide.gpkg`, `APD_cols.gpkg`, `APD_cols.csv` and `APD_flanks.gpkg`, except the pass names (section 4)
- `outputs/tables/M1_*.csv`, `outputs/kml/`
- `outputs/supp/FigS1_validation.*`, Tables S1.1 to S1.3 and S2.1
- `data/refs/M1_reference_*.csv`, `data/refs/M1_col_names.csv` (section 4 for names taken from OpenStreetMap)
- `data/interim/M1_*`, `data/interim/M2_flancos.csv`, `data/interim/M2_flank_routes.gpkg`

## 2. Files derived from WorldClim or EarthEnv: CC BY-NC 4.0

WorldClim 2.1 and its CMIP6 projections are free for academic and other non-commercial use, and their redistribution or commercial use needs prior permission. The EarthEnv cloud frequency is under CC BY-NC 4.0. Files that contain values or results computed from them are licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/), and they do not grant any right over the original rasters:

- `outputs/database/APD_climate_metrics.csv`, `APD_climate_future.csv`
- `outputs/tables/M2_*`, `outputs/tables/M3_futuro*.csv`, `M3_dTstar.csv`, `M3_corredor_montano*.csv`, `Table1_gateways.csv`, `Table2_*.csv`
- `outputs/supp/` other than the census files of section 1 and the WDPA files of section 3, `outputs/figs_qc/M2_*`
- `data/interim/M2_maha_vars.txt`, `data/interim/M3_futuro_por_lado.csv`, `data/interim/M3_corredor_bins.csv`, `data/interim/M3_lapse_rate.csv`

Maps made with WorldClim and GADM may be published in research articles, including under open licences, as both providers state. The figures in `outputs/figures/` and `outputs/figs_qc/` are therefore under CC BY 4.0, except Fig. 4 (section 3).

## 3. Files derived from the WDPA: UNEP-WCMC terms

The World Database on Protected Areas may not be redistributed, and its [terms](https://www.protectedplanet.net/en/legal) forbid sub-licensing it and using derived works commercially without written permission from UNEP-WCMC. The WDPA itself is not in this repository. The following files hold only derived attributes (inside or outside a protected area, protected share). They are distributed for non-commercial use, under no open licence of ours, and any further use is subject to those terms:

- `outputs/database/APD_conservation_gateways.csv`, `APD_montane_corridor.csv`
- `outputs/tables/M3_corredor_proteccion.csv`, `M3_corredor_umbral.csv`, `M3_perdida_pais.csv`, `M3_prioridad_puertas.csv`, `Table3_conservation_priority.csv`
- `outputs/figures/Fig4_conservation.*`
- `outputs/supp/tables/S4.1.csv`
- `data/interim/M3_corredor_bins_proteccion.csv`
- `docs/capas/pa/` (map tiles of the web viewer, excluded from the archived release through `.gitattributes`)

Citation required by the provider: UNEP-WCMC and IUCN (2026), Protected Planet: The World Database on Protected Areas (WDPA), November 2026, Cambridge, UK: UNEP-WCMC and IUCN. Available at: www.protectedplanet.net.

## 4. OpenStreetMap: ODbL

`data/refs/osm_mountain_pass_andes.json`, and the pass names in `osm_name` and `name` of `APD_cols` and of the tables that repeat them, come from OpenStreetMap, © OpenStreetMap contributors, and remain under the [Open Database Licence](https://opendatacommons.org/licenses/odbl/).

## 5. Required attributions

- **Copernicus GLO-30:** Copernicus Digital Elevation Model, https://doi.org/10.5270/ESA-c5d3d65, accessed from https://registry.opendata.aws/copernicus-dem. © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018, provided under COPERNICUS by the European Union and ESA; all rights reserved.
- **HydroSHEDS:** The Andean Passes Database incorporates data from the HydroSHEDS version 1 database which is © World Wildlife Fund, Inc. (2006-2022) and has been used herein under license. WWF has not evaluated the data as altered and incorporated within the Andean Passes Database, and therefore gives no warranty regarding its accuracy, completeness, currency or suitability for any particular purpose. Portions of the HydroSHEDS v1 database incorporate data which are the intellectual property rights of © USGS (2006-2008), NASA (2000-2005), ESRI (1992-1998), CIAT (2004-2006), UNEP-WCMC (1993), WWF (2004), Commonwealth of Australia (2007), and Her Royal Majesty and the British Crown and are used under license. The HydroSHEDS v1 database and more information are available at https://www.hydrosheds.org. The scientific citation for the HydroSHEDS version 1 database is: Lehner, B., Verdin, K., Jarvis, A. (2008): New global hydrography derived from spaceborne elevation data. Eos, Transactions, AGU, 89(10): 93-94.
- **HydroBASINS:** Lehner & Grill (2013), under the HydroSHEDS licence terms.
- **GMBA Mountain Inventory v2:** Snethlage et al. (2022), https://doi.org/10.48601/earthenv-t9k2-1407, CC BY 4.0.
- **WorldClim 2.1:** Fick & Hijmans (2017).
- **Global Aridity Index v3.1:** Zomer et al. (2022), https://doi.org/10.6084/m9.figshare.7504448, CC BY 4.0.
- **Cloud frequency:** Wilson & Jetz (2016), EarthEnv, CC BY-NC 4.0.
- **Global Human Modification:** Kennedy et al. (2019), https://doi.org/10.6084/m9.figshare.7283087, CC BY 4.0.
- **Model sensitivity values:** Hausfather et al. (2022), https://doi.org/10.5281/zenodo.6476375, CC BY 4.0.
- **GADM 4.1:** https://gadm.org, map outlines only.

## 6. Inputs not included

The third-party inputs keep their own terms and are not in this repository, either because their terms forbid it (WorldClim, GADM, WDPA) or because of their size. `data/inputs_manifest.csv` gives the path, size and SHA-256 of each file used, and the README gives the script that obtains it.
