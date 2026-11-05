# APD v1.0 · Andean Passes Database

Census of the cols (mountain passes) on the W–E continental divide of the Andes (11°N–55°S) and their climatic permeability under present (1970–2000) and future (CMIP6, 2041–2100) climates. Data product of the article *Climatic Permeability of Andean Mountain passes: a continental census under present and future climates* (submitted).

## Files

| File | Content |
|---|---|
| `APD_divide.gpkg` | W–E continental divide (Copernicus GLO-30) |
| `APD_cols.gpkg` / `APD_cols.csv` | 4579 cols with census attributes (height, depth, isolation, width, gateway flags, names) |
| `APD_flanks.gpkg` | Lowland flank points (500 / 1000 / 1500 m) and D8 descent routes (1000 m) |
| `APD_climate_metrics.csv` | Present crest–flank climate metrics per col × side × threshold |
| `APD_climate_future.csv` | Future metrics per gateway × SSP × period (8-GCM ensemble: mean, range, sign agreement) |
| `APD_conservation_gateways.csv` | Tropical gateways: protection, human pressure, E–W opening, montane bottleneck and priority class |
| `APD_montane_corridor.csv` | Tropical divide (0.5 km bins): montane gap today and under each SSP × period, robust loss, protection, gHM |
| `APD_data_dictionary.csv` | Every column, its meaning and units |

Join key: `col_id`. Coordinates: WGS84 (EPSG:4326).

## Main definitions

- **Crossing height:** the lowest W–E route through each 0.5 km divide segment. A **col** is a local minimum of that profile.
- **Depth:** inverted prominence. **Isolation:** distance along the divide to the nearest lower point.
- **Gateway:** a col with depth ≥ 100 m that is the lowest within ±100 km of divide.
- **Flank:** where the steepest-descent trace from the col reaches the threshold elevation (default 1000 m), on each side.
- **Col permeability:** the worst side (lower thermal overlap, larger absolute moisture gap).

## Sources

Copernicus GLO-30 DEM; HydroSHEDS (Lehner et al. 2008) and HydroBASINS (Lehner & Grill 2013); GMBA Mountain Inventory v2 (Snethlage et al. 2022); WorldClim 2.1 and its CMIP6 downscaled futures (Fick & Hijmans 2017); Global Aridity Index v3.1 (Zomer et al. 2022); EarthEnv MODIS cloud frequency (Wilson & Jetz 2016); OpenStreetMap pass names (ODbL); Global Human Modification (Kennedy et al. 2019, CC BY 4.0); WDPA, November 2026 (UNEP-WCMC & IUCN, Protected Planet; only derived attributes are included, the WDPA itself is not redistributed). The path, size and SHA-256 of every input are in `data/inputs_manifest.csv` of the repository.

## Code

Built by the scripts in `scripts/`: `06_validacion_huancabamba.R` writes the census layers and `28_APD.R` the flank, climate and conservation layers.

## Licence

Licences differ by file; `LICENSE-DATA.md` in the repository root gives the rule for each. In short: the census layers are CC BY 4.0; climate files derived from WorldClim and EarthEnv are CC BY-NC 4.0; files derived from the WDPA are for non-commercial use under the terms of UNEP-WCMC; pass names taken from OpenStreetMap remain under ODbL.

## Citation

See `CITATION.cff` in the repository root. The DOI is assigned by Zenodo at release.
