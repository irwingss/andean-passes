# Andean Passes Viewer

Web viewer of the Andean Passes Database (APD v1.0): https://irwingss.github.io/andean-passes/

It shows the 4579 cols of the W–E continental divide of the Andes (11°N to 55°S), the drainage domains, the montane gaps at 2000, 2500 and 3000 m, the descent routes of the gateways, the montane corridor under SSP3-7.0, protected areas and human modification, in a Leaflet viewer. The article (*Gateways through the Andes*) is in preparation; the results here are provisional until it is published.

| Path | Content |
|---|---|
| `docs/` | The viewer (GitHub Pages: branch main, folder /docs) and its layers in `docs/capas/` |
| `visor_datos.py` | Writes `docs/capas/` from `data/interim` and `outputs/` of the project |

```bash
python visor_datos.py tiles
```

Python with numpy, pandas, pyogrio, rasterio, shapely and Pillow. Without `tiles` it keeps the raster tiles and rewrites the rest.

Protected areas come from the WDPA, November 2026 (UNEP-WCMC & IUCN, Protected Planet) and appear only as a derived raster map; the WDPA itself is not redistributed. Layers are CC BY 4.0, and pass names taken from OpenStreetMap remain under ODbL.
