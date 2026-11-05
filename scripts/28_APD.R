# 28_APD.R - APD v1.0 climate layers (plan v2.0 section 8.4). Census layers (APD_divide, APD_cols) come from 06.
# APD_flanks.gpkg: layer flank_points (every col x side x threshold 500/1000/1500 m) and flank_routes (D8 descent,
#                  1000 m, sides that reached the threshold without the nearest-cell fallback).
# APD_climate_metrics.csv: present (WorldClim 2.1, 1970-2000) per col x side x threshold, cols with depth >= 100 m.
# APD_climate_future.csv: CMIP6 ensemble (8 GCMs) per col x SSP x period, gateways at 1000 m, worst side.
# Present and future are separate files because their unit differs (side vs col worst side).
# APD_conservation_gateways.csv, APD_montane_corridor.csv: plan v2.2 conservation layers (29, 30); WDPA only as derived values.
# Rscript scripts/28_APD.R
suppressPackageStartupMessages(library(sf))
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; tb <- file.path(cfg$paths$outputs, "tables")
db <- file.path(cfg$paths$outputs, "database")

fl <- read.csv(file.path(int, "M2_flancos.csv"))
pts <- st_as_sf(fl, coords = c("lon", "lat"), crs = 4326, remove = FALSE, na.fail = FALSE)   # sides with no flank: empty point
rt <- st_read(file.path(int, "M2_flank_routes.gpkg"), quiet = TRUE)
rt$flank_m <- cfg$params$flank_m[[1]]
f <- file.path(db, "APD_flanks.gpkg"); unlink(f)
st_write(pts, f, layer = "flank_points", quiet = TRUE)
st_write(rt, f, layer = "flank_routes", quiet = TRUE, append = TRUE)

m2 <- read.csv(file.path(tb, "M2_metricas_clima_collados.csv"))
keep <- c("col_id", "side_lab", "flank_m", "lon_flank", "lat_flank", "route_km", "fallback", "TO_ann", "TO_mon", "AI_gap", "CLD_gap",
          "maha_pct", "maha_n", grep("^(crest|flank)_", names(m2), value = TRUE))
write.csv(m2[, keep], file.path(db, "APD_climate_metrics.csv"), row.names = FALSE)

ens <- read.csv(file.path(tb, "M3_futuro_ensamble.csv"))
ens <- ens[, setdiff(names(ens), "gcm")]
write.csv(ens, file.path(db, "APD_climate_future.csv"), row.names = FALSE)
message(nrow(pts), " flank points, ", nrow(rt), " routes, ", nrow(m2), " present rows, ", nrow(ens), " future rows")

# conservation attributes (derived only; WDPA itself is not redistributed): 30_proteccion_M3.R
pr <- read.csv(file.path(tb, "M3_prioridad_puertas.csv"))
write.csv(pr[, setdiff(names(pr), c("name", "lon", "lat", "elev_cross", "d_km"))], file.path(db, "APD_conservation_gateways.csv"), row.names = FALSE)
bp <- read.csv(file.path(int, "M3_corredor_bins_proteccion.csv"), check.names = FALSE)
write.csv(bp, file.path(db, "APD_montane_corridor.csv"), row.names = FALSE)
