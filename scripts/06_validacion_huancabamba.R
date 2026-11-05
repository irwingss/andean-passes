# 06_validacion_huancabamba.R - M1.6, M1.8: validation against documented passes, key results (rank and isolation),
# depth / window sensitivity, APD v1.0 census layers. Huancabamba system (old M1.7): L2 scripts/06b_sistema_huancabamba.R.
# Inputs: GLO-30 profile (10_perfil_glo30.R) and its cols (04_collados.R).
# Rscript scripts/06_validacion_huancabamba.R
suppressPackageStartupMessages({ library(sf); library(terra) })
sf_use_s2(FALSE)
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; P <- cfg$params
tab <- "outputs/tables"; db <- "outputs/database"
dir.create(tab, recursive = TRUE, showWarnings = FALSE); dir.create(db, recursive = TRUE, showWarnings = FALSE)
km <- function(lon, lat, x, y) geosphere::distGeo(cbind(lon, lat), cbind(x, y)) / 1000

prof <- read.csv(file.path(int, "M1_profile.csv"))
cols <- read.csv(file.path(int, "M1_cols_raw.csv"))
ref <- read.csv(file.path(int, "M1_reference_passes_dd.csv"), encoding = "UTF-8")
cols <- cols[order(cols$d_km), ]
# Terminal cols: within 50 km of either end of the divide chain (cordillera tapering out, not a crossing)
cols$terminal <- cols$d_km - min(prof$d_km) < 50 | max(prof$d_km) - cols$d_km < 50
# Fuegian channels (Beagle, Magellan): sea-level "cols" where the W/E split is ambiguous; not crossings
cols$terminal <- cols$terminal | (cols$lat < -52 & cols$elev_cross < 50)
# Lara: NE end of the Merida Andes (Trujillo-Lara depression and beyond); flagged, not removed
cols$lara <- cols$lat > 9.3 & cols$lon > -70.4
# Names: nearest named OSM mountain_pass node within 3 km
osm <- jsonlite::fromJSON("data/refs/osm_mountain_pass_andes.json")$elements
osm <- osm[!is.na(osm$tags$name) & !grepl("^[0-9 .,]+m?$", osm$tags$name), ]   # drop names that are only a height ("2300 m")
cols$osm_name <- sapply(seq_len(nrow(cols)), function(j) {
  d <- km(osm$lon, osm$lat, cols$lon[j], cols$lat[j]); if (min(d) <= 3) osm$tags$name[which.min(d)] else NA })
# name: each named point goes to its single nearest col with depth >= 50 m (OSM <= 3 km; then the proper
# names in data/refs/M1_col_names.csv, <= 1.5 km, which win over OSM)
nmr <- read.csv("data/refs/M1_col_names.csv", encoding = "UTF-8")
cols$name <- NA_character_
put <- function(lon, lat, nm, maxkm) for (k in seq_along(nm)) {
  d <- km(cols$lon, cols$lat, lon[k], lat[k]); d[cols$depth < 50] <- Inf
  if (min(d) <= maxkm) cols$name[which.min(d)] <<- nm[k] else if (maxkm < 3) warning("no col for name ", nm[k])
}
put(osm$lon, osm$lat, osm$tags$name, 3)
put(nmr$lon, nmr$lat, nmr$name, 1.5)

# --- M1.6 Validation: nearest census col (any depth) within 5 km of each reference ---
# Compared value: topographic GLO-30 crest at the col (elev); elev_cross adds the I|E barrier on W|I cols.
v <- do.call(rbind, lapply(seq_len(nrow(ref)), function(k) {
  d <- km(cols$lon, cols$lat, ref$lon[k], ref$lat[k]); j <- which.min(d)
  dd <- min(km(prof$lon, prof$lat, ref$lon[k], ref$lat[k]))   # is the pass on the W|E divide at all?
  data.frame(ref[k, c("name", "country", "lat", "lon", "elev_ref_m", "elev_alt_m", "source")],
             dist_divide_km = dd, dist_col_km = d[j], col_bin = cols$bin[j], col_depth = cols$depth[j],
             elev_glo30 = cols$elev[j])
}))
v$matched <- v$dist_col_km <= 5
v$err_glo30 <- v$elev_glo30 - v$elev_ref_m
v$err_glo30_alt <- v$elev_glo30 - ifelse(is.na(v$elev_alt_m), v$elev_ref_m, v$elev_alt_m)
write.csv(v, file.path(tab, "M1_validacion.csv"), row.names = FALSE)
m <- v[v$matched, ]
val <- data.frame(n_refs = nrow(v), n_matched = nrow(m),
                  median_abs_err_glo30 = median(abs(m$err_glo30)), median_bias_glo30 = median(m$err_glo30),
                  mean_bias_glo30 = mean(m$err_glo30), median_abs_err_glo30_altsrc = median(abs(m$err_glo30_alt)))
val$criterion_met <- val$median_abs_err_glo30 <= 50 & abs(val$median_bias_glo30) <= 30
write.csv(val, file.path(tab, "M1_validacion_resumen.csv"), row.names = FALSE)
print(val)

# --- Key result: Porculla rank (elevation) and isolation among tropical cols ---
trop <- cols[cols$tropical & cols$depth >= P$depth_m[1] & !cols$terminal, ]
trop <- trop[order(trop$elev_cross), ]
trop$rank_trop <- seq_len(nrow(trop))
trop$rank_trop_sin_lara <- cumsum(!trop$lara) * !trop$lara
trop$rank_iso_trop <- rank(-trop$iso_km, ties.method = "min")   # H1 prediction: Porculla ranks 1
kc <- c("rank_trop", "rank_trop_sin_lara", "rank_iso_trop", "lat", "lon", "name", "osm_name", "elev", "elev_cross",
        "depth", "iso_km", "iso_N_km", "iso_S_km", "width_km", "gateway", "borders_I", "lara")
write.csv(head(trop[, kc], 40), file.path(tab, "M1_collados_tropicales_mas_bajos.csv"), row.names = FALSE)
write.csv(trop[order(trop$rank_iso_trop)[1:15], kc], file.path(tab, "M1_collados_tropicales_mas_aislados.csv"), row.names = FALSE)
pc <- trop[which.min(km(trop$lon, trop$lat, -79.505, -5.839)), ]
cat("n tropical cols:", nrow(trop), "\nPorculla:", pc$elev_cross, "m  rank:", pc$rank_trop, "/ sin Lara:", pc$rank_trop_sin_lara,
    " iso rank:", pc$rank_iso_trop, " iso N/S km:", pc$iso_N_km, pc$iso_S_km, " dist to ref km:",
    round(km(pc$lon, pc$lat, -79.505, -5.839), 2), "\n")

# --- M1.8 Sensitivity: depth thresholds x gateway windows ---
e <- prof$elev_cross; d <- prof$d_km
sens <- expand.grid(depth = P$depth_m, win = P$gateway_window_km)
sens <- cbind(sens, t(apply(sens, 1, function(s) {
  g <- cols[cols$depth >= s["depth"], ]
  isg <- sapply(seq_len(nrow(g)), function(j) g$elev_cross[j] <= min(e[abs(d - g$d_km[j]) <= s["win"]], na.rm = TRUE))
  gt <- g[isg & g$tropical & !g$terminal, ]; gt <- gt[order(gt$elev_cross), ]
  c(n_cols = nrow(g), n_gateways = sum(isg), n_gateways_trop = sum(isg & g$tropical),
    lowest_trop_lat = gt$lat[1], lowest_trop_elev = gt$elev_cross[1])
})))
write.csv(sens, file.path(tab, "M1_sensibilidad_umbral_ventana.csv"), row.names = FALSE)

# --- APD v1.0 census layers ---
invisible(file.copy(file.path(int, "M1_divide_glo30.gpkg"), file.path(db, "APD_divide.gpkg"), overwrite = TRUE))
cs <- cols; cs$col_id <- sprintf("APD%05d", seq_len(nrow(cs)))
st_write(st_as_sf(cs, coords = c("lon", "lat"), crs = 4326, remove = FALSE),
         file.path(db, "APD_cols.gpkg"), delete_dsn = TRUE, quiet = TRUE)
write.csv(cs, file.path(db, "APD_cols.csv"), row.names = FALSE)

# QC figures (profiles, elevation vs isolation): 11_figuras_M1.R
