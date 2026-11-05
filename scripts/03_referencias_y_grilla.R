# 03_referencias_y_grilla.R - two small inputs of the census: the validation passes in decimal degrees and the 15 s grid
# of the domains (W = 1, E or I = 2), which 08_marcadores.R and 14_flancos_M2.R use as their grid template.
# Rscript scripts/03_referencias_y_grilla.R
suppressPackageStartupMessages(library(terra))
cfg <- yaml::read_yaml("config.yml"); raw <- cfg$paths$data_dir; int <- cfg$paths$interim

refs <- read.csv("data/refs/M1_reference_passes.csv", encoding = "UTF-8")
dms <- function(s) { x <- as.numeric(regmatches(s, gregexpr("[0-9.]+", s))[[1]]); -(x[1] + x[2] / 60 + x[3] / 3600) }
fill <- is.na(refs$lat) & !is.na(refs$lat_dms) & nzchar(refs$lat_dms)
refs$lat[fill] <- sapply(refs$lat_dms[fill], dms); refs$lon[fill] <- sapply(refs$lon_dms[fill], dms)
write.csv(refs, file.path(int, "M1_reference_passes_dd.csv"), row.names = FALSE)

dem <- crop(rast(file.path(raw, "hyd_sa_dem_15s/hyd_sa_dem_15s.tif")), ext(-82, -60, -56.5, 13))
dom <- vect(file.path(int, "M1_domains.gpkg")); dom$code <- ifelse(dom$dom == "W", 1, 2)
rasterize(dom, dem, field = "code", filename = file.path(int, "M1_domains_15s.tif"), overwrite = TRUE)
