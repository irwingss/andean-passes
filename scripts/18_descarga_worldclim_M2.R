# 18_descarga_worldclim_M2.R - WorldClim 2.1 (1970-2000, 30") for the M2 sensitivity, cropped to the Andes box.
# Three 30-deg tiles (lon -90..-60; lat 30..0, 0..-30, -30..-60) of bio, tmax and tmin via geodata.
# Writes WC_bio*.tif and WC_tasmax/tasmin_MM.tif to data_dir/clima_wc. Rscript scripts/18_descarga_worldclim_M2.R
suppressPackageStartupMessages(library(terra))
cfg <- yaml::read_yaml("config.yml"); raw <- cfg$paths$data_dir
dir <- file.path(raw, "clima_wc"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 7200)
box <- ext(-82, -60, -56, 13)

get <- function(v) {
  t <- lapply(c(5, -15, -45), function(lat) geodata::worldclim_tile(v, -75, lat, path = dir))
  crop(merge(sprc(t)), box)
}
bio <- c(1, 5, 6, 12, 14, 15)
if (!all(file.exists(file.path(dir, paste0("WC_bio", bio, ".tif"))))) {
  b <- get("bio")
  for (k in bio) writeRaster(b[[k]], file.path(dir, paste0("WC_bio", k, ".tif")), overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
}
for (v in c("tmax", "tmin")) {
  f <- file.path(dir, sprintf("WC_tas%s_%02d.tif", sub("t", "", v), 1:12))
  if (all(file.exists(f))) next
  r <- get(v); for (k in 1:12) writeRaster(r[[k]], f[k], overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
}

prov_f <- file.path(raw, "PROVENANCE.csv"); prov <- read.csv(prov_f)
k <- "clima_wc/ (WorldClim 2.1 30s tiles 'bio', 'tmax', 'tmin' x lat 5/-15/-45, lon -75; cropped)"
if (!k %in% prov$file) {
  tiles <- list.files(file.path(dir, "climate", "tiles_2.1"), "\\.tif$", full.names = TRUE)
  prov[nrow(prov) + 1, ] <- list(k, "https://geodata.ucdavis.edu/climate/worldclim/2_1/tiles/tile/",
                                 format(Sys.Date()), digest::digest(paste(sapply(tiles, digest::digest, file = TRUE, algo = "sha256")), algo = "sha256"),
                                 "WorldClim 2.1 (free for academic use, attribution: Fick & Hijmans 2017)")
  write.csv(prov, prov_f, row.names = FALSE)
}
message("done")
