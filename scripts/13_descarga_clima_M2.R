# 13_descarga_clima_M2.R - M2 hydric and cloud inputs cropped to the Andes box (-82..-60, -56..13), 30".
# Global Aridity Index v3.1 (Zomer et al. 2022; built on WorldClim 2.1) and EarthEnv MODIS cloud frequency
# (Wilson & Jetz 2016). Temperature and bio layers: 18_descarga_worldclim_M2.R.
# Idempotent: skips layers already present. Rscript scripts/13_descarga_clima_M2.R
suppressPackageStartupMessages(library(terra))
cfg <- yaml::read_yaml("config.yml"); raw <- cfg$paths$data_dir
dir <- file.path(raw, "clima"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 7200)
box <- ext(-82, -60, -56, 13)

prov_f <- file.path(raw, "PROVENANCE.csv"); prov <- read.csv(prov_f)
log_prov <- function(f, u, lic) if (!f %in% prov$file)
  prov[nrow(prov) + 1, ] <<- list(f, u, format(Sys.Date()), digest::digest(file = file.path(raw, f), algo = "sha256"), lic)

# EarthEnv refuses GDAL range reads: fetch the whole file (0.7 GB), then crop
cu <- "https://data.earthenv.org/cloud/MODCF_meanannual.tif"; cf <- file.path(raw, basename(cu))
if (!file.exists(file.path(dir, "MODCF_meanannual.tif"))) {
  if (!file.exists(cf)) download.file(cu, cf, mode = "wb", quiet = TRUE)
  writeRaster(crop(rast(cf), box), file.path(dir, "MODCF_meanannual.tif"), datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
}
log_prov("clima/MODCF_meanannual.tif", paste(cu, "(cropped to -82,-60,-56,13)"), "EarthEnv cloud (CC BY-NC 4.0)")

# Global Aridity Index v3.1, annual; values x 10000 in the file
zf <- file.path(raw, "Global-AI_ET0__annual_v3_1.zip"); ai <- file.path(dir, "GlobalAI_v31_annual.tif")
if (!file.exists(ai)) {
  if (!file.exists(zf)) download.file("https://ndownloader.figshare.com/files/56300327", zf, mode = "wb", quiet = TRUE)
  tif <- grep("/ai_v31_yr\\.tif$", unzip(zf, list = TRUE)$Name, value = TRUE)   # exact name: the folder is also "AI"
  r <- crop(rast(paste0("/vsizip/", normalizePath(zf), "/", tif)), box) / 10000
  writeRaster(r, ai, datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
}
log_prov("Global-AI_ET0__annual_v3_1.zip", "https://ndownloader.figshare.com/files/56300327 (doi:10.6084/m9.figshare.7504448)",
         "CC BY 4.0")
write.csv(prov, prov_f, row.names = FALSE)
message("done")
