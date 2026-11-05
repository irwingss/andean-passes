# 19_descarga_M3.R - M3 inputs for future permeability, as deltas on the WorldClim 2.1 30" baseline (plan v1.4).
#   WorldClim 2.1 30" monthly prec tiles (tmax/tmin already in 18_descarga_worldclim_M2.R)
#   WorldClim 2.1 10' baseline tmax/tmin/prec + CMIP6 10' for the 8 screened GCMs x SSP126/370/585 x 2041-2100
# Rscript scripts/19_descarga_M3.R
suppressPackageStartupMessages(library(terra))
cfg <- yaml::read_yaml("config.yml"); raw <- cfg$paths$data_dir
options(timeout = 7200)
wd <- file.path(raw, "clima_wc"); fd <- file.path(raw, "cmip6")
dir.create(fd, recursive = TRUE, showWarnings = FALSE)

# 30" monthly precipitation, cropped like the M2 layers
f <- file.path(wd, sprintf("WC_prec_%02d.tif", 1:12))
if (!all(file.exists(f))) {
  t <- lapply(c(5, -15, -45), function(lat) geodata::worldclim_tile("prec", -75, lat, path = wd))
  r <- crop(merge(sprc(t)), ext(-82, -60, -56, 13))
  for (k in 1:12) writeRaster(r[[k]], f[k], overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
}

# 10' baseline and CMIP6 (geodata caches; skips existing files)
GCM <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-ESM2", "CNRM-ESM2-1", "GISS-E2-1-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")
SSP <- c("126", "370", "585"); PER <- c("2041-2060", "2061-2080", "2081-2100")
for (v in c("tmax", "tmin", "prec")) geodata::worldclim_global(v, 10, path = fd)
for (m in GCM) for (s in SSP) for (p in PER) for (v in c("tmax", "tmin", "prec"))
  if (is.null(geodata::cmip6_world(m, s, p, v, 10, path = fd))) message("MISSING ", m, " ", s, " ", p, " ", v)

prov_f <- file.path(raw, "PROVENANCE.csv"); prov <- read.csv(prov_f)
add <- function(k, url, lic, files) if (!k %in% prov$file) prov[nrow(prov) + 1, ] <<- list(k, url, format(Sys.Date()),
  digest::digest(paste(sapply(files, digest::digest, file = TRUE, algo = "sha256")), algo = "sha256"), lic)
add("clima_wc/WC_prec_MM.tif (WorldClim 2.1 30s 'prec' tiles; cropped)", "https://geodata.ucdavis.edu/climate/worldclim/2_1/tiles/tile/",
    "WorldClim 2.1 (Fick & Hijmans 2017)", f)
add("cmip6/ (WorldClim 2.1 10m baseline + CMIP6 10m, 8 GCM x 3 SSP x 3 periods)", "https://geodata.ucdavis.edu/cmip6/10m/",
    "WorldClim 2.1 (Fick & Hijmans 2017)", list.files(fd, "\\.tif$", recursive = TRUE, full.names = TRUE))
write.csv(prov, prov_f, row.names = FALSE)
message("done")
