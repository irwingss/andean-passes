# 01_download_M1.R - inputs of the census, of the protection layers and of the supplement, into data_dir, with PROVENANCE.csv.
# Idempotent: skips files already present. Climate inputs: 13, 18, 19. GLO-30 tiles: 07. Rscript scripts/01_download_M1.R
cfg <- yaml::read_yaml("config.yml")
raw <- cfg$paths$data_dir
dir.create(raw, recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)

src <- data.frame(
  file = c("hybas_sa_lev01-12_v1c.zip", "hyd_sa_dem_15s.zip", "hyd_sa_dir_15s.zip", "GMBA_Inventory_v2.0_standard_300.zip",
           "ghm/gHM.zip", "refs/Hausfather2022_supporting.xlsx"),
  url = c("https://data.hydrosheds.org/file/HydroBASINS/standard/hybas_sa_lev01-12_v1c.zip",
          "https://data.hydrosheds.org/file/hydrosheds-v1-dem/hyd_sa_dem_15s.zip",
          "https://data.hydrosheds.org/file/hydrosheds-v1-dir/hyd_sa_dir_15s.zip",
          "https://data.earthenv.org/mountains/standard/GMBA_Inventory_v2.0_standard_300.zip",
          "https://ndownloader.figshare.com/files/13448294",
          "https://zenodo.org/api/records/6476375/files/Hausfather%20et%20al%20Nature%202022%20supporting%20data.xlsx/content"),
  license = c(rep("HydroSHEDS License (free, attribution)", 3), "CC BY 4.0", "CC BY 4.0", "CC BY 4.0"))

prov_f <- file.path(raw, "PROVENANCE.csv")
prov <- if (file.exists(prov_f)) read.csv(prov_f) else
  data.frame(file = character(), url = character(), date = character(), sha256 = character(), license = character())

for (i in seq_len(nrow(src))) {
  dest <- file.path(raw, src$file[i])
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(dest)) download.file(src$url[i], dest, mode = "wb", quiet = TRUE)
  if (!src$file[i] %in% prov$file) {
    prov[nrow(prov) + 1, ] <- list(src$file[i], src$url[i], format(Sys.Date()),
                                   digest::digest(file = dest, algo = "sha256"), src$license[i])
  }
  if (grepl("\\.zip$", dest)) unzip(dest, exdir = if (grepl("gHM", dest)) dirname(dest) else sub("\\.zip$", "", dest), overwrite = FALSE)
}
write.csv(prov, prov_f, row.names = FALSE)

# GADM 4.1 level 1 (map outlines only): GADM forbids redistribution, so it is fetched here; files land in gadm/gadm/
for (cc in c("BOL", "CHL", "COL", "ECU", "PER", "VEN")) geodata::gadm(cc, level = 1, path = file.path(raw, "gadm"), version = "4.1")

# WDPA: its terms must be accepted on the Protected Planet site, and they forbid redistribution
w <- file.path(raw, "wdpa", sprintf("WDPA_Nov2026_Public_shp_%d.zip", 0:2))
if (!all(file.exists(w)))
  message("MANUAL STEP. Download the global WDPA shapefile, November 2026 release, from https://www.protectedplanet.net/en/thematic-areas/wdpa ",
          "and extract its parts WDPA_Nov2026_Public_shp_0.zip, _1.zip and _2.zip into ", file.path(raw, "wdpa"), ". ",
          "A later release gives slightly different protection results; 30_proteccion_M3.R reads this file name.")
print(prov[, c("file", "date", "license")])
