# exp_paises.R - EXPERIMENT: country outlines for the maps, from the GADM 4.1 level-1 files already in data/raw.
# Rscript scripts/exp_paises.R   ->  data/interim/M1_paises.gpkg
suppressMessages(library(terra))
f <- list.files("data/raw/gadm/gadm", pattern="_pk.rds$", full.names=TRUE)
rd <- function(p) { x <- readRDS(p); if (inherits(x, "PackedSpatVector")) x <- vect(x); aggregate(x, by="COUNTRY") }
v <- do.call(rbind, lapply(f, rd))
v <- crop(v, ext(-85, -55, -30, 15))
writeVector(v, "data/interim/M1_paises.gpkg", overwrite=TRUE)
cat(nrow(v), paste(v$COUNTRY, collapse=", "), "\n")
