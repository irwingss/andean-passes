# 08_marcadores.R - M1 (GLO-30 divide): marker / corridor raster for the marker-controlled flooding
# in 09_divisoria_glo30.py, and the list of GLO-30 tiles it needs.
# HydroBASINS only says which ocean each side drains to. Codes on the 15 s grid:
#   0 = outside the corridor (barrier), 1 = W seed, 2 = E seed, 3 = free (the DEM decides).
# Corridor = within 25 km of the HydroBASINS W|EI or I|E divide, plus the endorheic I basins connected
# to E (left free, so the flooding resolves the full W -> I -> E crossing). Seeds lie > 12 km from both divides.
# Rscript scripts/08_marcadores.R [--force]
suppressPackageStartupMessages({ library(sf); library(terra) })
sf_use_s2(FALSE)
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim
force <- "--force" %in% commandArgs(TRUE)
out <- file.path(int, "M1_markers_15s.tif")
if (file.exists(out) && !force) { message("exists, skip"); quit(save = "no") }
aea <- "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +datum=WGS84 +units=m"

dom <- vect(file.path(int, "M1_domains.gpkg"))
E <- dom[dom$dom == "E"]; Iall <- disagg(dom[dom$dom == "I"])
conn <- relate(Iall, E, "intersects")[, 1]
dom15 <- rast(file.path(int, "M1_domains_15s.tif"))            # grid template (-82..-60, -56..13)
code <- rasterize(dom[dom$dom == "W"], dom15, field = 1)
code <- cover(code, rasterize(E, dom15, field = 2))
code <- cover(code, rasterize(Iall[conn], dom15, field = 3))
code <- cover(code, rasterize(Iall[!conn], dom15, field = 1))   # coastal sinks (Sechura, Pisco): Pacific slope

# Distance to the nearest domain boundary (any W/E/I contact) on the raster. The HydroBASINS divide lines
# are clipped to GMBA + 20 km, which drops the W|I contact where it runs along the Coastal Cordillera
# (19.5-20.8 S, 23.6-25 S); without it the W side had no seeds there.
andes <- vect(st_transform(st_buffer(st_transform(read_sf(file.path(int, "M1_gmba_andes.gpkg")), aea), 100e3), 4326))
code <- mask(crop(code, andes), andes)
bnd <- ifel(focal(code, 3, "min", na.rm = TRUE) != focal(code, 3, "max", na.rm = TRUE), 1, NA)
# ponytail: distance on a 1' grid (terra::distance segfaults at 15 s); 1.85 km steps are fine for 12 / 25 km thresholds
dkm <- resample(distance(aggregate(bnd, 4, "max", na.rm = TRUE)) / 1000, code, method = "near")
corr <- !is.na(code) & (dkm <= 25 | code == 3)
m <- ifel(!corr, 0,
     ifel(code == 1 & dkm > 12, 1,
     ifel(code == 2 & dkm > 12, 2, 3)))
m <- extend(m, dom15, fill = 0)
writeRaster(m, out, datatype = "INT1U", overwrite = TRUE)

# GLO-30 tiles (1 x 1 deg) touched by the corridor
a <- aggregate(m > 0, fact = 240, fun = "max")                   # 240 x 15 s = 1 deg
xy <- xyFromCell(a, which(values(a) == 1))
tiles <- data.frame(la = floor(xy[, 2]), lo = floor(xy[, 1]))
tiles$name <- sprintf("Copernicus_DSM_COG_10_%s%02d_00_%s%03d_00_DEM", ifelse(tiles$la < 0, "S", "N"), abs(tiles$la),
                      ifelse(tiles$lo < 0, "W", "E"), abs(tiles$lo))
write.csv(tiles, file.path(int, "M1_glo30_tiles.csv"), row.names = FALSE)
f <- freq(m); print(f)
cat("tiles:", nrow(tiles), "\n")
