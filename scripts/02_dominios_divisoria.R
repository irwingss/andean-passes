# 02_dominios_divisoria.R - Study area and sea labels of the census (PRD O1 and O4).
# Writes data/interim/M1_gmba_andes.gpkg, the union of the GMBA v2 units under "South America > Andes", and
# data/interim/M1_domains.gpkg, the level-12 HydroBASINS dissolved by terminal outlet into W (Pacific and Caribbean),
# E (Atlantic) and I (endorheic). 41_etiquetas_passes.py reads these domains only to label the side each pass drains to.
# Rscript scripts/02_dominios_divisoria.R [--force]
suppressPackageStartupMessages({ library(sf); library(terra) })
sf_use_s2(FALSE)
cfg <- yaml::read_yaml("config.yml"); raw <- cfg$paths$data_dir; int <- cfg$paths$interim
force <- "--force" %in% commandArgs(TRUE)
out_dom <- file.path(int, "M1_domains.gpkg")
if (file.exists(out_dom) && !force) { message("exists, skip (use --force)"); quit(save = "no") }

# --- Andes mask: GMBA v2 units under "South America > Andes", + margin ---
g <- read_sf(file.path(raw, "GMBA_Inventory_v2.0_standard_300/GMBA_Inventory_v2.0_standard_300.shp"))
andes <- st_union(st_make_valid(g[grepl("> Andes", g$Path), ]))
aea <- "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +datum=WGS84 +units=m"  # SA Albers
andes_m <- st_transform(andes, aea)
mask200 <- st_transform(st_buffer(andes_m, 200e3), 4326)  # working extent for basins

# --- Domains ---
b <- read_sf(file.path(raw, "hybas_sa_lev01-12_v1c/hybas_sa_lev12_v1c.shp"))
outl <- b[b$HYBAS_ID == b$MAIN_BAS, ]
xy <- st_coordinates(st_point_on_surface(outl))
# ponytail: ocean assignment by outlet coordinates (rules logged in BITACORA); fine because only
# the divide inside the Andes matters and all Andean outlets are far from the rule boundaries.
cls <- ifelse(outl$ENDO == 2, "I",
       ifelse(xy[, 2] > 7, ifelse(xy[, 1] < -62, "W", "E"),
       ifelse(xy[, 2] >= -41, ifelse(xy[, 1] > -66, "E", "W"),
                              ifelse(xy[, 1] > -69.8, "E", "W"))))
b$dom <- cls[match(b$MAIN_BAS, outl$HYBAS_ID)]
stopifnot(!anyNA(b$dom))
b <- b[lengths(st_intersects(b, mask200)) > 0, ]
dom <- aggregate(vect(b)[, "dom"], by = "dom")           # terra dissolve: 3 polygons W, E, I
writeVector(dom, out_dom, overwrite = TRUE)
st_write(st_sf(geometry = andes), file.path(int, "M1_gmba_andes.gpkg"), delete_dsn = TRUE, quiet = TRUE)
