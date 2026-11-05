# 02_dominios_divisoria.R - M1.1-1.2: classify HydroBASINS lev12 into W / E / I domains
# by terminal-basin outlet, then extract the W | (E u I) and I | E divides clipped to GMBA Andes + 20 km.
# Rscript scripts/02_dominios_divisoria.R [--force]
suppressPackageStartupMessages({ library(sf); library(terra); library(ggplot2) })
sf_use_s2(FALSE)
cfg <- yaml::read_yaml("config.yml"); raw <- cfg$paths$data_dir; int <- cfg$paths$interim
force <- "--force" %in% commandArgs(TRUE)
out_dom <- file.path(int, "M1_domains.gpkg"); out_div <- file.path(int, "M1_divide_raw.gpkg")
if (file.exists(out_div) && !force) { message("exists, skip (use --force)"); quit(save = "no") }

# --- Andes mask: GMBA v2 units under "South America > Andes", + margin ---
g <- read_sf(file.path(raw, "GMBA_Inventory_v2.0_standard_300/GMBA_Inventory_v2.0_standard_300.shp"))
andes <- st_union(st_make_valid(g[grepl("> Andes", g$Path), ]))
aea <- "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +datum=WGS84 +units=m"  # SA Albers
andes_m <- st_transform(andes, aea)
mask20 <- st_transform(st_buffer(andes_m, cfg$params$gmba_margin_km * 1000), 4326)
mask200 <- st_transform(st_buffer(andes_m, 200e3), 4326)  # working extent for basins

# --- Domains ---
if (file.exists(out_dom) && !force) dom <- vect(out_dom) else {
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
}

# --- Divides: shared boundaries between dissolved domains ---
W <- dom[dom$dom == "W"]; E <- dom[dom$dom == "E"]
# Only endorheic basins connected to E belong to the barrier (drops coastal sinks: Sechura, Pisco)
I <- disagg(dom[dom$dom == "I"]); I <- I[relate(I, E, "intersects")[, 1]]
I <- aggregate(I); I$dom <- "I"; EI <- aggregate(rbind(E[, 0], I[, 0]))
shared <- function(a, z) {
  l <- sharedPaths(rbind(a[, 0], z[, 0]))                 # edges common to both polygons
  st_line_merge(st_union(st_as_sf(l)))
}
to_lines <- function(x) st_cast(st_cast(st_collection_extract(x, "LINESTRING"), "MULTILINESTRING"), "LINESTRING")
div_w <- to_lines(st_intersection(shared(W, EI), mask20))
div_ie <- to_lines(st_intersection(shared(I, E), mask20))
div <- rbind(st_sf(type = "W|EI", geometry = div_w), st_sf(type = "I|E", geometry = div_ie))
div$len_km <- as.numeric(st_length(st_transform(div, aea))) / 1000
st_write(div, out_div, delete_dsn = TRUE, quiet = TRUE)
st_write(st_sf(geometry = andes), file.path(int, "M1_gmba_andes.gpkg"), delete_dsn = TRUE, quiet = TRUE)
print(aggregate(len_km ~ type, st_drop_geometry(div), function(x) c(n = length(x), total = sum(x), max = max(x))))

# --- QC map ---
dsf <- st_as_sf(dom)
p <- ggplot() + geom_sf(data = dsf, aes(fill = dom), colour = NA) +
  scale_fill_manual(values = c(W = "#E69F00", E = "#56B4E9", I = "#CC79A7")) +
  geom_sf(data = st_sf(geometry = andes), fill = NA, colour = "grey30", linewidth = 0.2) +
  geom_sf(data = div, aes(colour = type), linewidth = 0.4) +
  scale_colour_manual(values = c("W|EI" = "black", "I|E" = "red")) +
  coord_sf(xlim = c(-82, -62), ylim = c(-56, 12)) + theme_bw(8)
dir.create("outputs/figs_qc", showWarnings = FALSE, recursive = TRUE)
ggsave("outputs/figs_qc/M1_dominios_hydrobasins.png", p, width = 120, height = 260, units = "mm", dpi = 200)
pn <- p + coord_sf(xlim = c(-80, -66), ylim = c(-8, 12))
ggsave("outputs/figs_qc/M1_dominios_hydrobasins_norte.png", pn, width = 160, height = 200, units = "mm", dpi = 200)
