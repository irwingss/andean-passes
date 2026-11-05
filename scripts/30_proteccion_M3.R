# 30_proteccion_M3.R - M3.6: protection and human pressure on gateways and on the montane corridor (plan v2.2).
# Rules fixed in BITACORA 2026-11-03 before computing.
# WDPA Nov 2026 (UNEP-WCMC & IUCN, Protected Planet): polygons, SITE_TYPE = PA, STATUS Designated / Inscribed /
#   Established, REALM not Marine, UNESCO-MAB Biosphere Reserves excluded (as in Protected Planet coverage statistics).
#   Licence forbids redistribution: only derived attributes leave this script.
# Global Human Modification (Kennedy et al. 2019, 1 km Mollweide, CC BY 4.0), projected to the WorldClim 30" grid.
# Gateway (30 tropical, flank 1000 m): crest in PA; % of route length (both sides, D8, no fallback) in PA; gHM mean
#   within 1 km of the crest and along the routes. Divide bin: in PA, gHM.
# Priority: critical if (a) fraction of dT* >= 0.5 in 2081-2100 under SSP3-7.0 (E-W lowland corridor opening) or
#   (b) robust montane-corridor loss (SSP3-7.0, 2081-2100) within +-25 km of divide; unprotected if the crest is
#   outside every PA; high pressure if crest gHM >= 0.4.
# Rscript scripts/30_proteccion_M3.R
suppressPackageStartupMessages({ library(sf); library(terra) })
sf_use_s2(FALSE)
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; raw <- cfg$paths$data_dir
tab <- file.path(cfg$paths$outputs, "tables")
bb <- st_as_sfc(st_bbox(c(xmin = -82, ymin = -25, xmax = -60, ymax = 13), crs = 4326))   # tropical Andes + margin

# ---- clipped layers (cached in interim) ----
pa_f <- file.path(int, "M3_wdpa_andes.gpkg")
if (!file.exists(pa_f)) {
  pa <- do.call(rbind, lapply(0:2, function(i) {
    x <- st_read(sprintf("/vsizip/%s/wdpa/WDPA_Nov2026_Public_shp_%d.zip/WDPA_Nov2026_Public_shp-polygons.shp", raw, i),
                 wkt_filter = st_as_text(bb), quiet = TRUE)
    x[x$SITE_TYPE == "PA" & x$STATUS %in% c("Designated", "Inscribed", "Established") & x$REALM != "Marine" &
        x$DESIG_ENG != "UNESCO-MAB Biosphere Reserve", c("SITE_ID", "NAME", "DESIG_ENG", "IUCN_CAT", "ISO3", "STATUS_YR")]
  }))
  pa <- st_make_valid(pa)
  st_write(pa, pa_f, delete_dsn = TRUE, quiet = TRUE)
}
pa <- st_read(pa_f, quiet = TRUE)
gh_f <- file.path(int, "M3_ghm_andes.tif")
if (!file.exists(gh_f)) {
  tm <- crop(rast(file.path(raw, "clima_wc", "WC_bio1.tif")), ext(-82, -60, -25, 13))
  g <- rast(file.path(raw, "ghm", "gHM", "gHM.tif"))
  g <- crop(g, project(ext(tm), "EPSG:4326", crs(g)) * 1.1)
  writeRaster(project(g, tm, method = "bilinear"), gh_f, overwrite = TRUE)
}
gh <- rast(gh_f)
message(nrow(pa), " protected areas")
pa_u <- st_union(st_geometry(pa))

# ---- gateways ----
h3 <- read.csv(file.path(tab, "M2_H3_tipologia.csv"))
cols <- read.csv("outputs/database/APD_cols.csv", encoding = "UTF-8")
g <- merge(h3[, c("col_id", "lon", "lat", "elev_cross", "h3_type")], cols[, c("col_id", "name", "d_km")])   # APD names
gp <- st_as_sf(g, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
g$crest_in_PA <- lengths(st_intersects(gp, pa)) > 0
b5 <- buffer(vect(gp), 5000); b5$col_id <- g$col_id                                      # sensitivity: share of 5 km disc in PA
i5 <- intersect(b5, vect(pa_u))
g$pct_PA_5km <- 0; if (nrow(i5)) g$pct_PA_5km[match(i5$col_id, g$col_id)] <- 100 * expanse(i5) / expanse(b5)[match(i5$col_id, g$col_id)]
g$ghm_crest <- extract(gh, buffer(vect(gp), 1000), mean, na.rm = TRUE)[, 2]               # terra: geodesic metres
rt <- st_read(file.path(int, "M2_flank_routes.gpkg"), quiet = TRUE)
rt <- rt[rt$col_id %in% g$col_id, ]; rt$col_id_side <- paste(rt$col_id, rt$side)
rv <- vect(rt); rt$len <- perim(rv)                                                        # geodesic length (m)
ri <- intersect(rv, vect(pa_u))
rt$len_pa <- 0
if (nrow(ri)) { a <- tapply(perim(ri), ri$col_id_side, sum); rt$len_pa <- ifelse(is.na(a[rt$col_id_side]), 0, a[rt$col_id_side]) }
rt$ghm <- extract(gh, rv, mean, na.rm = TRUE)[, 2]
ra <- aggregate(cbind(len, len_pa) ~ col_id, st_drop_geometry(rt), sum)
g$routes_pct_PA <- 100 * ra$len_pa[match(g$col_id, ra$col_id)] / ra$len[match(g$col_id, ra$col_id)]
ga <- aggregate(ghm ~ col_id, st_drop_geometry(rt), mean); g$ghm_routes <- ga$ghm[match(g$col_id, ga$col_id)]

# ---- divide bins (tropical) ----
b <- read.csv(file.path(int, "M3_corredor_bins.csv"), check.names = FALSE)
b <- b[b$tropical, ]
bp <- st_as_sf(b, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
b$in_PA <- lengths(st_intersects(bp, pa)) > 0
b$ghm <- extract(gh, cbind(b$lon, b$lat))[, 1]
write.csv(b, file.path(int, "M3_corredor_bins_proteccion.csv"), row.names = FALSE)

# ---- priority ----
ens <- read.csv(file.path(tab, "M3_futuro_ensamble.csv"))
e90 <- ens[ens$ssp == 370 & ens$period == "2081-2100", ]
g$frac_dTstar_370_2090 <- e90$frac_dTstar.mean[match(g$col_id, e90$col_id)]
g$opening <- !is.na(g$frac_dTstar_370_2090) & g$frac_dTstar_370_2090 >= 0.5
g$montane_loss_km <- vapply(g$d_km, function(d) sum(b$loss_robust_370_2090[abs(b$d_km - d) <= 25]) * 0.5, 0)
g$montane_bottleneck <- g$montane_loss_km > 0
g$critical <- g$opening | g$montane_bottleneck
g$critical_unprotected <- g$critical & !g$crest_in_PA
g$high_pressure <- g$ghm_crest >= 0.4
write.csv(g[order(g$lat, decreasing = TRUE), ], file.path(tab, "M3_prioridad_puertas.csv"), row.names = FALSE)

km <- function(x) sum(x) * 0.5
cor_sum <- data.frame(
  stretch = c("montane corridor today", "robust loss by 2081-2100 (SSP3-7.0)", "montane corridor kept"),
  km = c(km(!b$gap_now), km(b$loss_robust_370_2090), km(!b$gap_now & !b$loss_robust_370_2090)),
  pct_in_PA = c(100 * mean(b$in_PA[!b$gap_now]), 100 * mean(b$in_PA[b$loss_robust_370_2090]),
                100 * mean(b$in_PA[!b$gap_now & !b$loss_robust_370_2090])),
  ghm_mean = c(mean(b$ghm[!b$gap_now], na.rm = TRUE), mean(b$ghm[b$loss_robust_370_2090], na.rm = TRUE),
               mean(b$ghm[!b$gap_now & !b$loss_robust_370_2090], na.rm = TRUE)),
  pct_high_pressure = c(100 * mean(b$ghm[!b$gap_now] >= 0.4, na.rm = TRUE), 100 * mean(b$ghm[b$loss_robust_370_2090] >= 0.4, na.rm = TRUE),
                        100 * mean(b$ghm[!b$gap_now & !b$loss_robust_370_2090] >= 0.4, na.rm = TRUE)))
write.csv(cor_sum, file.path(tab, "M3_corredor_proteccion.csv"), row.names = FALSE)
print(cor_sum, digits = 3)

# sensitivity to the montane threshold; robust loss by country (divide bins on a border go to the country map.where returns)
lt <- read.csv(file.path(int, "M3_perdida_por_umbral.csv"))
lt <- lt[match(b$bin, lt$bin), ]
sens <- do.call(rbind, lapply(grep("^loss_", names(lt), value = TRUE), function(v) {
  l <- lt[[v]]; lk <- vapply(g$d_km, function(d) sum(l[abs(b$d_km - d) <= 25]), 0) > 0; crit <- g$opening | lk
  data.frame(threshold_m = as.numeric(sub("loss_", "", v)), loss_km = km(l), n_critical = sum(crit), n_both = sum(g$opening & lk),
             critical = paste(ifelse(is.na(g$name), "Táchira depression", g$name)[crit][order(-g$lat[crit])], collapse = "; "), pct_lost_in_PA = 100 * mean(b$in_PA[l]),
             pct_lost_high = 100 * mean(b$ghm[l] >= 0.4, na.rm = TRUE))
}))
write.csv(sens, file.path(tab, "M3_corredor_umbral.csv"), row.names = FALSE)
pais <- maps::map.where("world", b$lon, b$lat); na <- is.na(pais); pais[na] <- maps::map.where("world", b$lon[na] + 0.05, b$lat[na])
b$pais <- sub(":.*", "", pais)
pz <- do.call(rbind, lapply(split(b, b$pais), function(x) data.frame(country = x$pais[1], corridor_km = km(!x$gap_now), lost_km = km(x$loss_robust_370_2090))))
pz$pct_lost <- 100 * pz$lost_km / pz$corridor_km
write.csv(pz, file.path(tab, "M3_perdida_pais.csv"), row.names = FALSE)
print(sens[, -5], digits = 3); print(pz, digits = 3)
print(g[order(g$lat, decreasing = TRUE), c("name", "lat", "elev_cross", "h3_type", "crest_in_PA", "routes_pct_PA", "ghm_crest",
                                            "frac_dTstar_370_2090", "montane_loss_km", "critical", "critical_unprotected", "high_pressure")], digits = 2)
