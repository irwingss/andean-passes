# 15_metricas_clima_M2.R - M2.2-M2.5: crest-flank climate metrics for every col with flanks (14_flancos_M2.R).
# Climate: WorldClim 2.1 (1970-2000, 30"; 18_descarga_worldclim_M2.R), plus Global-AI v3.1 (built on WorldClim 2.1)
# and MODIS cloud frequency (13_descarga_clima_M2.R), all on the WorldClim grid.
# Crest = 30" cell of the col. Flank = mean of the 30" cells <= threshold within 5 km of the flank point
# (the flank cell itself if none). Metrics per col x side x threshold:
#   TO_ann   Janzen overlap of annual ranges [bio6, bio5]: d / sqrt(R1 R2), d = interval overlap (<= 1, < 0 if disjoint)
#   TO_mon   Zuloaga & Kerr (2016) Eq. 1 with monthly [tasmin, tasmax]: sum_i d_i / sqrt(R1i R2i) (<= 12)
#   AI_gap   log(AI crest / AI flank), Global-AI v3.1
#   CLD_gap  MODIS cloud frequency crest - flank (%)
# Gateways only (default threshold): Mahalanobis D2 of the crest climate against the lowland cells (<= threshold,
# same side) within 100 km of the flank, on variables with |r| < 0.7; reported as the percentile of the crest D2
# among the lowland cells (1 = more extreme than every lowland cell).
# Rscript scripts/15_metricas_clima_M2.R [--force]
suppressPackageStartupMessages(library(terra))
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; raw <- cfg$paths$data_dir; P <- cfg$params
force <- "--force" %in% commandArgs(TRUE)
out <- file.path(cfg$paths$outputs, "tables", "M2_metricas_clima_collados.csv")
if (file.exists(out) && !force) { message("exists, skip"); quit(save = "no") }
set.seed(cfg$seed)

cd <- file.path(raw, "clima"); wd <- file.path(raw, "clima_wc")
bio <- c("bio1", "bio5", "bio6", "bio12", "bio14", "bio15")
mon <- sprintf("%s_%02d", rep(c("tasmax", "tasmin"), each = 12), 1:12)
cl <- rast(c(file.path(wd, paste0("WC_", c(bio, mon), ".tif")),
             file.path(cd, c("GlobalAI_v31_annual.tif", "MODCF_meanannual.tif"))))
names(cl) <- c(bio, mon, "ai", "cloud")
tmpl <- cl[[1]]
cl$cloud <- cl$cloud / 100                                  # MODCF stored as % x 100

# 30" elevation (HydroSHEDS 15 s averaged) and side (1 W, 2 E) on the WorldClim grid
dem <- resample(crop(rast(file.path(raw, "hyd_sa_dem_15s", "hyd_sa_dem_15s.tif")), tmpl), tmpl, "average")
dom <- vect(file.path(int, "M1_domains.gpkg"))
E <- dom[dom$dom == "E"]; I <- disagg(dom[dom$dom == "I"]); conn <- relate(I, E, "intersects")[, 1]
side <- rasterize(dom[dom$dom == "W"], tmpl, field = 1)
side <- cover(side, rasterize(E, tmpl, field = 2))
side <- cover(side, rasterize(I[conn], tmpl, field = 2))
side <- cover(side, rasterize(I[!conn], tmpl, field = 1))
ev <- values(dem)[, 1]; sv <- values(side)[, 1]; nc <- ncol(tmpl); nr <- nrow(tmpl)

# Cells of side `sd` with elevation <= thr within `km` of (x, y)
near_cells <- function(x, y, km, thr, sd) {
  k <- ceiling(km / (res(tmpl)[2] * 111.32)); rc <- rowColFromCell(tmpl, cellFromXY(tmpl, cbind(x, y)))
  kx <- ceiling(k / cos(y * pi / 180))
  g <- expand.grid(dr = -k:k, dc = -kx:kx); rr <- rc[1] + g$dr; cc <- rc[2] + g$dc
  ok <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc
  cell <- cellFromRowCol(tmpl, rr[ok], cc[ok]); xy <- xyFromCell(tmpl, cell)
  dk <- sqrt(((xy[, 1] - x) * cos(y * pi / 180))^2 + (xy[, 2] - y)^2) * 111.32
  cell[dk <= km & !is.na(ev[cell]) & ev[cell] <= thr & sv[cell] %in% sd]
}

cols <- read.csv(file.path(cfg$paths$outputs, "database", "APD_cols.csv"))
fl <- read.csv(file.path(int, "M2_flancos.csv")); fl <- fl[!is.na(fl$lon), ]
fl$cells <- lapply(seq_len(nrow(fl)), function(j) {
  c5 <- near_cells(fl$lon[j], fl$lat[j], P$flank_buffer_km, fl$flank_m[j], fl$side[j])
  if (length(c5)) c5 else cellFromXY(tmpl, cbind(fl$lon[j], fl$lat[j]))
})
fl$n_cells <- lengths(fl$cells)
cc <- cols[cols$col_id %in% fl$col_id, ]
crest_cell <- setNames(cellFromXY(tmpl, cbind(cc$lon, cc$lat)), cc$col_id)

uc <- unique(c(crest_cell, unlist(fl$cells)))
V <- as.matrix(extract(cl, uc)); rownames(V) <- uc
crest <- V[as.character(crest_cell[fl$col_id]), , drop = FALSE]
flank <- t(vapply(fl$cells, function(x) colMeans(V[as.character(x), , drop = FALSE], na.rm = TRUE), numeric(ncol(V))))

ovl <- function(lo1, hi1, lo2, hi2) (pmin(hi1, hi2) - pmax(lo1, lo2)) / sqrt((hi1 - lo1) * (hi2 - lo2))
fl$TO_ann <- ovl(crest[, "bio6"], crest[, "bio5"], flank[, "bio6"], flank[, "bio5"])
mx <- sprintf("tasmax_%02d", 1:12); mn <- sprintf("tasmin_%02d", 1:12)
fl$TO_mon <- rowSums(ovl(crest[, mn], crest[, mx], flank[, mn], flank[, mx]))
fl$AI_gap <- log(crest[, "ai"] / flank[, "ai"])
fl$CLD_gap <- crest[, "cloud"] - flank[, "cloud"]
for (v in c("bio1", "bio5", "bio6", "bio12", "ai", "cloud")) {
  fl[[paste0("crest_", v)]] <- crest[, v]; fl[[paste0("flank_", v)]] <- flank[, v]
}

# Mahalanobis for gateways at the default threshold
g <- which(fl$col_id %in% cols$col_id[cols$gateway] & fl$flank_m == P$flank_m[1])
mvars0 <- c("bio5", "bio6", "bio12", "bio15", "cloud", "ai", "bio1", "bio14")   # priority order
low <- lapply(g, function(j) near_cells(fl$lon[j], fl$lat[j], 100, fl$flank_m[j], fl$side[j]))
pool <- as.matrix(extract(cl[[mvars0]], unique(unlist(lapply(low, function(x) x[sample.int(length(x), min(2000, length(x)))])))))
r <- abs(cor(pool, use = "pairwise.complete.obs")); mv <- character(0)
for (v in mvars0) if (all(r[v, mv] < 0.7)) mv <- c(mv, v)
message("Mahalanobis variables: ", paste(mv, collapse = ", "))
fl$maha_pct <- NA_real_; fl$maha_n <- NA_integer_
for (k in seq_along(g)) {
  j <- g[k]; X <- na.omit(as.matrix(extract(cl[[mv]], low[[k]])))
  if (nrow(X) < 10 * length(mv)) next
  mu <- colMeans(X); S <- cov(X)
  d0 <- mahalanobis(crest[j, mv], mu, S); dx <- mahalanobis(X, mu, S)
  fl$maha_pct[j] <- mean(dx < d0); fl$maha_n[j] <- nrow(X)
}

res_ <- fl[, setdiff(names(fl), "cells")]
res_ <- merge(cols[, c("col_id", "name", "lon", "lat", "elev_cross", "depth", "d_km", "tropical", "gateway", "terminal", "lara")],
              res_, by = "col_id", suffixes = c("", "_flank"))
write.csv(res_, out, row.names = FALSE)
writeLines(mv, file.path(int, "M2_maha_vars.txt"))
print(res_[res_$name %in% "Abra de Porculla", c("side_lab", "flank_m", "TO_ann", "TO_mon", "AI_gap", "CLD_gap", "maha_pct")])
