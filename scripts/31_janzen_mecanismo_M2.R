# 31_janzen_mecanismo_M2.R - M2.6: does the latitude interaction appear when the lowland tolerance is the flank's ANNUAL envelope
# and the exposure is the crest's MONTHLY regime (Janzen's premise), and is it absent only with month-to-month matching?
suppressPackageStartupMessages({library(terra); library(mgcv)})
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; raw <- cfg$paths$data_dir
wd <- file.path(raw, "clima_wc")
mon <- sprintf("%s_%02d", rep(c("tasmax", "tasmin"), each = 12), 1:12)
cl <- rast(file.path(wd, paste0("WC_", c("bio5", "bio6", mon), ".tif"))); names(cl) <- c("bio5", "bio6", mon)
tmpl <- cl[[1]]
dem <- resample(crop(rast(file.path(raw, "hyd_sa_dem_15s", "hyd_sa_dem_15s.tif")), tmpl), tmpl, "average")
dom <- vect(file.path(int, "M1_domains.gpkg"))
E <- dom[dom$dom == "E"]; I <- disagg(dom[dom$dom == "I"]); conn <- relate(I, E, "intersects")[, 1]
side <- rasterize(dom[dom$dom == "W"], tmpl, field = 1)
side <- cover(side, rasterize(E, tmpl, field = 2)); side <- cover(side, rasterize(I[conn], tmpl, field = 2))
side <- cover(side, rasterize(I[!conn], tmpl, field = 1))
nc <- ncol(tmpl); nr <- nrow(tmpl)
near_cells <- function(x, y, km, thr, sd) {
  k <- ceiling(km / (res(tmpl)[2] * 111.32)); rc <- rowColFromCell(tmpl, cellFromXY(tmpl, cbind(x, y)))
  kx <- ceiling(k / cos(y * pi / 180)); g <- expand.grid(dr = -k:k, dc = -kx:kx); rr <- rc[1] + g$dr; cc <- rc[2] + g$dc
  ok <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc; cell <- cellFromRowCol(tmpl, rr[ok], cc[ok]); xy <- xyFromCell(tmpl, cell)
  dk <- sqrt(((xy[, 1] - x) * cos(y * pi / 180))^2 + (xy[, 2] - y)^2) * 111.32
  e <- extract(dem, cell)[, 1]; s <- extract(side, cell)[, 1]
  cell[dk <= km & !is.na(e) & e <= thr & s %in% sd]
}
m2 <- read.csv("outputs/tables/M2_metricas_clima_collados.csv")
g <- subset(m2, gateway & !terminal & flank_m == 1000)
fl <- read.csv(file.path(int, "M2_flancos.csv")); fl <- fl[!is.na(fl$lon) & fl$flank_m == 1000, ]
fl <- merge(fl, g[, c("col_id", "side_lab", "lon", "lat", "elev_cross", "d_km", "TO_ann", "TO_mon")], by = "col_id", suffixes = c("_fl", ""))
fl <- fl[fl$side_lab == ifelse(fl$side == 1, "W", "E"), ]
V <- function(cells) colMeans(as.matrix(extract(cl, cells)), na.rm = TRUE)
res <- do.call(rbind, lapply(seq_len(nrow(fl)), function(j) {
  c5 <- near_cells(fl$lon_fl[j], fl$lat_fl[j], 5, 1000, fl$side[j]); if (!length(c5)) c5 <- cellFromXY(tmpl, cbind(fl$lon_fl[j], fl$lat_fl[j]))
  f <- V(c5); cr <- V(cellFromXY(tmpl, cbind(fl$lon[j], fl$lat[j])))
  mx <- sprintf("tasmax_%02d", 1:12); mn <- sprintf("tasmin_%02d", 1:12)
  # Janzen-premise index: crest monthly regime against the flank's annual envelope [bio6, bio5]
  d <- pmin(cr[mx], f["bio5"]) - pmax(cr[mn], f["bio6"])
  data.frame(fl[j, c("col_id", "side_lab", "lat", "elev_cross", "d_km", "TO_ann", "TO_mon")],
             TO_env = sum(d / sqrt((cr[mx] - cr[mn]) * (f["bio5"] - f["bio6"]))),   # <= 12, < 0 if disjoint
             months_in = sum(cr[mn] >= f["bio6"] & cr[mx] <= f["bio5"]),
             crest_seas = max(cr[mx]) - min(cr[mn]) - mean(cr[mx] - cr[mn]),
             flank_diurnal = mean(f[mx] - f[mn]), flank_annual = f["bio5"] - f["bio6"])
}))
w <- do.call(rbind, lapply(split(res, res$col_id), function(x) x[which.min(x$TO_env), ]))   # worst side
w$elev_km <- w$elev_cross / 1000; w$alat <- abs(w$lat)
fit <- function(y) { m <- gam(as.formula(paste(y, "~ elev_km * alat + s(d_km)")), data = w, method = "REML")
  s <- summary(m)$p.table["elev_km:alat", ]; ci <- s[1] + c(-1, 1) * 1.96 * s[2]
  sprintf("%-8s interaction %+.4f [%+.4f, %+.4f] p=%.2g dev=%.2f", y, s[1], ci[1], ci[2], s[4], summary(m)$dev.expl) }
# months_in is floored at 0 in the tropics, so a linear interaction on it is not interpretable; it is kept as a descriptor
res_fit <- sapply(c("TO_ann", "TO_mon", "TO_env"), fit); cat(res_fit, sep = "\n")
writeLines(res_fit, file.path(cfg$paths$outputs, "tables", "M2_janzen_mecanismo_modelos.txt"))
w$band <- cut(w$alat, c(0, 12, 23.5, 36, 60))
print(aggregate(cbind(flank_annual, flank_diurnal, months_in) ~ band, w, median))
t <- subset(w, alat <= 23.5 & lat <= 11)
cat("tropical cor(TO_ann, elev) =", round(cor(t$TO_ann, t$elev_cross, method = "spearman"), 2), " n =", nrow(t), "\n")
write.csv(w, file.path(cfg$paths$outputs, "tables", "M2_janzen_mecanismo.csv"), row.names = FALSE)
# Rscript scripts/31_janzen_mecanismo_M2.R
