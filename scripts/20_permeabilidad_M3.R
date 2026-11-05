# 20_permeabilidad_M3.R - M3: crest-flank permeability under CMIP6 futures (H4), delta method on WorldClim 2.1 30".
# Units: tropical + extratropical gateways (no terminal channels) and the Huancabamba system, flank threshold 1000 m,
# same crest cell / flank cells as 15_metricas_clima_M2.R.
# Delta per scenario (BITACORA 2026-11-02): WorldClim 2.1 CMIP6 10' minus WorldClim 2.1 10' baseline, monthly additive
# dTmax and dTmin, monthly ratio (P_f + 1) / (P_0 + 1); bilinear from 10' to each 30" cell (coastal NA filled by focal means).
# Metrics per col x side x scenario (bio5 = max tmax, bio6 = min tmin, P = annual sum):
#   TO_fix  Janzen overlap crest(future) vs flank(present) (fixed lowland niche)
#   Pg_fix  log(P crest(future) / P flank(present))
#   dTstar       uniform warming that centres the crest range on the present flank range (max TO_fix); per side,
#                scenario-free (M3_dTstar.csv). frac_dTstar = crest warming in the scenario / dTstar.
#   over_T       crest bio5 or bio6 above every present lowland cell (<= threshold, same side, 100 km): TO_fix falls by
#                overshoot, not closure.
# crest(future) vs flank(future) is not computed: a shared coarse delta leaves it near-invariant by construction.
# Col level = worst side (min TO_fix and frac_dTstar, max |Pg_fix| and dTstar, any flag), as in 17_janzen_M2.R: at col level
# Pg_fix is a magnitude, so dPg > 0 means the precipitation contrast widens, whichever side is wetter.
# Rscript scripts/20_permeabilidad_M3.R [--force]
suppressPackageStartupMessages(library(terra))
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; raw <- cfg$paths$data_dir; P <- cfg$params
tab <- file.path(cfg$paths$outputs, "tables"); force <- "--force" %in% commandArgs(TRUE)
outs <- file.path(tab, c("M3_futuro.csv", "M3_futuro_ensamble.csv", "M3_dTstar.csv"))
if (all(file.exists(outs)) && !force) { message("exists, skip"); quit(save = "no") }
wd <- file.path(raw, "clima_wc"); fd <- file.path(raw, "cmip6")
box <- ext(-82, -60, -56, 13)

# ---- units and cells (as in 15_metricas_clima_M2.R) ----
m2 <- read.csv(file.path(tab, "M2_metricas_clima_collados.csv"))
H <- P$huancabamba
m2$huanca <- m2$lat >= H$lat[1] & m2$lat <= H$lat[2] & m2$elev_cross < H$max_elev_m & m2$depth >= H$min_depth_m
u <- m2[!m2$terminal & (m2$gateway | m2$huanca) & m2$flank_m == P$flank_m[1] & !is.na(m2$TO_ann), ]

tx <- sprintf("tasmax_%02d", 1:12); tn <- sprintf("tasmin_%02d", 1:12); pr <- sprintf("prec_%02d", 1:12)
base <- rast(file.path(wd, paste0("WC_", c(tx, tn, pr), ".tif"))); names(base) <- c(tx, tn, pr)
tmpl <- base[[1]]
dem <- resample(crop(rast(file.path(raw, "hyd_sa_dem_15s", "hyd_sa_dem_15s.tif")), tmpl), tmpl, "average")
dom <- vect(file.path(int, "M1_domains.gpkg"))
E <- dom[dom$dom == "E"]; I <- disagg(dom[dom$dom == "I"]); conn <- relate(I, E, "intersects")[, 1]
side <- rasterize(dom[dom$dom == "W"], tmpl, field = 1)
side <- cover(side, rasterize(E, tmpl, field = 2))
side <- cover(side, rasterize(I[conn], tmpl, field = 2))
side <- cover(side, rasterize(I[!conn], tmpl, field = 1))
ev <- values(dem)[, 1]; sv <- values(side)[, 1]; nc <- ncol(tmpl); nr <- nrow(tmpl)
near_cells <- function(x, y, km, thr, sd) {
  k <- ceiling(km / (res(tmpl)[2] * 111.32)); rc <- rowColFromCell(tmpl, cellFromXY(tmpl, cbind(x, y)))
  kx <- ceiling(k / cos(y * pi / 180))
  g <- expand.grid(dr = -k:k, dc = -kx:kx); rr <- rc[1] + g$dr; cc <- rc[2] + g$dc
  ok <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc
  cell <- cellFromRowCol(tmpl, rr[ok], cc[ok]); xy <- xyFromCell(tmpl, cell)
  dk <- sqrt(((xy[, 1] - x) * cos(y * pi / 180))^2 + (xy[, 2] - y)^2) * 111.32
  cell[dk <= km & !is.na(ev[cell]) & ev[cell] <= thr & sv[cell] %in% sd]
}
fcells <- lapply(seq_len(nrow(u)), function(j) {
  c5 <- near_cells(u$lon_flank[j], u$lat_flank[j], P$flank_buffer_km, u$flank_m[j], u$side[j])
  if (length(c5)) c5 else cellFromXY(tmpl, cbind(u$lon_flank[j], u$lat_flank[j]))
})
ccell <- cellFromXY(tmpl, cbind(u$lon, u$lat))
uc <- unique(c(ccell, unlist(fcells))); key <- function(x) match(x, uc)
xy <- xyFromCell(tmpl, uc)
B <- as.matrix(extract(base, uc))
message(nrow(u), " col-sides, ", length(uc), " cells")

# present lowland envelope: cells <= threshold, same side, within 100 km of the flank point (the M2 Mahalanobis pool)
env <- rast(file.path(wd, paste0("WC_bio", c(5, 6), ".tif"))); names(env) <- c("b5", "b6")
rng <- t(vapply(seq_len(nrow(u)), function(j) {
  x <- extract(env, near_cells(u$lon_flank[j], u$lat_flank[j], 100, u$flank_m[j], u$side[j]))
  if (!nrow(x)) return(rep(NA_real_, 4))
  c(apply(x, 2, min, na.rm = TRUE), apply(x, 2, max, na.rm = TRUE))
}, numeric(4)))
colnames(rng) <- paste0(rep(c("lo_", "hi_"), each = 2), c("b5", "b6"))

# ---- metrics for one scenario (dTx, dTn, rP: cells x 12) ----
ovl <- function(lo1, hi1, lo2, hi2) (pmin(hi1, hi2) - pmax(lo1, lo2)) / sqrt((hi1 - lo1) * (hi2 - lo2))
cellclim <- function(dTx = 0, dTn = 0, rP = 1) {
  Tx <- B[, tx] + dTx; Tn <- B[, tn] + dTn; Pm <- B[, pr] * rP
  cbind(b5 = apply(Tx, 1, max), b6 = apply(Tn, 1, min), P = rowSums(Pm))
}
agg <- function(C) {
  cr <- C[key(ccell), , drop = FALSE]
  fl <- t(vapply(fcells, function(x) colMeans(C[key(x), , drop = FALSE], na.rm = TRUE), numeric(3)))
  list(cr = cr, fl = fl)
}
A0 <- agg(cellclim())
stopifnot(cor(A0$cr[, "b5"], u$crest_bio5, use = "complete") > 0.99,           # recomputed bio5/bio6 match M2
          cor(A0$fl[, "b6"], u$flank_bio6, use = "complete") > 0.99)
# dT* = uniform warming that centres the crest range [bio6, bio5] on the present flank range, i.e. maximises TO_fix
mid <- function(A) (A[, "b5"] + A[, "b6"]) / 2
dTstar <- mid(A0$fl) - mid(A0$cr)
to_at <- function(s) ovl(A0$cr[, "b6"] + s, A0$cr[, "b5"] + s, A0$fl[, "b6"], A0$fl[, "b5"])
stopifnot(all(to_at(dTstar) >= pmax(to_at(dTstar - 0.5), to_at(dTstar + 0.5)) - 1e-9, na.rm = TRUE))
metrics <- function(C) {
  A <- agg(C)
  data.frame(u[, c("col_id", "side_lab")],
             TO_fix = ovl(A$cr[, "b6"], A$cr[, "b5"], A0$fl[, "b6"], A0$fl[, "b5"]),
             Pg_fix = log(A$cr[, "P"] / A0$fl[, "P"]),
             # >= 1: crest centred on today's lowland; ponytail: NA when dT* < 0.5 C (crest already centred, ratio explodes)
             dTstar = dTstar, frac_dTstar = ifelse(dTstar < 0.5, NA, (mid(A$cr) - mid(A0$cr)) / dTstar),
             over_T = A$cr[, "b5"] > rng[, "hi_b5"] | A$cr[, "b6"] > rng[, "hi_b6"],  # hotter than any present lowland cell
             crest_b5 = A$cr[, "b5"], crest_b6 = A$cr[, "b6"], crest_P = A$cr[, "P"],
             flank_b5 = A$fl[, "b5"], flank_b6 = A$fl[, "b6"], flank_P = A$fl[, "P"])
}

# coarse delta raster -> cells x 12 (coastal NA filled outward, then bilinear)
fill <- function(r, n = 4) { for (i in 1:n) r <- cover(r, focal(r, 3, mean, na.policy = "only", na.rm = TRUE)); r }
at_cells <- function(r) {                                     # extract() on a coordinate matrix returns no ID column
  m <- as.matrix(extract(fill(r), xy, method = "bilinear")); stopifnot(dim(m) == c(nrow(xy), 12)); m
}

# ---- future ----
GCM <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-ESM2", "CNRM-ESM2-1", "GISS-E2-1-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")
SSP <- c("126", "370", "585"); PER <- c("2041-2060", "2061-2080", "2081-2100")
w10 <- lapply(c(tmax = "tmax", tmin = "tmin", prec = "prec"), function(v) crop(geodata::worldclim_global(v, 10, path = fd), box))
fut <- do.call(rbind, lapply(GCM, function(g) do.call(rbind, lapply(SSP, function(s) do.call(rbind, lapply(PER, function(p) {
  f <- lapply(c(tmax = "tmax", tmin = "tmin", prec = "prec"), function(v) crop(geodata::cmip6_world(g, s, p, v, 10, path = fd), box))
  C <- cellclim(at_cells(f$tmax - w10$tmax), at_cells(f$tmin - w10$tmin), at_cells((f$prec + 1) / (w10$prec + 1)))
  data.frame(gcm = g, ssp = s, period = p, metrics(C))
}))))))

# ---- worst side per col and write ----
worst <- function(d, by) {
  mn <- function(x) if (all(is.na(x))) NA else min(x, na.rm = TRUE)                # frac is NA on an already-centred side
  a <- aggregate(d[c("TO_fix", "frac_dTstar")], d[, c("col_id", by)], mn)
  b <- merge(aggregate(d["Pg_fix"], d[, c("col_id", by)], function(x) max(abs(x))), aggregate(d["dTstar"], d[, c("col_id", by)], max))
  f <- aggregate(d["over_T"], d[, c("col_id", by)], any)
  k <- unique(u[, c("col_id", "name", "lon", "lat", "elev_cross", "d_km", "tropical", "gateway", "huanca", "lara")])
  merge(k, merge(merge(a, b), f))
}
fw <- worst(fut, c("gcm", "ssp", "period"))
by <- fw[, c("col_id", "ssp", "period")]
ens <- do.call(data.frame, aggregate(fw[, c("TO_fix", "Pg_fix", "frac_dTstar")], by,
                                     function(x) c(mean = mean(x), min = min(x), max = max(x))))
ens <- merge(ens, aggregate(cbind(p_over_T = fw$over_T), by, mean))   # share of GCMs flagged
# H4 change vs present (worst side, no delta) and GCM sign agreement = share of GCMs with the sign of the ensemble mean;
# robust if >= 0.8 (7 of 8; AR6 simple approach, BITACORA 2026-11-03)
now <- worst(data.frame(gcm = "now", ssp = "now", period = "now", metrics(cellclim())), c("gcm", "ssp", "period"))
fw$dTO <- fw$TO_fix - now$TO_fix[match(fw$col_id, now$col_id)]
fw$dPg <- fw$Pg_fix - now$Pg_fix[match(fw$col_id, now$col_id)]
agree <- function(x) mean(sign(x) == sign(mean(x)))
ens <- merge(ens, do.call(data.frame, aggregate(fw[, c("dTO", "dPg")], by,
                                                function(x) c(mean = mean(x), min = min(x), max = max(x), agree = agree(x)))))
ens$TO_now <- now$TO_fix[match(ens$col_id, now$col_id)]; ens$Pg_now <- now$Pg_fix[match(ens$col_id, now$col_id)]
ens$gcm <- "ensemble"
ds <- data.frame(u[, c("col_id", "name", "side_lab", "lat", "elev_cross", "gateway", "tropical", "huanca")], dTstar = dTstar,
                 TO_now = to_at(0), TO_max = to_at(dTstar), rng)
write.csv(fw, outs[1], row.names = FALSE)
write.csv(ens, outs[2], row.names = FALSE)
write.csv(ds, outs[3], row.names = FALSE)
write.csv(fut, file.path(int, "M3_futuro_por_lado.csv"), row.names = FALSE)
por <- u$col_id[u$name %in% "Abra de Porculla"][1]
print(subset(ens, col_id == por))
