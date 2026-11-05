# 29_corredor_montano_M3.R - M3.5: climatic montane gaps along the divide (N-S montane corridor), plan v2.2.
# Lower limit of montane climate along the tropical divide: 2000 m (Josse et al. 2009); 1500 and 2500 m as sensitivity.
#   present gap = crossing height < threshold (2000 m; 1500 / 2500 m sensitivity), i.e. the M1 montane gaps;
#   future gap  = crossing height < threshold + dT / G, where dT = annual mean warming of the bin ((dTmax + dTmin) / 2 over
#   12 months, same 10' WorldClim CMIP6 delta as 20_permeabilidad_M3.R) and G = local lapse rate, slope of WorldClim bio1
#   on the HydroSHEDS 15 s DEM (30" mean) over cells <= 50 km from the divide, 1-degree latitude windows (< 200 cells:
#   nearest window). The isotherm of today's threshold moves up by dT / |G|.
# Metrics (tropical Andes, divide from the N end to the last bin at >= 23.5 S; and the whole chain): km of divide in gap,
# gaps >= 10 km, longest gap, montane fragments >= 10 km, longest montane stretch. Per bin: share of GCMs in gap;
# robust loss = montane today and in gap in >= 7 of 8 GCMs.
# Rscript scripts/29_corredor_montano_M3.R
suppressPackageStartupMessages(library(terra))
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; raw <- cfg$paths$data_dir; P <- cfg$params
tab <- file.path(cfg$paths$outputs, "tables"); fd <- file.path(raw, "cmip6"); box <- ext(-82, -60, -56, 13)
TH <- unlist(P$montane_gap_m); seg <- P$segment_km; MINK <- 10

prof <- read.csv(file.path(int, "M1_profile.csv"))
prof <- prof[order(prof$bin), ]
trop_end <- max(prof$d_km[prof$lat >= -23.5])

# ---- local lapse rate G (deg C per m) per 1-degree window ----
bio1 <- rast(file.path(raw, "clima_wc", "WC_bio1.tif"))
dem <- resample(crop(rast(file.path(raw, "hyd_sa_dem_15s", "hyd_sa_dem_15s.tif")), bio1), bio1, "average")
pv <- vect(prof[seq(1, nrow(prof), 10), c("lon", "lat")], geom = c("lon", "lat"), crs = "EPSG:4326")
near <- rasterize(aggregate(buffer(pv, 50000)), bio1)
v <- na.omit(data.frame(lat = crds(bio1, na.rm = FALSE)[, 2], t = values(bio1)[, 1], z = values(dem)[, 1], m = values(near)[, 1]))
v$w <- floor(v$lat)
lap <- do.call(rbind, lapply(split(v, v$w), function(d) data.frame(w = d$w[1], n = nrow(d),
                                                                   G = if (nrow(d) >= 200) coef(lm(t ~ z, d))[2] else NA)))
ok <- which(!is.na(lap$G)); lap$G <- lap$G[ok[vapply(seq_len(nrow(lap)), function(i) which.min(abs(ok - i)), 1L)]]
stopifnot(all(lap$G < 0))
prof$G <- lap$G[match(floor(prof$lat), lap$w)]
message("lapse rate (C/km): ", paste(round(range(-1000 * prof$G), 2), collapse = " to "))

# ---- annual mean warming per bin and scenario ----
fill <- function(r, n = 4) { for (i in 1:n) r <- cover(r, focal(r, 3, mean, na.policy = "only", na.rm = TRUE)); r }
w10 <- lapply(c(tmax = "tmax", tmin = "tmin"), function(x) crop(geodata::worldclim_global(x, 10, path = fd), box))
xy <- as.matrix(prof[, c("lon", "lat")])
GCM <- c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-ESM2", "CNRM-ESM2-1", "GISS-E2-1-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")
SSP <- c("126", "370", "585"); PER <- c("2041-2060", "2061-2080", "2081-2100")
sc <- expand.grid(gcm = GCM, ssp = SSP, period = PER, stringsAsFactors = FALSE)
dT <- vapply(seq_len(nrow(sc)), function(i) {
  f <- lapply(c(tmax = "tmax", tmin = "tmin"), function(x) crop(geodata::cmip6_world(sc$gcm[i], sc$ssp[i], sc$period[i], x, 10, path = fd), box))
  d <- (mean(f$tmax - w10$tmax) + mean(f$tmin - w10$tmin)) / 2
  extract(fill(d), xy, method = "bilinear")[, 1]
}, numeric(nrow(prof)))
stopifnot(!anyNA(dT))

# ---- gaps and fragments ----
runs <- function(g, d) {                                         # g: logical gap per bin, d: d_km
  r <- rle(g); en <- cumsum(r$lengths); st <- en - r$lengths + 1; len <- d[en] - d[st] + seg
  c(gap_km = sum(len[r$values]), n_gap10 = sum(r$values & len >= MINK), max_gap_km = max(c(0, len[r$values])),
    n_frag10 = sum(!r$values & len >= MINK), max_frag_km = max(c(0, len[!r$values])))
}
zones <- list(tropical = prof$d_km <= trop_end, all = rep(TRUE, nrow(prof)))
shift <- dT / -prof$G                                            # isotherm rise (m), bins x scenarios
met <- do.call(rbind, lapply(TH, function(th) do.call(rbind, lapply(names(zones), function(z) {
  k <- zones[[z]]; d <- prof$d_km[k]; e <- prof$elev_cross[k]
  now <- data.frame(threshold_m = th, zone = z, gcm = "present", ssp = "present", period = "1970-2000", t(runs(e < th, d)))
  fut <- data.frame(threshold_m = th, zone = z, sc, t(vapply(seq_len(nrow(sc)), function(i) runs(e < th + shift[k, i], d), numeric(5))))
  rbind(now, fut)
}))))
fm <- met[met$gcm != "present", ]
ens <- do.call(data.frame, aggregate(fm[, c("gap_km", "n_gap10", "max_gap_km", "n_frag10", "max_frag_km")],
                                     fm[, c("threshold_m", "zone", "ssp", "period")], function(x) c(mean = mean(x), min = min(x), max = max(x))))

# per bin (default threshold): share of GCMs in gap per SSP x period, robust loss
th <- TH[1]; key <- paste(sc$ssp, sc$period, sep = "_")
fr <- sapply(unique(key), function(k) rowMeans(prof$elev_cross < th + shift[, key == k, drop = FALSE]))
colnames(fr) <- paste0("pgap_", unique(key))
bins <- data.frame(prof[, c("bin", "d_km", "lon", "lat", "elev_cross")], tropical = zones$tropical, G_C_per_km = -1000 * prof$G,
                   gap_now = prof$elev_cross < th,
                   rise_m_370_2090 = rowMeans(shift[, key == "370_2081-2100"]), fr, check.names = FALSE)
bins$loss_robust_370_2090 <- !bins$gap_now & bins$`pgap_370_2081-2100` >= 0.875
# robust loss at every threshold (sensitivity of the corridor and of the critical gateways, 30_proteccion_M3.R)
k370 <- key == "370_2081-2100"
lt <- sapply(TH, function(t) prof$elev_cross >= t & rowMeans(prof$elev_cross < t + shift[, k370, drop = FALSE]) >= 0.875)
colnames(lt) <- paste0("loss_", TH)
write.csv(data.frame(bin = prof$bin, lt), file.path(int, "M3_perdida_por_umbral.csv"), row.names = FALSE)
write.csv(met, file.path(tab, "M3_corredor_montano.csv"), row.names = FALSE)
write.csv(ens, file.path(tab, "M3_corredor_montano_ensamble.csv"), row.names = FALSE)
write.csv(bins, file.path(int, "M3_corredor_bins.csv"), row.names = FALSE)
write.csv(lap, file.path(int, "M3_lapse_rate.csv"), row.names = FALSE)
print(subset(met, gcm == "present" & zone == "tropical"))
print(subset(ens, zone == "tropical" & threshold_m == th), digits = 3)
