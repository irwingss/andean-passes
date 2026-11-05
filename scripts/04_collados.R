# 04_collados.R - M1.4: cols on the crest profile with depth, isolation (N/S), width, gateway flags
# and montane gaps (plan section 4). Profile is cut at breaks (gaps >= 2 km) only for width;
# depth and isolation run along the whole N->S profile, as defined.
# Rscript scripts/04_collados.R [--force]
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; P <- cfg$params
force <- "--force" %in% commandArgs(TRUE)
out <- file.path(int, "M1_cols_raw.csv")
if (file.exists(out) && !force) { message("exists, skip"); quit(save = "no") }

prof <- read.csv(file.path(int, "M1_profile.csv"))
prof <- prof[!is.na(prof$elev_cross), ]
e <- prof$elev_cross; d <- prof$d_km; n <- length(e)
# summits between cols: highest crossing height per bin when the profile has it (GLO-30 profile)
pk <- if (is.null(prof$elev_peak)) e else prof$elev_peak

# Local minima: runs of equal values lower than both neighbouring runs (ends excluded)
r <- rle(e); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
nr <- length(r$values)
ismin <- c(FALSE, r$values[-c(1, nr)] < r$values[-c(nr - 1, nr)] & r$values[-c(1, nr)] < r$values[-c(1, 2)], FALSE)
ci <- (starts[ismin] + ends[ismin]) %/% 2                  # centre of the flat run

# Nearest strictly lower point to the N (left) and S (right): monotonic stack, O(n)
nearest_lower <- function(v) {
  res <- rep(NA_integer_, length(v)); st <- integer(0)
  for (i in seq_along(v)) {
    while (length(st) && v[st[length(st)]] >= v[i]) st <- st[-length(st)]
    if (length(st)) res[i] <- st[length(st)]
    st <- c(st, i)
  }
  res
}
Ln <- nearest_lower(e); Rs <- n + 1L - rev(nearest_lower(rev(e)))
Rs[is.na(Rs)] <- NA

cols <- do.call(rbind, lapply(ci, function(i) {
  l <- Ln[i]; s <- Rs[i]
  maxN <- max(pk[(if (is.na(l)) 1 else l + 1):i]); maxS <- max(pk[i:(if (is.na(s)) n else s - 1)])
  # width: contiguous profile below elev + 100 m around the col
  hi <- e > e[i] + 100
  wl <- suppressWarnings(max(which(hi[1:i]))); wr <- suppressWarnings(min(which(hi[i:n]))) + i - 1
  wl <- if (is.finite(wl)) d[wl] else d[1]; wr <- if (is.finite(wr)) d[wr] else d[n]
  data.frame(bin = prof$bin[i], d_km = d[i], lon = prof$lon[i], lat = prof$lat[i],
             elev = prof$elev[i], elev_cross = e[i], borders_I = prof$borders_I[i],
             depth = min(maxN, maxS) - e[i],
             iso_N_km = if (is.na(l)) NA else d[i] - d[l],
             iso_S_km = if (is.na(s)) NA else d[s] - d[i],
             iso_N_toend = is.na(l), iso_S_toend = is.na(s),
             width_km = wr - wl)
}))
cols$iso_km <- pmin(cols$iso_N_km, cols$iso_S_km, na.rm = TRUE)
cols$iso_km[is.na(cols$iso_N_km) & is.na(cols$iso_S_km)] <- Inf   # global minimum
cols$tropical <- cols$lat <= P$tropics_lat[1] & cols$lat >= P$tropics_lat[2]

# Gateways: depth >= threshold and lowest point of the profile within +-window km
gate <- function(dep, win) sapply(seq_len(nrow(cols)), function(j) {
  if (cols$depth[j] < dep) return(FALSE)
  w <- abs(d - cols$d_km[j]) <= win
  cols$elev_cross[j] <= min(e[w])
})
for (win in P$gateway_window_km) cols[[paste0("gateway_", win)]] <- gate(P$depth_m[1], win)
cols$gateway <- cols[[paste0("gateway_", P$gateway_window_km[1])]]
write.csv(cols, out, row.names = FALSE)

# Montane gaps: profile runs below each threshold
mg <- do.call(rbind, lapply(P$montane_gap_m, function(th) {
  rr <- rle(e < th); en <- cumsum(rr$lengths); st <- en - rr$lengths + 1
  k <- which(rr$values)
  if (!length(k)) return(NULL)
  data.frame(threshold_m = th, lat_start = prof$lat[st[k]], lat_end = prof$lat[en[k]],
             length_km = d[en[k]] - d[st[k]] + cfg$params$segment_km,
             min_elev = sapply(k, function(z) min(e[st[z]:en[z]])))
}))
write.csv(mg, file.path(int, "M1_montane_gaps.csv"), row.names = FALSE)
cat("cols:", nrow(cols), " depth>=100:", sum(cols$depth >= 100), " gateways:", sum(cols$gateway), "\n")
