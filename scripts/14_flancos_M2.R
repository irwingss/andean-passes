# 14_flancos_M2.R - M2.1: lowland flanks of each col (plan section 4).
# Cols: depth >= 100 m everywhere, plus the Huancabamba system (depth >= 50 m, < 3000 m, 4-7.5 S).
# For each side (W; E = E + endorheic I connected to E) the start is the nearest HydroSHEDS 15 s cell of
# that domain within 5 km of the col; the steepest-descent trace (HydroSHEDS D8) runs downstream until it
# first reaches each flank threshold (500/1000/1500 m) or 150 km of path. If it does not, the flank is the
# nearest cell of that side below the threshold (flagged fallback). Traces of all cols advance together.
# Rscript scripts/14_flancos_M2.R [--force]
suppressPackageStartupMessages({ library(sf); library(terra) })
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; raw <- cfg$paths$data_dir; P <- cfg$params
force <- "--force" %in% commandArgs(TRUE)
out <- file.path(int, "M2_flancos.csv")
if (file.exists(out) && !force) { message("exists, skip"); quit(save = "no") }

cols <- read.csv(file.path(cfg$paths$outputs, "database", "APD_cols.csv"))
H <- P$huancabamba
cols$huanca <- cols$lat >= H$lat[1] & cols$lat <= H$lat[2] & cols$elev_cross < H$max_elev_m & cols$depth >= H$min_depth_m
cols <- cols[cols$depth >= P$depth_m[1] | cols$huanca, ]
message(nrow(cols), " cols")

# Side raster on the 15 s grid: 1 = W (incl. coastal endorheic sinks), 2 = E (incl. I connected to E)
dom <- vect(file.path(int, "M1_domains.gpkg"))
tmpl <- rast(file.path(int, "M1_domains_15s.tif"))
E <- dom[dom$dom == "E"]; I <- disagg(dom[dom$dom == "I"]); conn <- relate(I, E, "intersects")[, 1]
side <- rasterize(dom[dom$dom == "W"], tmpl, field = 1)
side <- cover(side, rasterize(E, tmpl, field = 2))
side <- cover(side, rasterize(I[!conn], tmpl, field = 1))
# Fallback targets exclude I connected to E: that polygon runs from the Altiplano down to the Pacific-slope
# sinks (Pampa del Tamarugal, Atacama), so an "E" fallback could land on the Pacific coast
fbv <- values(side)[, 1]
side <- cover(side, rasterize(I[conn], tmpl, field = 2))
dir <- crop(rast(file.path(raw, "hyd_sa_dir_15s", "hyd_sa_dir_15s.tif")), tmpl)
dem <- crop(rast(file.path(raw, "hyd_sa_dem_15s", "hyd_sa_dem_15s.tif")), tmpl)
nc <- ncol(tmpl); nr <- nrow(tmpl)
sv <- values(side)[, 1]; dv <- values(dir)[, 1]; ev <- values(dem)[, 1]

# Start cells: nearest cell of each side within 5 km of the col
win <- 5 / 111.32 / res(tmpl)[1]                     # 5 km in cells (lat direction), ~12
starts <- do.call(rbind, lapply(seq_len(nrow(cols)), function(j) {
  c0 <- cellFromXY(tmpl, cbind(cols$lon[j], cols$lat[j]))
  rc <- rowColFromCell(tmpl, c0); k <- ceiling(win)
  g <- expand.grid(dr = -k:k, dc = -k:k)
  rr <- rc[1] + g$dr; cc <- rc[2] + g$dc; ok <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc
  cell <- cellFromRowCol(tmpl, rr[ok], cc[ok]); xy <- xyFromCell(tmpl, cell)
  dk <- sqrt(((xy[, 1] - cols$lon[j]) * cos(cols$lat[j] * pi / 180))^2 + (xy[, 2] - cols$lat[j])^2) * 111.32
  s <- sv[cell]
  do.call(rbind, lapply(1:2, function(sd) {
    i <- which(s == sd & dk <= 5)
    if (!length(i)) return(data.frame(col_id = cols$col_id[j], side = sd, cell = NA, start_km = NA))
    i <- i[which.min(dk[i])]
    data.frame(col_id = cols$col_id[j], side = sd, cell = cell[i], start_km = dk[i])
  }))
}))
message(sum(is.na(starts$cell)), " col-sides without a start cell within 5 km")

# D8 offsets (HydroSHEDS: 1 E, 2 SE, 4 S, 8 SW, 16 W, 32 NW, 64 N, 128 NE)
code <- c(1, 2, 4, 8, 16, 32, 64, 128)
dr <- c(0, 1, 1, 1, 0, -1, -1, -1); dc <- c(1, 1, 0, -1, -1, -1, 0, 1)
off <- setNames(dr * nc + dc, code); diag <- setNames(dr != 0 & dc != 0, code)
thr <- sort(P$flank_m); maxkm <- P$flank_max_trace_km
cs <- res(tmpl)[1] * 111.32                          # cell size in km (lat direction)

n <- nrow(starts); cur <- starts$cell; dist <- rep(0, n); live <- !is.na(cur)
hit_cell <- matrix(NA_real_, n, length(thr)); hit_km <- hit_cell
hist <- list(cur)                                    # cell of every trace at every step
while (any(live)) {
  i <- which(live)
  for (t in seq_along(thr)) {
    h <- i[is.na(hit_cell[i, t]) & !is.na(ev[cur[i]]) & ev[cur[i]] <= thr[t]]
    hit_cell[h, t] <- cur[h]; hit_km[h, t] <- dist[h]
  }
  d <- as.character(dv[cur[i]])
  stop_ <- !(d %in% names(off)) | !is.na(hit_cell[i, 1]) | dist[i] > maxkm
  live[i[stop_]] <- FALSE; i <- i[!stop_]; d <- d[!stop_]
  if (!length(i)) break
  lat <- yFromCell(tmpl, cur[i])
  step <- ifelse(diag[d], sqrt(cs^2 + (cs * cos(lat * pi / 180))^2), ifelse(dr[match(d, code)] == 0, cs * cos(lat * pi / 180), cs))
  cur[i] <- cur[i] + off[d]; dist[i] <- dist[i] + step
  nxt <- hist[[length(hist)]]; nxt[i] <- cur[i]; nxt[-i] <- NA; hist[[length(hist) + 1]] <- nxt
}

# Fallback: nearest exorheic cell of that side below the threshold (ponytail: searched within +-3 deg (Altiplano: coast > 2 deg from the divide))
fall <- function(j, t) {
  c0 <- starts$cell[j]; if (is.na(c0)) {
    cl <- cols[cols$col_id == starts$col_id[j], ]; c0 <- cellFromXY(tmpl, cbind(cl$lon, cl$lat))
  }
  xy0 <- xyFromCell(tmpl, c0)
  for (r in if (starts$side[j] == 2) c(3, 6) else 3) {  # E retry at +-6 deg: southern Altiplano to the Chaco
    cell <- cells(tmpl, ext(xy0[1] - r, xy0[1] + r, xy0[2] - r, xy0[2] + r))
    cell <- cell[fbv[cell] %in% starts$side[j] & !is.na(ev[cell]) & ev[cell] <= thr[t]]
    if (length(cell)) break
  }
  if (!length(cell)) return(c(NA, NA))
  xy <- xyFromCell(tmpl, cell)
  dk <- sqrt(((xy[, 1] - xy0[1]) * cos(xy0[2] * pi / 180))^2 + (xy[, 2] - xy0[2])^2) * 111.32
  c(cell[which.min(dk)], min(dk))
}
res_ <- do.call(rbind, lapply(seq_along(thr), function(t) {
  fb <- is.na(hit_cell[, t])
  for (j in which(fb)) { f <- fall(j, t); hit_cell[j, t] <- f[1]; hit_km[j, t] <- f[2] }
  xy <- xyFromCell(tmpl, hit_cell[, t])
  data.frame(starts[, c("col_id", "side", "start_km")], side_lab = c("W", "E")[starts$side], flank_m = thr[t],
             lon = xy[, 1], lat = xy[, 2], route_km = hit_km[, t], fallback = fb, elev15 = ev[hit_cell[, t]])
}))
write.csv(res_, out, row.names = FALSE)

# Routes (default threshold) for maps and the dry-forest gap
t1 <- match(P$flank_m[1], thr)
hm <- do.call(cbind, hist)
geo <- lapply(seq_len(n), function(j) {
  p <- hm[j, ]; p <- p[!is.na(p)]; if (!length(p) || is.na(hit_cell[j, t1]) || res_$fallback[(t1 - 1) * n + j]) return(NULL)
  p <- p[seq_len(match(hit_cell[j, t1], p))]
  if (length(p) < 2) p <- c(p, p)
  st_linestring(xyFromCell(tmpl, p))
})
keep <- !vapply(geo, is.null, TRUE)
routes <- st_sf(starts[keep, c("col_id", "side")], side_lab = c("W", "E")[starts$side[keep]],
                geometry = st_sfc(geo[keep], crs = 4326))
st_write(routes, file.path(int, "M2_flank_routes.gpkg"), delete_dsn = TRUE, quiet = TRUE)

s <- res_[res_$flank_m == P$flank_m[1], ]
message(sprintf("flank %d m: %d col-sides, %d fallback, median route %.1f km",
                P$flank_m[1], nrow(s), sum(s$fallback), median(s$route_km, na.rm = TRUE)))
print(s[s$col_id %in% cols$col_id[cols$name %in% "Abra de Porculla"], ])
