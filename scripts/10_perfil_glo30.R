# 10_perfil_glo30.R - M1: ordered N->S divide chain and 0.5 km crossing-height profile from the
# GLO-30 flooding (09_divisoria_glo30.py). Writes M1_profile.csv / M1_divide_bins.gpkg for 04_collados.R.
# Per bin: elev_cross = lowest crossing height (cols), elev_peak = highest (summits for depth),
# elev = topographic crest at the lowest crossing; lon/lat = that vertex.
# Rscript scripts/10_perfil_glo30.R
suppressPackageStartupMessages({ library(sf); library(terra) })
sf_use_s2(FALSE)
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; seg <- cfg$params$segment_km

fs <- list.files(file.path(int, "M1_glo30_divide"), "^blk_.*\\.csv$", full.names = TRUE)
v <- do.call(rbind, lapply(fs, function(f) { x <- read.csv(f); if (nrow(x)) x$blk <- basename(f); x }))
# Clip to GMBA Andes + 20 km (the corridor reaches 100 km out to seed the Coastal Cordillera; the divide
# beyond the Andes, e.g. the Venezuelan Coastal Range, is not part of the census). A piece cut by the mask
# splits into separate runs.
aea <- "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +datum=WGS84 +units=m"
m20 <- st_transform(st_buffer(st_transform(read_sf(file.path(int, "M1_gmba_andes.gpkg")), aea), cfg$params$gmba_margin_km * 1000), 4326)
mr <- rasterize(vect(m20), rast(file.path(int, "M1_markers_15s.tif")), field = 1, background = 0)
v <- v[order(v$blk, v$contour, v$run, v$seq), ]
v$inm <- mr[cellFromXY(mr, cbind(v$lon, v$lat))][, 1] == 1
v$run <- paste(v$run, cumsum(c(TRUE, diff(v$inm) != 0)))
v <- v[v$inm, ]
v$piece <- as.integer(factor(paste(v$blk, v$contour, v$run)))
v <- v[order(v$piece, v$seq), ]
dk <- function(m) c(0, geosphere::distGeo(m[-nrow(m), ], m[-1, ]) / 1000)
pl <- split(v[, c("lon", "lat", "cross", "crest")], v$piece)
len <- sapply(pl, function(p) sum(dk(as.matrix(p[, 1:2]))))
cat("pieces:", length(pl), " km total:", round(sum(len)), " pieces >= 2 km:", sum(len >= 2), "km:", round(sum(len[len >= 2])), "\n")
# ponytail: pieces < 2 km dropped (salar / flat-floor fragments); a real divide piece only gets that short
# where a block edge clips it, and the < 2 km gap closure below rejoins those.
pl <- pl[len >= 2]

# --- Greedy N->S chaining (as in 03_perfil_collados.R) ---
ends <- t(sapply(pl, function(p) c(p$lon[1], p$lat[1], p$lon[nrow(p)], p$lat[nrow(p)])))
start <- which.max(pmax(ends[, 2], ends[, 4]))
cur <- pl[[start]]; if (cur$lat[1] < cur$lat[nrow(cur)]) cur <- cur[nrow(cur):1, ]
chain <- list(cur); gaps <- 0; left <- setdiff(seq_along(pl), start); dropped <- 0; edge_jump <- FALSE
while (length(left)) {
  e <- c(cur$lon[nrow(cur)], cur$lat[nrow(cur)])
  d1 <- geosphere::distGeo(e, ends[left, 1:2, drop = FALSE]); d2 <- geosphere::distGeo(e, ends[left, 3:4, drop = FALSE])
  # prefer pieces that reach up to the current latitude (same or next block), so a long jump across the
  # plateau to a piece further south never wins over the piece that continues the line
  ok <- pmax(ends[left, 2], ends[left, 4]) >= e[2] - 0.05
  if (!any(ok)) ok[] <- TRUE
  j <- which(ok)[which.min(pmin(d1, d2)[ok])]; k <- left[j]
  # a piece > 500 km from the chain end is a fragment the greedy pass skipped: drop it, do not jump back
  if (min(d1[j], d2[j]) > 500e3) { dropped <- dropped + len[names(pl)[k]]; left <- left[-j]; next }
  y0 <- e[2]
  # orientation: N->S for pieces spanning > 0.2 deg of latitude, otherwise the end nearest the chain
  rev_k <- if (abs(ends[k, 2] - ends[k, 4]) > 0.2) ends[k, 2] < ends[k, 4] else d1[j] > d2[j]
  cur <- if (!rev_k) pl[[k]] else pl[[k]][nrow(pl[[k]]):1, ]
  dj <- if (!rev_k) d1[j] else d2[j]
  # Jump that starts and ends on the same block edge: across the Altiplano-Puna plateau the W|E contact is
  # not unique (the whole plateau sits on the crest; where the two floods touch depends on the block), so
  # the line continues there with no lateral distance and no break. Crossing heights are unaffected.
  edge <- abs(y0 - round(y0)) < 0.001 && abs(cur$lat[1] - round(y0)) < 0.001
  chain[[length(chain) + 1]] <- cur; gaps <- c(gaps, dj / 1000); left <- left[-j]
  edge_jump <- c(edge_jump, edge)
}
geff <- ifelse(edge_jump, 0, gaps)
pid <- rep(seq_along(chain), sapply(chain, nrow))
xy <- do.call(rbind, chain)
# Chainage measured with a 0.5 km ruler (the bin size): 30 m vertices trace every zigzag of the crest, so a
# vertex-by-vertex length would grow with DEM resolution. Each vertex takes the chainage of the last ruler mark.
ruler <- function(lon, lat) {
  n <- length(lon); d <- numeric(n); a <- 1; da <- 0; kx <- 111.32 * cos(lat * pi / 180)
  for (i in seq_len(n)[-1]) {
    s <- sqrt(((lon[i] - lon[a]) * kx[i])^2 + ((lat[i] - lat[a]) * 110.57)^2)
    if (s >= seg) { da <- da + s; a <- i }
    d[i] <- da
  }
  d
}
dl <- lapply(chain, function(p) ruler(p$lon, p$lat))
off <- cumsum(c(0, sapply(dl, max)[-length(dl)]) + geff)   # pieces end to end, real gaps count as distance
newp <- c(TRUE, diff(pid) != 0)
brk <- newp & geff[pid] >= 2; brk[1] <- TRUE
xy$d <- unlist(dl) + off[pid]; xy$grp <- cumsum(brk); xy$gapb <- ifelse(brk, geff[pid], 0)

# --- 0.5 km bins within continuous stretches ---
xy$bin0 <- floor((xy$d - ave(xy$d, xy$grp, FUN = min)) / seg)
xy$key <- paste(xy$grp, xy$bin0)
xy <- xy[!is.na(xy$cross), ]
b <- split(seq_len(nrow(xy)), factor(xy$key, levels = unique(xy$key)))
prof <- do.call(rbind, lapply(b, function(i) {
  j <- i[which.min(xy$cross[i])]
  data.frame(d_km = mean(range(xy$d[i])), lon = xy$lon[j], lat = xy$lat[j], elev = xy$crest[j],
             elev_cross = xy$cross[j], elev_peak = max(xy$cross[i]), grp = xy$grp[j],
             gap_before_km = xy$gapb[i[1]])
}))
prof$bin <- seq_len(nrow(prof))
# W|I stretches: bin within 2 km of the endorheic basins connected to E
dom <- vect(file.path(int, "M1_domains.gpkg")); I <- disagg(dom[dom$dom == "I"])
I <- I[relate(I, dom[dom$dom == "E"], "intersects")[, 1]]
pts <- vect(prof[, c("lon", "lat")], geom = c("lon", "lat"), crs = "EPSG:4326")
prof$borders_I <- is.related(buffer(pts, 2000), I, "intersects")
prof <- prof[, c("bin", "d_km", "lon", "lat", "elev", "grp", "gap_before_km", "borders_I", "elev_cross", "elev_peak")]
write.csv(prof, file.path(int, "M1_profile.csv"), row.names = FALSE)

ln <- lapply(split(xy, xy$grp), function(z) as.matrix(z[, c("lon", "lat")]))
writeVector(vect(ln, type = "lines", crs = "EPSG:4326"), file.path(int, "M1_divide_glo30.gpkg"), overwrite = TRUE)
writeVector(cbind(pts, prof), file.path(int, "M1_divide_bins.gpkg"), overwrite = TRUE)
chs <- data.frame(lat_start = sapply(chain, function(p) p$lat[1]), lat_end = sapply(chain, function(p) p$lat[nrow(p)]),
                  lon_start = sapply(chain, function(p) p$lon[1]), km = sapply(dl, max), gap_before_km = gaps, edge_jump = edge_jump)
write.csv(chs, file.path(int, "M1_glo30_chain_pieces.csv"), row.names = FALSE)
cat("dropped fragments km:", round(dropped, 1), "
")
cat("bins:", nrow(prof), " length km:", round(max(xy$d)), " edge jumps:", sum(edge_jump & gaps >= 2), " breaks >= 2 km:", sum(geff >= 2),
    " largest gaps:", paste(round(sort(geff, decreasing = TRUE)[1:5], 1), collapse = ", "), "\n")
