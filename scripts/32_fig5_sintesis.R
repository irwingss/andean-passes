# 32_fig5_sintesis.R - Fig. 5 (synthesis), drawn on real relief: (A) Janzen at two gateways of similar height,
# (B) the moisture axis, (C) one warming mechanism with two effects (E-W opening, N-S montane break).
# Step 1 renders 3D blocks (rayshader, HydroSHEDS 15" DEM, WorldClim / Global-AI drapes) and projects the marked points
# to image pixels; step 2 writes the SVG (vectors on top of the renders) and exports PDF + 600 dpi PNG with Inkscape.
# Renders are cached in data/interim/fig5_renders (--force to redo).
# Rscript scripts/32_fig5_sintesis.R [--force]
suppressPackageStartupMessages({ library(terra); library(rayshader); library(data.table) })
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; raw <- cfg$paths$data_dir; tb <- file.path(cfg$paths$outputs, "tables")
P <- cfg$params; force <- "--force" %in% commandArgs(TRUE)
rd <- file.path(int, "fig5_renders"); dir.create(rd, showWarnings = FALSE)
dem0 <- rast(file.path(raw, "hyd_sa_dem_15s", "hyd_sa_dem_15s.tif"))
m2 <- fread(file.path(tb, "M2_metricas_clima_collados.csv"))[flank_m == 1000]
bt <- fread(file.path(int, "M3_corredor_bins.csv"))
thr <- P$montane_gap_m[[1]]

# ---- blocks: extent, aggregation, drape, view, points to project ----
pt <- function(id, sides = c("W", "E")) {
  x <- m2[col_id == id & side_lab %in% sides]
  rbind(x[1, .(key = "crest", lon, lat)], x[, .(key = paste0("flank_", side_lab), lon = lon_flank, lat = lat_flank)])
}
wc <- function(v) rast(file.path(raw, "clima_wc", paste0("WC_", v, ".tif")))
ai <- rast(file.path(raw, "clima", "GlobalAI_v31_annual.tif"))
belt <- function(lift) function(d) classify(d, rbind(c(-Inf, thr + lift, 1), c(thr + lift, Inf, 2)))   # 1 lowland, 2 montane
cbins <- bt[tropical == TRUE & lat < 0.5 & lat > -10]
rise_c <- median(bt[tropical == TRUE]$rise_m_370_2090)                       # same median as the text (whole tropical divide)
lost <- cbins[loss_robust_370_2090 == TRUE][order(-elev_cross)][1]                   # highest crest that loses montane climate
high <- cbins[which.max(elev_cross)]
B <- list(
  trop  = list(ext = c(-79.25, -77.75, -1.35, 0.45), fact = 1, ex = 2.2, theta = 0, phi = 22, drape = function(d) wc("bio1"),
               pts = pt("APD00868", "W")),
  temp  = list(ext = c(-70.85, -69.20, -34.65, -33.35), fact = 1, ex = 2.2, theta = 0, phi = 22, drape = function(d) wc("bio1"),
               pts = pt("APD03317", "W")),
  alti  = list(ext = c(-72.0, -64.8, -19.0, -15.7), fact = 2, ex = 4, theta = 0, phi = 32, drape = function(d) log10(ai),
               pts = pt("APD02335")),
  tach  = list(ext = c(-73.2, -71.6, 7.0, 8.4), fact = 1, ex = 3, theta = 0, phi = 32, drape = function(d) log10(ai),
               pts = pt("APD00203")),
  now   = list(ext = c(-81.0, -76.6, -10.0, 0.6), fact = 2, ex = 5, theta = 270, phi = 34, drape = belt(0),
               pts = rbind(lost[, .(key = "lost", lon, lat)], high[, .(key = "high", lon, lat)],
                     setnames(data.table(c("Colombia", "Ecuador", "Peru"), c(-77.0, -78.3, -77.0), c(0.35, -1.6, -7.6)), c("key", "lon", "lat"))), borders = TRUE),
  fut   = list(ext = c(-81.0, -76.6, -10.0, 0.6), fact = 2, ex = 5, theta = 270, phi = 34, drape = belt(rise_c),
               pts = rbind(lost[, .(key = "lost", lon, lat)], high[, .(key = "high", lon, lat)],
                     setnames(data.table(c("Colombia", "Ecuador", "Peru"), c(-77.0, -78.3, -77.0), c(0.35, -1.6, -7.6)), c("key", "lon", "lat"))), borders = TRUE))

# country borders (GADM 4.1, level 1 dissolved to level 0) for the panel-c blocks
borders_sf <- do.call(rbind, lapply(c("COL", "ECU", "PER"), function(i)
  aggregate(unwrap(readRDS(file.path(raw, "gadm", "gadm", sprintf("gadm41_%s_1_pk.rds", i))))))) |>
  sf::st_as_sf() |> sf::st_simplify(dTolerance = 0.005) |> sf::st_boundary() |> sf::st_geometry() |> sf::st_sf()

pal <- list(bio1 = list(lim = c(-8, 28), col = c("#2b5b84", "#8fb8d4", "#ebe6d6", "#e5aa63", "#b0492a")),
            ai   = list(lim = c(-2.7, 0.4), col = c("#8a5524", "#d8b779", "#eee6c8", "#93c189", "#2d7a48")),
            belt = list(lim = c(1, 2), col = c("#e9e2cf", "#3f8f55")))
pal_of <- function(k) if (k %in% c("trop", "temp")) "bio1" else if (k %in% c("alti", "tach")) "ai" else "belt"

render_block <- function(k) {
  b <- B[[k]]; png_f <- file.path(rd, paste0(k, ".png")); xy_f <- file.path(rd, paste0(k, "_pts.csv"))
  if (file.exists(png_f) && !force) return(invisible())
  d <- crop(dem0, ext(b$ext)); if (b$fact > 1) d <- aggregate(d, b$fact, mean, na.rm = TRUE)
  d[is.na(d)] <- 0; d <- clamp(d, lower = 0)                                           # sea at 0 m
  v <- resample(crop(b$drape(d), ext(b$ext), snap = "out"), d, "bilinear")
  m <- raster_to_matrix(d); vm <- raster_to_matrix(v)
  p <- pal[[pal_of(k)]]; vm <- pmin(pmax(vm, p$lim[1]), p$lim[2]); vm[is.na(vm)] <- p$lim[1]
  vm[1, 1] <- p$lim[1]; vm[nrow(vm), ncol(vm)] <- p$lim[2]                             # pin the colour scale to lim
  tex <- if (pal_of(k) == "belt") rep(p$col, each = 256) else colorRampPalette(p$col)(512)
  zs <- res(d)[2] * 111320 / b$ex
  img <- height_shade(vm, texture = tex) |>
    add_shadow(lamb_shade(m, zscale = zs / 3, sunangle = 300, sunaltitude = 40), 0.35) |>
    add_shadow(ambient_shade(m, zscale = zs), 0.2)
  if (isTRUE(b$borders)) img <- add_overlay(img, generate_line_overlay(borders_sf, extent = ext(d), heightmap = m, color = "grey15", linewidth = 1.2))
  plot_3d(img, m, zscale = zs, theta = b$theta, phi = b$phi, zoom = 0.62, fov = 25, windowsize = c(1600, 900),
          solid = TRUE, soliddepth = -thr / zs * 0.25, solidcolor = "#b9b4aa", solidlinecolor = "#b9b4aa", shadow = FALSE,
          background = "white")
  # project the marked points: let rayshader place them, read their 3D vertices, map to window fractions
  ex_ <- ext(d); e <- extract(d, as.matrix(b$pts[, .(lon, lat)]))[, 1]
  render_points(extent = ex_, lat = b$pts$lat, long = b$pts$lon, altitude = e, zscale = zs, heightmap = m, size = 1, color = "red")
  id <- max(rgl::ids3d()$id); w <- rgl::rgl.user2window(rgl::rgl.attrib(id, "vertices")); rgl::pop3d(id = id)
  render_snapshot(png_f, software_render = TRUE, width = 3600, height = 2025); rgl::close3d()   # high-res, same camera
  # trim white margins, keep pixel offsets for the projected points
  a <- png::readPNG(png_f); nw <- apply(a[, , 1:3], c(1, 2), min) < 0.985
  r <- range(which(rowSums(nw) > 0)); c <- range(which(colSums(nw) > 0))
  H <- dim(a)[1]; W <- dim(a)[2]
  png::writePNG(a[r[1]:r[2], c[1]:c[2], ], png_f)
  fwrite(data.table(b$pts, elev = e, px = w[, 1] * W - c[1] + 1, py = (1 - w[, 2]) * H - r[1] + 1,
                    w = c[2] - c[1] + 1, h = r[2] - r[1] + 1), xy_f)
}
for (k in names(B)) { message(k); render_block(k) }

# =====================================================================================================================
# Step 2 - SVG composition (mm units, 170 mm wide), then Inkscape export
# =====================================================================================================================
W <- 170
# type scale matched to Figs 2-5 (ggplot theme_classic, base_size 7, ragg): 1 pt = 0.3528 mm
pt <- 0.3528
s_tag <- 9 * pt; s_strip <- 7 * pt; s_title <- 7 * pt; s_sub <- 6 * pt; s_axis <- 5.6 * pt; s_lab <- 1.9   # s_lab = geom_text size
fnt <- "font-family:Arial,Helvetica,sans-serif"
col_fl <- "#D55E00"; col_cr <- "#0072B2"; ink <- "black"; grey <- "#4D4D4D"                    # grey30, as in the other figures
halo <- "paint-order:stroke;stroke:white;stroke-width:0.5px;stroke-linejoin:round"
lw <- 0.19                                                                                   # ggplot axis line (0.5 linewidth), mm
svg <- character()
add <- function(...) svg <<- c(svg, sprintf(...))
txt <- function(x, y, s, size = s_lab, anchor = "start", weight = "normal", fill = ink, style = "", rot = 0)
  add('<text x="%.2f" y="%.2f" transform="rotate(%g %.2f %.2f)" style="%s;font-size:%.2fpx;font-weight:%s;fill:%s;text-anchor:%s;%s">%s</text>',
      x, y, rot, x, y, fnt, size, weight, fill, anchor, style, s)
ln <- function(x1, y1, x2, y2, col = ink, w = lw)
  add('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" style="stroke:%s;stroke-width:%.2f"/>', x1, y1, x2, y2, col, w)
rc <- function(x, y, w, h, fill, op = 1)
  add('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" style="fill:%s;fill-opacity:%.2f"/>', x, y, w, h, fill, op)
dot <- function(x, y, fill = "white", r = 0.8)
  add('<circle cx="%.2f" cy="%.2f" r="%.2f" style="fill:%s;stroke:%s;stroke-width:0.3"/>', x, y, r, fill, ink)
head_ <- function(y, tag, s) { txt(1, y, tag, s_tag, weight = "bold"); txt(5, y - 0.15, s, s_sub, fill = grey) }   # tag + subtitle
strip <- function(x, y, s, anchor = "middle") txt(x, y, s, s_strip, anchor, "bold")                        # facet-strip style
badge <- function(x, y, n) {                                     # numbered marker shared by the strip titles and the locator map
  add('<circle cx="%.2f" cy="%.2f" r="1.25" style="fill:white;stroke:%s;stroke-width:0.3"/>', x, y, ink)
  txt(x, y + 0.72, n, s_axis, "middle", "bold") }
strip_n <- function(x, y, s, n) { strip(x, y, s); badge(x - nchar(s) * s_strip * 0.27 - 2, y - 0.85, n) }   # ~ half text width
img <- function(k, x0, y0, wmm) {                                # place a render; returns its marked points in mm
  f <- file.path(rd, paste0(k, ".png")); p <- fread(file.path(rd, paste0(k, "_pts.csv")))
  hmm <- wmm * p$h[1] / p$w[1]; s <- wmm / p$w[1]
  add('<image x="%.2f" y="%.2f" width="%.2f" height="%.2f" preserveAspectRatio="none" xlink:href="data:image/png;base64,%s"/>',
      x0, y0, wmm, hmm, jsonlite::base64_enc(readBin(f, "raw", file.size(f))))
  p[, `:=`(X = x0 + px * s, Y = y0 + py * s)]
  list(p = split(p, by = "key"), h = hmm)
}
cbar <- function(x, y, w, cols, lab, at, labs) {                 # horizontal colour bar, ggplot legend sizes
  n <- 60; g <- colorRampPalette(cols)(n)
  for (i in seq_len(n)) rc(x + (i - 1) * w / n, y, w / n + 0.05, 1.3, g[i])
  for (j in seq_along(at)) { ln(x + at[j] * w, y + 1.3, x + at[j] * w, y + 1.8); txt(x + at[j] * w, y + 3.6, labs[j], s_axis, "middle", fill = grey) }
  txt(x, y - 0.9, lab, s_sub)
}
# Janzen bars: annual ranges [bio6, bio5] of flank and crest on one temperature axis (ggplot-like axis)
therm <- function(x, y, h, fl, cr) {
  T0 <- -12; T1 <- 32; ty <- function(t) y + h * (T1 - t) / (T1 - T0)
  for (t in seq(-10, 30, 10)) { ln(x, ty(t), x - 0.6, ty(t)); txt(x - 1, ty(t) + 0.7, t, s_axis, "end", fill = grey) }
  ln(x, ty(T0), x, ty(T1)); ln(x, ty(T0), x + 14, ty(T0))
  txt(x - 5.2, y + h / 2, "Temperature (°C)", s_title, "middle", rot = -90)
  lo <- max(fl[1], cr[1]); hi <- min(fl[2], cr[2])
  if (hi > lo) rc(x + 1, ty(hi), 13, ty(lo) - ty(hi), "#999999", 0.3)
  rc(x + 2, ty(fl[2]), 3, ty(fl[1]) - ty(fl[2]), col_fl, 0.9)
  rc(x + 9, ty(cr[2]), 3, ty(cr[1]) - ty(cr[2]), col_cr, 0.9)
  txt(x + 3.5, ty(T0) + 2.4, "Flank", s_axis, "middle", fill = grey); txt(x + 10.5, ty(T0) + 2.4, "Crest", s_axis, "middle", fill = grey)
  txt(x + 7, ty(T0) + 4.8, if (hi > lo) sprintf("overlap %.1f °C", hi - lo) else sprintf("gap %.1f °C", lo - hi), s_lab, "middle")
}
v <- function(id, sd, col) m2[col_id == id & side_lab == sd][[col]]   # "side" is a column name
latl <- function(y) ifelse(y < 0, sprintf("%.1f°S", -y), sprintf("%.1f°N", y))

# ---- a: Janzen ----
head_(5, "a", "Janzen: at similar crossing height, the annual range of the lowland flank reaches the temperate crest but not the tropical one")
for (j in 1:2) {
  k <- c("trop", "temp")[j]; id <- c("APD00868", "APD03317")[j]; x0 <- (j - 1) * 85; g <- m2[col_id == id][1]
  strip_n(x0 + 31, 11, sprintf("%s, %s, col %.0f m", c("Tropical", "Temperate")[j], latl(g$lat), g$elev_cross), j)
  r <- img(k, x0 + 1, 13.5, 58)
  for (q in c("crest", "flank_W")) {
    pq <- r$p[[q]]; dot(pq$X, pq$Y, ifelse(q == "crest", col_cr, col_fl))
    txt(pq$X + 1.4, pq$Y - 1, ifelse(q == "crest", "crest", "flank, 1000 m"), style = halo)
  }
  therm(x0 + 70, 13, 29, c(v(id, "W", "flank_bio6"), v(id, "W", "flank_bio5")), c(v(id, "W", "crest_bio6"), v(id, "W", "crest_bio5")))
}
cbar(4, 48, 30, pal$bio1$col, "Mean annual temperature (°C)", (c(-5, 5, 15, 25) + 8) / 36, c("-5", "5", "15", "25"))
txt(127, 49.5, "Bars: annual range, coldest-month minimum to warmest-month maximum", s_sub, "middle", fill = grey)

# ---- b: moisture ----
yB <- 58
head_(yB, "b", "Moisture: the Altiplano crest is far wetter than its Pacific flank; Táchira is humid throughout")
strip_n(55, yB + 6, "Altiplano, 17.7°S, col 4257 m", 3); strip_n(141, yB + 6, "Táchira depression, 7.7°N, col 1023 m", 4)
aiv <- function(id, s) if (s == "crest") v(id, "W", "crest_ai") else v(id, sub("flank_", "", s), "flank_ai")
lab_ai <- function(r, id, keys) for (q in keys) {
  pq <- r$p[[q]]; dot(pq$X, pq$Y); ln(pq$X, pq$Y - 0.8, pq$X, pq$Y - 3)
  txt(pq$X, pq$Y - 3.6, sprintf("%s, AI %s", c(crest = "crest", flank_W = "W flank", flank_E = "E flank")[[q]],
                                formatC(aiv(id, q), digits = 2, format = "fg")), anchor = "middle", style = halo)
}
ra <- img("alti", 1, yB + 8.5, 106); lab_ai(ra, "APD02335", c("flank_W", "crest", "flank_E"))
rt <- img("tach", 112, yB + 8.5, 56); lab_ai(rt, "APD00203", "crest")
txt(rt$p$crest$X, rt$p$crest$Y + 4, "flanks, AI 1.1 and 1.2", anchor = "middle", style = halo)
cbar(4, yB + 8.5 + max(ra$h, rt$h) + 3, 40, pal$ai$col, "Aridity index, P / PET (Global-AI v3.1; log scale)",
     (log10(c(0.003, 0.03, 0.3, 2)) + 2.7) / 3.1, c("0.003", "0.03", "0.3", "2"))

# ---- c: warming ----
yC <- yB + 8.5 + max(ra$h, rt$h) + 12
head_(yC, "c", sprintf("Warming: the isotherm now at %.0f m rises %.0f m; low crests open west–east and break the montane belt north–south, high crests do not", thr, rise_c))
xC <- 37; wC <- 131
strip(xC, yC + 5.5, "Present (1970–2000)", "start")
rn <- img("now", xC, yC + 7, wC)
yF <- yC + 7 + rn$h + 8.5
add('<path d="M %.2f %.2f l 0 4.2 m -1.2 -1.5 l 1.2 1.5 l 1.2 -1.5" style="fill:none;stroke:%s;stroke-width:0.4"/>', xC + wC / 2, yF - 8, ink)
txt(xC + wC / 2 + 2.2, yF - 5, sprintf("isotherm +%.0f m", rise_c), fill = grey)
strip(xC, yF - 1.5, "2081–2100, SSP3-7.0 (mean of 8 GCMs)", "start")
rf <- img("fut", xC, yF, wC)
call_ <- function(pq, s) { dot(pq$X, pq$Y); ln(pq$X, pq$Y - 0.8, pq$X, pq$Y - 5); txt(pq$X, pq$Y - 5.6, s, anchor = "middle", style = halo) }
call_(rn$p$lost, "montane belt continuous")
call_(rf$p$lost, "west–east crossing opens; montane belt breaks")
call_(rf$p$high, "high crests keep montane climate")
for (r in list(rn, rf)) for (q in c("Colombia", "Ecuador", "Peru"))
  txt(r$p[[q]]$X, r$p[[q]]$Y, toupper(q), s_axis, "middle", fill = grey, style = paste0(halo, ";letter-spacing:0.25px"))
yL <- yF + rf$h + 2.5; H <- yL + 3.5
txt(xC, yL + 1, "N (0.6°N)", s_axis, fill = grey); txt(xC + wC, yL + 1, "S (10°S)", s_axis, "end", fill = grey)
rc(xC + 30, yL - 0.9, 3, 2.2, pal$belt$col[1]); txt(xC + 34, yL + 1, "Lowland climate (below the isotherm)", s_axis)
rc(xC + 80, yL - 0.9, 3, 2.2, pal$belt$col[2]); txt(xC + 84, yL + 1, "Montane climate (above it)", s_axis)

# ---- locator map ----
# locator map spans the two blocks of c; numbered markers match the strip titles of a and b
k_ <- (yF + rf$h - (yC + 7)) / 48                                # mm per degree
mx <- function(lon) 1 + (lon + 82) * k_; my <- function(lat) yC + 7 + (12 - lat) * k_
sa1 <- sf::st_as_sf(maps::map("world", fill = TRUE, plot = FALSE)) |> sf::st_make_valid() |>
  sf::st_crop(xmin = -82, xmax = -60, ymin = -36, ymax = 12) |> sf::st_geometry() |> sf::st_cast("MULTIPOLYGON") |> sf::st_cast("POLYGON")
for (pg in sa1) { xy <- pg[[1]]; add('<path d="M %s Z" style="fill:#e2e2e2;stroke:white;stroke-width:0.25;stroke-linejoin:round"/>', paste(sprintf("%.2f %.2f", mx(xy[, 1]), my(xy[, 2])), collapse = " L ")) }
e <- B$now$ext
add('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" style="fill:none;stroke:%s;stroke-width:0.35"/>', mx(e[1]), my(e[4]), mx(e[2]) - mx(e[1]), my(e[3]) - my(e[4]), ink)
txt(mx(e[1]) + (mx(e[2]) - mx(e[1])) / 2, my(e[3]) + 2.6, "area in c", s_axis, "middle", "bold")
for (q in list(list(-78.6, -0.65, 1, "Tropical (a)"), list(-69.8, -34.24, 2, "Temperate (a)"), list(-69.67, -17.67, 3, "Altiplano (b)"),
               list(-72.42, 7.69, 4, "Táchira (b)"))) {
  badge(mx(q[[1]]), my(q[[2]]), q[[3]]); txt(mx(q[[1]]) + 1.8, my(q[[2]]) + 0.7, q[[4]], s_axis, style = halo) }

out <- file.path(cfg$paths$outputs, "figures", "Fig5_synthesis")
writeLines(c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="%gmm" height="%gmm" viewBox="0 0 %g %g">', W, H, W, H),
             sprintf('<rect width="%g" height="%g" fill="white"/>', W, H), svg, "</svg>"), paste0(out, ".svg"), useBytes = TRUE)
ink_exe <- Sys.getenv("INKSCAPE", Sys.which("inkscape"))
if (!nzchar(ink_exe)) ink_exe <- "C:/Program Files/Inkscape/bin/inkscape.com"   # Windows default; or set INKSCAPE
system2(ink_exe, c(shQuote(paste0(out, ".svg")), "--export-type=pdf", "--export-text-to-path", paste0("--export-filename=", shQuote(paste0(out, ".pdf")))))
system2(ink_exe, c(shQuote(paste0(out, ".svg")), "--export-type=png", "--export-dpi=600", "--export-background=white",
                   paste0("--export-filename=", shQuote(paste0(out, ".png")))))
message("written ", out, ".{svg,pdf,png}")
