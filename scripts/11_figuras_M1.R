# 11_figuras_M1.R - M1 figures (English): continental crossing-height profile with gateways (+ South America
# inset), the Huancabamba system zoom with a relief map, and elevation vs isolation of tropical cols (+ inset).
# Okabe-Ito colours plus shape (never colour alone); the profile line is broken at real gaps (>= 2 km);
# only key cols are labelled (ggrepel). Point names come from APD_cols$name (proper names, then OSM).
# Outputs PNG (600 dpi) and vector PDF in outputs/figs_qc/.
# Rscript scripts/11_figuras_M1.R
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })
suppressMessages(sf::sf_use_s2(FALSE))
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; P <- cfg$params
tr_lat <- unlist(P$tropics_lat)   # yaml gives a list for mixed int / float
fq <- "outputs/figs_qc"; dir.create(fq, showWarnings = FALSE)
prof <- read.csv(file.path(int, "M1_profile.csv"))
cols <- read.csv("outputs/database/APD_cols.csv", encoding = "UTF-8")
oi <- c(verm = "#D55E00", blue = "#0072B2", grey = "#999999", black = "#000000", yell = "#E69F00")
th <- theme_classic(base_size = 8) +
  theme(legend.position = "top", legend.title = element_blank(), legend.key.size = unit(3, "mm"),
        plot.title = element_text(size = 9, face = "bold"), axis.title = element_text(size = 8))
save2 <- function(p, name, w, h) {
  ggsave(file.path(fq, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(fq, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
}
latlab <- function(y, d = 1) ifelse(y > 0, sprintf(paste0("%.", d, "f°N"), y), ifelse(y < 0, sprintf(paste0("%.", d, "f°S"), -y), "0°"))
# Regional descriptors for key gateways that have no proper name (lat, lon of the col)
reg <- data.frame(n = c("Táchira", "Altiplano"), lat = c(7.686, -20.367), lon = c(-72.416, -68.710))
cols$lab <- sub("^Abra (de |el )?", "", cols$name)     # short form for labels: "Abra de Porculla" -> "Porculla"
for (k in seq_len(nrow(reg))) {
  j <- which.min(abs(cols$lat - reg$lat[k]) + abs(cols$lon - reg$lon[k])); if (is.na(cols$lab[j])) cols$lab[j] <- reg$n[k]
}
sa <- sf::st_as_sf(maps::map("world", fill = TRUE, plot = FALSE))
sa <- suppressWarnings(sf::st_crop(sf::st_make_valid(sa), xmin = -82, xmax = -34, ymin = -56, ymax = 13))

# --- Continental profile ---
gw <- cols[cols$gateway, ]
lv <- c("Tropical gateway (11°N–23.5°S)", "Extratropical gateway", "Gateway at the end of the chain")
gw$cls <- factor(ifelse(gw$terminal, lv[3], ifelse(gw$tropical, lv[1], lv[2])), levels = lv)
key <- gw[gw$tropical & !gw$terminal, ]
key <- key[order(key$elev_cross), ][1:4, ]
key$lab <- sprintf("%s\n%.0f m", key$lab, key$elev_cross)
p1 <- ggplot() +
  annotate("rect", xmin = tr_lat[2], xmax = tr_lat[1], ymin = -Inf, ymax = Inf, fill = "grey93") +
  annotate("text", x = mean(tr_lat), y = 6900, label = "Tropical Andes", size = 2.6, colour = "grey40") +
  geom_line(data = prof, aes(lat, elev_cross, group = grp), colour = "grey55", linewidth = 0.2) +
  geom_point(data = gw, aes(lat, elev_cross, colour = cls, shape = cls, fill = cls), size = 1.6, stroke = 0.5) +
  geom_text_repel(data = key, aes(lat, elev_cross, label = lab), size = 2.3, lineheight = 0.9,
                  ylim = c(-650, 600), direction = "x", force = 3, box.padding = 0.6,
                  segment.size = 0.25, segment.colour = "grey30", min.segment.length = 0, seed = 1) +
  scale_x_reverse(breaks = seq(10, -55, -5), labels = function(x) latlab(x, 0)) +
  scale_y_continuous(limits = c(-700, 7000), breaks = seq(0, 7000, 1000), expand = c(0, 0)) +
  scale_colour_manual(values = unname(oi[c("verm", "blue", "grey")]), drop = FALSE) +
  scale_fill_manual(values = c(unname(oi["verm"]), unname(oi["blue"]), "white"), drop = FALSE) +
  scale_shape_manual(values = c(21, 24, 21), drop = FALSE) +
  labs(x = "Latitude (north to south)", y = "W–E crossing height (m)",
       title = "W|E divide of the Andes (Copernicus GLO-30): crossing-height profile and gateways") +
  th
# Inset: South America with the divide coloured tropical / extratropical; no frame, transparent background
prof$zone <- ifelse(prof$lat >= tr_lat[2], "trop", "extra")
prof$seg <- paste(prof$grp, cumsum(c(TRUE, diff(prof$zone == "trop") != 0)))
map <- ggplot() +
  geom_sf(data = sa, fill = "grey88", colour = NA) +
  geom_hline(yintercept = tr_lat[2], colour = "grey55", linewidth = 0.2, linetype = "22") +
  geom_path(data = prof, aes(lon, lat, group = seg, colour = zone), linewidth = 0.45, show.legend = FALSE) +
  scale_colour_manual(values = c(trop = unname(oi["verm"]), extra = unname(oi["blue"]))) +
  coord_sf(xlim = c(-82, -34), ylim = c(-56, 13), expand = FALSE) +
  theme_void() + theme(plot.background = element_blank(), panel.background = element_blank())
# legend to the left so the map sits in the free top-right corner beside title and legend
p1 <- p1 + theme(legend.justification = "left") +
  patchwork::inset_element(map, left = 0.83, bottom = 0.50, right = 1.00, top = 1.00, align_to = "full")
save2(p1, "M1_perfil_continental", 180, 85)

# --- Gateway systems, north to south: divide profile + relief map (A Táchira, B S Colombia, C Huancabamba) ---
nmr <- read.csv("data/refs/M1_col_names.csv", encoding = "UTF-8")
adm_all <- sf::st_make_valid(do.call(rbind, lapply(c("VEN", "COL", "ECU", "PER"), function(cc)
  sf::st_as_sf(terra::unwrap(readRDS(file.path(cfg$paths$data_dir, "gadm", "gadm", sprintf("gadm41_%s_1_pk.rds", cc)))))[, c("GID_0", "NAME_1")])))
ctry_nm <- c(VEN = "VENEZUELA", COL = "COLOMBIA", ECU = "ECUADOR", PER = "PERU")
cu_all <- lapply(split(adm_all, adm_all$GID_0), sf::st_union)
# international borders only: shared edges of each pair of countries (not coasts or crop frames)
borders <- do.call(c, lapply(combn(names(cu_all), 2, simplify = FALSE), function(pr)
  suppressWarnings(sf::st_intersection(sf::st_boundary(cu_all[[pr[1]]]), sf::st_boundary(cu_all[[pr[2]]])))))
borders <- suppressWarnings(sf::st_collection_extract(sf::st_sf(geometry = borders), "LINESTRING"))
aea <- "+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42"
hill <- function(id, bb) {             # GLO-30 hillshade aggregated to ~240 m, cached; fetches missing tiles
  f <- file.path(int, sprintf("M1_hillshade_%s.tif", id))
  if (!file.exists(f)) {
    tl <- expand.grid(la = floor(bb[["ymin"]]):floor(bb[["ymax"]]), lo = floor(bb[["xmin"]]):floor(bb[["xmax"]]))
    nm <- sprintf("Copernicus_DSM_COG_10_%s%02d_00_%s%03d_00_DEM", ifelse(tl$la < 0, "S", "N"), abs(tl$la),
                  ifelse(tl$lo < 0, "W", "E"), abs(tl$lo))
    fs <- file.path(cfg$paths$data_dir, "glo30", paste0(nm, ".tif"))
    for (k in which(!file.exists(fs)))
      try(download.file(sprintf("%s/%s/%s.tif", cfg$remote$glo30_bucket, nm[k], nm[k]), fs[k], mode = "wb", quiet = TRUE), silent = TRUE)
    unlink(fs[file.exists(fs) & file.size(fs) == 0])
    dem <- terra::aggregate(terra::crop(terra::vrt(fs[file.exists(fs)]), terra::ext(bb[c("xmin", "xmax", "ymin", "ymax")])), 8, "mean")
    dem[dem <= 0] <- NA                                                                  # sea
    terra::writeRaster(terra::shade(terra::terrain(dem, "slope", unit = "radians"),
                                    terra::terrain(dem, "aspect", unit = "radians"), 40, 315), f, overwrite = TRUE)
  }
  terra::rast(f)
}
lv2 <- c("Gateway (lowest within ±100 km)", "Col (depth ≥ 50 m)")
system_fig <- function(id, d_from, d_to, anchor, title, pad = 0.3, ctry_xy = NULL) {
  z <- prof[prof$d_km >= d_from & prof$d_km <= d_to, ]
  k <- cols[cols$d_km >= d_from & cols$d_km <= d_to & cols$depth >= 50, ]
  d0 <- cols$d_km[which(cols$lab %in% anchor)][1]
  z$x <- z$d_km - d0; k$x <- k$d_km - d0
  k$cls <- factor(ifelse(k$gateway, lv2[1], lv2[2]), levels = lv2)
  # labels: gateways (proper / regional name, else latitude) and the proper names of the names table
  k$lb <- ifelse(is.na(k$lab), latlab(k$lat, 2), k$lab)
  lb <- k[k$gateway | k$name %in% nmr$name, ]; lb$txt <- sprintf("%s\n%.0f m", lb$lb, lb$elev_cross)
  # latitude on the top axis: first crossing of each half degree
  lt <- seq(ceiling(min(z$lat) * 2) / 2, floor(max(z$lat) * 2) / 2, 0.5)
  lx <- sapply(lt, function(y) { j <- which.min(abs(z$lat - y)); if (abs(z$lat[j] - y) < 0.02) z$x[j] else NA })
  lt <- lt[!is.na(lx)]; lx <- lx[!is.na(lx)]
  yl <- c(floor((min(z$elev_cross, na.rm = TRUE) - 500) / 250) * 250, max(z$elev_cross, na.rm = TRUE) + 150)
  p <- ggplot() +
    geom_line(data = z, aes(x, elev_cross, group = grp), colour = "grey55", linewidth = 0.3) +
    geom_point(data = k[k$cls == lv2[2], ], aes(x, elev_cross, colour = cls, shape = cls), size = 0.9) +
    geom_point(data = k[k$cls == lv2[1], ], aes(x, elev_cross, colour = cls, shape = cls), size = 2, fill = oi["verm"]) +
    geom_text_repel(data = lb, aes(x, elev_cross, label = txt), size = 2.3, lineheight = 0.9,
                    ylim = c(yl[1], min(z$elev_cross, na.rm = TRUE) - 120), direction = "x", force = 2, box.padding = 0.4,
                    segment.size = 0.25, segment.colour = "grey30", min.segment.length = 0, seed = 1) +
    scale_colour_manual(values = unname(oi[c("verm", "black")]), drop = FALSE) +
    scale_shape_manual(values = c(24, 1), drop = FALSE) +
    scale_x_continuous(sec.axis = dup_axis(breaks = lx, labels = latlab(lt, 1), name = NULL)) +
    scale_y_continuous(limits = yl, breaks = seq(0, 7000, 500)) +
    labs(x = sprintf("Distance along the divide from %s (km; negative = north)", anchor), y = "W–E crossing height (m)",
         title = title) +
    th
  # map
  bb <- c(xmin = min(z$lon) - pad, xmax = max(z$lon) + pad, ymin = min(z$lat) - pad, ymax = max(z$lat) + pad)
  adm <- suppressWarnings(sf::st_crop(adm_all, bb))
  bd <- suppressWarnings(sf::st_crop(borders, bb))
  big <- as.numeric(sf::st_area(sf::st_transform(adm, aea))) > 2e9                     # skip slivers at the map edge
  al <- suppressWarnings(sf::st_point_on_surface(adm[big, ])); al <- cbind(al, sf::st_coordinates(al))
  if (is.null(ctry_xy)) {                                                              # country label at the centre of its part
    cc <- suppressWarnings(sf::st_point_on_surface(aggregate(adm["GID_0"], list(g = adm$GID_0), function(v) v[1])))
    ctry_xy <- data.frame(sf::st_coordinates(cc), lab = ctry_nm[cc$g])
  }
  hs <- terra::as.data.frame(terra::mask(hill(id, bb), terra::vect(adm)), xy = TRUE); names(hs)[3] <- "v"
  hg <- k[k$gateway, ]
  ml <- rbind(data.frame(X = al$X, Y = al$Y, t = al$NAME_1, ff = "italic", col = "grey15", sz = 1.7, fl = NA),
              data.frame(X = ctry_xy$X, Y = ctry_xy$Y, t = ctry_xy$lab, ff = "bold", col = "grey10", sz = 2.3, fl = alpha("white", 0.7)),
              data.frame(X = hg$lon, Y = hg$lat, t = hg$lb, ff = "bold", col = unname(oi["verm"]), sz = 1.9, fl = alpha("white", 0.85)))
  m <- ggplot() +
    geom_raster(data = hs, aes(x, y, fill = v), show.legend = FALSE) +
    scale_fill_gradient(low = "grey35", high = "grey97", na.value = NA) +
    geom_sf(data = adm, fill = NA, colour = "white", linewidth = 0.15) +
    geom_sf(data = bd, colour = "grey15", linewidth = 0.25, linetype = "22") +
    geom_path(data = z, aes(lon, lat, group = grp), colour = "black", linewidth = 0.6) +
    geom_point(data = hg, aes(lon, lat), shape = 24, size = 1.8, colour = oi["verm"], fill = oi["verm"]) +
    geom_label_repel(data = ml, aes(X, Y, label = t, fontface = ff, colour = I(col), size = I(sz), fill = I(fl)),
                     linewidth = 0, label.padding = unit(0.4, "mm"), box.padding = 0.25, point.padding = 0.15,
                     xlim = bb[1:2] + c(0.04, -0.04), ylim = bb[3:4] + c(0.04, -0.04), max.overlaps = Inf,
                     segment.size = 0.2, segment.colour = "grey30", min.segment.length = 0.3, seed = 2, show.legend = FALSE) +
    coord_sf(xlim = bb[1:2], ylim = bb[3:4], expand = FALSE, label_axes = "-N--") +
    scale_y_continuous(breaks = lt, labels = latlab(lt, 1)) +
    theme_void(base_size = 7) +
    theme(axis.text.y.right = element_text(size = 6, colour = "grey20", margin = margin(l = 2)),
          plot.background = element_blank(), panel.background = element_blank())
  list(p = p, m = m)
}
# country labels placed by hand on free ground (auto placement fell under the gateway labels)
sysA <- system_fig("tachira", 360, 860, "Táchira", "A  Táchira depression (Venezuela–Colombia)",
                   ctry_xy = data.frame(X = c(-70.9, -73.05), Y = c(7.45, 6.95), lab = c("VENEZUELA", "COLOMBIA")))
sysB <- system_fig("scolombia", 1380, 1980, "Aguas Claras", "B  Meta-Huila and Aguas Claras (Colombian Andes)")
d0C <- cols$d_km[which(cols$lab %in% "Porculla")][1]
sysC <- system_fig("huancabamba", d0C - 300, d0C + 340, "Porculla", "C  Huancabamba system (northern Peru–southern Ecuador)",
                   ctry_xy = data.frame(X = c(-78.6, -78.45), Y = c(-4.25, -5.35), lab = c("ECUADOR", "PERU")))
for (s in list(list(sysA, "M1_perfil_Tachira"), list(sysB, "M1_perfil_SurColombia"), list(sysC, "M1_perfil_Huancabamba")))
  save2(s[[1]]$p + s[[1]]$m + patchwork::plot_layout(widths = c(3.7, 1.3)), s[[2]], 200, 100)
row <- function(s) s$p + s$m + patchwork::plot_layout(widths = c(3.7, 1.3))
abc <- patchwork::wrap_plots(row(sysA), row(sysB), row(sysC), ncol = 1) + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Low W–E gateways of the tropical Andes, north to south (Copernicus GLO-30)",
                             theme = theme(plot.title = element_text(size = 10, face = "bold"))) &
  theme(legend.position = "top")
save2(abc, "M1_perfiles_ABC", 200, 290)

# --- Elevation vs isolation, tropical cols (depth >= 100 m, chain ends excluded) ---
tr <- cols[cols$tropical & cols$depth >= P$depth_m[[1]] & !cols$terminal, ]
lv3 <- c("Gateway", "Col", "Lara (NE end of the Mérida Andes)")
tr$cls <- factor(ifelse(tr$lara, lv3[3], ifelse(tr$gateway, lv3[1], lv3[2])), levels = lv3)
# key cols only: named / regional gateways that are low (< 2500 m) or very isolated (>= 1000 km)
tr$key <- !is.na(tr$lab) & tr$gateway & (tr$elev_cross < 2500 | tr$iso_km >= 1000) & !tr$lara
tr$txt <- ifelse(tr$key, sprintf("%s (%s)\n%.0f m · %.0f km", tr$lab, latlab(tr$lat, 1), tr$elev_cross, tr$iso_km), "")
p3 <- ggplot(tr, aes(iso_km, elev_cross)) +
  geom_point(data = tr[tr$cls == lv3[2], ], aes(colour = cls, shape = cls, fill = cls), size = 0.8, alpha = 0.5) +
  geom_point(data = tr[tr$cls != lv3[2], ], aes(colour = cls, shape = cls, fill = cls), size = 1.8, stroke = 0.4) +
  # every col is passed (empty label for most) so the labels steer clear of all points
  geom_text_repel(aes(label = txt), size = 2.1, lineheight = 0.9, box.padding = 0.5, point.padding = 0.2, force = 3,
                  max.overlaps = Inf, segment.size = 0.25, segment.colour = "grey30", min.segment.length = 0, seed = 1) +
  scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1000")) +
  annotation_logticks(sides = "b", linewidth = 0.2, short = unit(0.6, "mm"), mid = unit(1, "mm"), long = unit(1.4, "mm")) +
  scale_y_continuous(breaks = seq(1000, 6000, 1000)) +
  scale_colour_manual(values = c(unname(oi["verm"]), "grey45", unname(oi["yell"])), drop = FALSE) +
  scale_fill_manual(values = c(unname(oi["verm"]), "grey45", unname(oi["yell"])), drop = FALSE) +
  scale_shape_manual(values = c(21, 16, 24), drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  labs(x = "Isolation along the divide (km, log scale)", y = "W–E crossing height (m)",
       title = "Tropical Andean cols: crossing height vs isolation",
       caption = sprintf(paste0("n = %d cols with depth ≥ 100 m (11°N–23.5°S, chain ends excluded). Isolation = distance along the divide",
                                "\nto the nearest lower point (minimum of N and S), measured with a 0.5 km ruler."), nrow(tr))) +
  th + theme(plot.caption = element_text(size = 6, hjust = 0, colour = "grey30"))
# Side map: the tropical Andes only, with the gateways (orange), Lara (yellow) and the labelled cols
tz <- prof[prof$lat >= tr_lat[2] & prof$lat <= tr_lat[1], ]
kk <- tr[tr$key, ]
m3 <- ggplot() +
  geom_sf(data = sa, fill = "grey88", colour = NA) +
  geom_path(data = tz, aes(lon, lat, group = grp), colour = "grey40", linewidth = 0.35) +
  geom_point(data = tr[tr$cls == lv3[1], ], aes(lon, lat), shape = 21, size = 1.5, stroke = 0.3, colour = "black", fill = oi["verm"]) +
  geom_point(data = tr[tr$lara, ], aes(lon, lat), shape = 24, size = 1.5, stroke = 0.3, colour = "black", fill = oi["yell"]) +
  # labels in a tidy column east of the cordillera
  geom_text_repel(data = kk, aes(lon, lat, label = lab), size = 2.1, direction = "y", nudge_x = pmax(-71, kk$lon + 1.5) - kk$lon, hjust = 0,
                  segment.size = 0.2, segment.colour = "grey30", min.segment.length = 0, box.padding = 0.2, seed = 1) +
  coord_sf(xlim = c(-82, -60), ylim = c(-24.5, 12.5), expand = FALSE) +
  theme_void() + theme(plot.background = element_blank(), panel.background = element_blank())
p3 <- p3 + m3 + patchwork::plot_layout(widths = c(3, 1.15))
save2(p3, "M1_elev_vs_aislamiento", 190, 120)
cat("ok\n")

# --- Drainage domains with the final GLO-30 W|E divide (replaces the HydroBASINS QC map of 02) ---
# A: whole Andes; B-D: zooms N to S, labelling all 30 non-terminal tropical gateways. Divide drawn from M1_profile.csv, broken at real gaps (grp).
dom <- sf::read_sf(file.path(int, "M1_domains.gpkg"))
dom$dom <- factor(dom$dom, levels = c("W", "E", "I"),
                  labels = c("Pacific + Caribbean (W)", "Atlantic (E)", "Endorheic (I)"))
ie <- sf::read_sf(file.path(int, "M1_divide_raw.gpkg")); ie <- ie[ie$type == "I|E", ]
andes <- sf::read_sf(file.path(int, "M1_gmba_andes.gpkg"))
stp <- c(Inf, sqrt(diff(prof$lon)^2 + diff(prof$lat)^2))
prof$mgrp <- cumsum(c(TRUE, diff(prof$grp) != 0) | stp > 0.03)
llv <- c("W|E divide (GLO-30)", "I|E divide (HydroBASINS)", "Andes (GMBA v2)")
dmap <- function(xl, yl) {
  p <- ggplot() +
    geom_sf(data = sa, fill = "grey92", colour = "white", linewidth = 0.15) +
    geom_sf(data = dom, aes(fill = dom), colour = NA, alpha = 0.55) +
    geom_sf(data = andes, aes(colour = llv[3], linetype = llv[3]), fill = NA, linewidth = 0.12, key_glyph = "path") +
    geom_sf(data = ie, aes(colour = llv[2], linetype = llv[2]), linewidth = 0.35, key_glyph = "path") +
    geom_path(data = prof, aes(lon, lat, group = mgrp, colour = llv[1], linetype = llv[1]), linewidth = 0.45) +
    geom_hline(yintercept = tr_lat[2], colour = "grey40", linewidth = 0.2, linetype = "22") +
    scale_fill_manual(values = setNames(c("#E69F00", "#56B4E9", "#CC79A7", oi[["verm"]], oi[["blue"]], "white"),
                                        c(levels(dom$dom), lv)), breaks = levels(dom$dom), name = "Drainage domain") +
    scale_colour_manual(values = setNames(c("black", "#882255", "grey45"), llv), breaks = llv, name = NULL) +
    scale_linetype_manual(values = setNames(c("solid", "22", "solid"), llv), breaks = llv, name = NULL) +
    coord_sf(xlim = xl, ylim = yl, expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 7) + theme(panel.grid = element_line(colour = "grey85", linewidth = 0.15),
                                    panel.background = element_rect(fill = "#F4F8FB"))
  p
}
# Zoom extents; each gateway is labelled once, in the panel whose latitude band (lab) holds it
zm <- list(B = list(x = c(-81, -67), y = c(-7, 11.5), lab = c(-7, 12), t = "B  Northern Andes", nx = 3),
           C = list(x = c(-80.5, -68.5), y = c(-15.6, -6.5), lab = c(-15.5, -7), t = "C  Central Andes of Peru", nx = 1.5),
           D = list(x = c(-72, -63.5), y = c(-24.5, -15.5), lab = c(-24.5, -15.5), t = "D  Altiplano–Puna plateau", nx = 1.5))
gw$txt <- sprintf("%s%s %.0f m", toupper(substr(gw$lab, 1, 1)), substring(gw$lab, 2), gw$elev_cross)  # "la Viuda" -> "La Viuda"
zoom <- function(z) {
  g <- gw[gw$lon >= z$x[1] & gw$lon <= z$x[2] & gw$lat >= z$y[1] & gw$lat <= z$y[2] & !gw$terminal & gw$tropical, ]
  l <- g[g$lat > z$lab[1] & g$lat <= z$lab[2], ]
  dmap(z$x, z$y) +
    geom_point(data = g, aes(lon, lat, shape = cls, fill = cls), size = 1.5, stroke = 0.35,
               show.legend = c(shape = TRUE, fill = FALSE)) +
    geom_label_repel(data = l, aes(lon, lat, label = txt), size = 1.9, label.size = 0.1, label.padding = 0.1,
                     fill = alpha("white", 0.85), min.segment.length = 0, segment.size = 0.25, box.padding = 0.25,
                     nudge_x = z$nx, direction = "y", force = 2, max.overlaps = Inf, seed = 1) +
    scale_shape_manual(values = setNames(c(21, 24, 21), lv), name = "Gateway", limits = lv[1],
                       guide = guide_legend(override.aes = list(fill = oi[["verm"]], size = 1.8))) +
    ggtitle(z$t)
}
box <- function(z, t) list(
  annotate("rect", xmin = z$x[1], xmax = z$x[2], ymin = z$y[1], ymax = z$y[2], fill = NA, colour = "black", linewidth = 0.3),
  annotate("text", x = z$x[2] - 0.4, y = z$y[2] - 0.4, label = t, hjust = 1, vjust = 1, size = 2.4, fontface = "bold"))
pA <- dmap(c(-82, -60), c(-56, 12.5)) + box(zm$B, "B") + box(zm$C, "C") + box(zm$D, "D") +
  annotate("text", x = -60.6, y = tr_lat[2] - 0.8, label = "23.5°S", hjust = 1, vjust = 1, size = 2, colour = "grey30") +
  ggtitle("A  Andes")
pD <- zoom(zm$D) + scale_x_continuous(breaks = seq(-72, -64, 2))
pm <- patchwork::wrap_plots(A = pA, B = zoom(zm$B), C = zoom(zm$C), D = pD,
                            design = "ABC
ABD", guides = "collect",
                            # absolute panel sizes (mm) matched to each map's aspect: no slack, titles sit on the frames
                            widths = unit(c(44, 108, 72), "mm"), heights = unit(c(56, 81), "mm")) &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.key.width = unit(6, "mm"), plot.title = element_text(size = 8, face = "bold"))
save2(pm, "M1_dominios", 252, 170)
