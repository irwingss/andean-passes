# 26_figuras_manuscrito.R - Article 1 figures (plan v2.0 section 8.1), English, from the M1-M3 outputs only.
# Style: Okabe-Ito for categories (always paired with shape or line type), the config terrain ramp for elevation and viridis
# for other continuous values, one sans font,
# 170 mm double-column width, PNG 600 dpi + vector PDF (cairo, fonts embedded). No gateway gets an accent colour (plan v2.0):
# the main tropical gateways are labelled alike. Writes outputs/figures/Fig1-Fig4 (.png, .pdf); Fig. 5 (synthesis) is 32_fig5_sintesis.R. A2 figures: L2 26_figuras_A2.R.
# Rscript scripts/26_figuras_manuscrito.R
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(patchwork); library(data.table); library(terra); library(mgcv) })
suppressMessages(sf::sf_use_s2(FALSE))
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; raw <- cfg$paths$data_dir; tb <- file.path(cfg$paths$outputs, "tables")
P <- cfg$params; tr_lat <- unlist(P$tropics_lat)
er <- cfg$style$elev_ramp                                         # one terrain ramp for every elevation colour (config.yml)
scale_elev <- function(aesthetics, name) scale_colour_gradientn(name, colours = er$col, values = scales::rescale(er$m),
  limits = range(er$m), oob = scales::squish, aesthetics = aesthetics)
fd <- file.path(cfg$paths$outputs, "figures"); dir.create(fd, showWarnings = FALSE)
oi <- c(verm = "#D55E00", blue = "#0072B2", sky = "#56B4E9", green = "#009E73", yell = "#E69F00", pink = "#CC79A7", grey = "grey55")
th <- theme_classic(base_size = 7) +
  theme(legend.position = "top", legend.title = element_blank(), legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), plot.tag = element_text(size = 9, face = "bold"),
        strip.background = element_blank(), strip.text = element_text(size = 7, face = "bold"))
save2 <- function(p, name, w, h) {
  ggsave(file.path(fd, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(fd, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
}
latlab <- function(y, d = 0) ifelse(y > 0, sprintf(paste0("%.", d, "f°N"), y), ifelse(y < 0, sprintf(paste0("%.", d, "f°S"), -y), "0°"))
sa <- sf::st_as_sf(maps::map("world", fill = TRUE, plot = FALSE))
sa <- suppressWarnings(sf::st_crop(sf::st_make_valid(sa), xmin = -82, xmax = -34, ymin = -56, ymax = 13))

# ---- shared data ----
prof <- fread(file.path(int, "M1_profile.csv"))
cols <- fread("outputs/database/APD_cols.csv", encoding = "UTF-8")
cols[, lab := sub("^Abra (de |el )?", "", name)]
cols[gateway & is.na(lab) & abs(lat - 7.686) < 0.01 & abs(lon + 72.416) < 0.01, lab := "Táchira"]   # regional name, as in Table 1
m2 <- fread(file.path(tb, "M2_metricas_clima_collados.csv"))[flank_m == P$flank_m[1] & !terminal]
ax <- function(v) if (all(is.na(v))) NA_real_ else max(abs(v), na.rm = TRUE)         # moisture gaps: larger absolute value
cw <- m2[, .(name = name[1], lat = lat[1], lon = lon[1], elev_cross = elev_cross[1], tropical = tropical[1], gateway = gateway[1],
             TO_ann = min(TO_ann), AI_gap = ax(AI_gap), CLD_gap = ax(CLD_gap), Pg0 = ax(log(crest_bio12 / flank_bio12))), by = col_id]
tg <- cw[gateway & tropical]                                                    # the 30 tropical gateways (flank 1000 m)

# =====================================================================================================================
# Fig. 1 - the Andean barrier and its gateways
# =====================================================================================================================
k200 <- cols[depth >= 200]
low10 <- cols[gateway & tropical & !terminal][order(elev_cross)][1:10]
low10[, txt := sprintf("%s %.0f m", ifelse(is.na(lab), latlab(lat, 1), lab), elev_cross)]
f1a <- ggplot() +
  geom_sf(data = sa, fill = "grey90", colour = NA) +
  geom_point(data = k200, aes(lon, lat, size = depth), shape = 21, fill = NA, colour = "grey35", stroke = 0.12, alpha = 0.7) +
  geom_path(data = prof, aes(lon, lat, group = grp, colour = elev), linewidth = 0.7) +
  geom_text_repel(data = low10, aes(lon, lat, label = txt), size = 1.9, nudge_x = 6, direction = "y", hjust = 0,
                  segment.size = 0.2, segment.colour = "grey30", min.segment.length = 0, seed = 1, max.overlaps = Inf) +
  geom_hline(yintercept = tr_lat[2], colour = "grey45", linewidth = 0.2, linetype = "22") +
  scale_elev("colour", "Crest elevation (m)") +
  scale_size("Col depth (m)", range = c(0.2, 3), breaks = c(200, 500, 1000, 2000)) +
  coord_sf(xlim = c(-82, -55), ylim = c(-56, 13), expand = FALSE) +
  labs(x = NULL, y = NULL) + th +
  theme(legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 6),
        legend.key.width = unit(6, "mm"), axis.line = element_blank())
env <- prof[, .(elev_cross = min(elev_cross)), by = .(lat = round(lat * 4) / 4)]   # lower envelope, 0.25° bins
gaps <- fread(file.path(int, "M1_montane_gaps.csv"))[threshold_m == P$montane_gap_m[[1]]]
f1b <- ggplot() +
  geom_rect(data = gaps, aes(xmin = lat_start, xmax = lat_end, ymin = -Inf, ymax = Inf), fill = oi[["sky"]], alpha = 0.25) +
  geom_line(data = prof, aes(lat, elev, group = grp), colour = "grey70", linewidth = 0.12) +
  geom_line(data = env, aes(lat, elev_cross), colour = "black", linewidth = 0.4) +
  geom_point(data = cols[gateway == TRUE], aes(lat, elev_cross), shape = 21, size = 0.9, fill = oi[["blue"]], stroke = 0.15) +
  scale_x_reverse(breaks = seq(10, -55, -10), labels = latlab) +
  scale_y_continuous(breaks = seq(0, 7000, 2000)) +
  labs(x = "Latitude", y = "Elevation (m)",
       subtitle = paste0("Grey: crest; black: lowest crossing per 0.25°; dots: gateways; bands: gaps < ", P$montane_gap_m[[1]], " m")) +
  th + theme(plot.subtitle = element_text(size = 6, colour = "grey30"))
tr <- cols[tropical & depth >= P$depth_m[[1]] & !terminal]
tr[, cls := factor(ifelse(gateway, "Gateway", "Col"), c("Gateway", "Col"))]
ex <- tr[gateway == TRUE][unique(c(order(elev_cross)[1:5], order(-iso_km)[1:5]))]       # 5 lowest + 5 most isolated
ex[, txt := ifelse(is.na(lab), latlab(lat, 1), lab)]
# Two shaded bands: the five lowest and the five most isolated gateways; their crossing holds the gateways that meet both criteria
gw <- tr[gateway == TRUE]
y_cut <- mean(sort(gw$elev_cross)[5:6]); x_cut <- exp(mean(log(sort(gw$iso_km, decreasing = TRUE)[5:6])))
band <- data.frame(xmin = c(0, x_cut), xmax = c(Inf, Inf), ymin = c(-Inf, -Inf), ymax = c(y_cut, Inf),
                   crit = factor(c("Five lowest gateways", "Five most isolated gateways")))
f1c <- ggplot(tr, aes(iso_km, elev_cross)) +
  geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = crit), inherit.aes = FALSE, alpha = 0.22) +
  scale_fill_manual(values = c(oi[["yell"]], oi[["sky"]])) +
  geom_point(aes(colour = cls, shape = cls, size = cls), alpha = 0.8) +
  geom_text_repel(data = ex, aes(label = txt), size = 1.9, segment.size = 0.2, min.segment.length = 0, seed = 1, max.overlaps = Inf) +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) +
  annotation_logticks(sides = "b", linewidth = 0.2, short = unit(0.4, "mm"), mid = unit(0.8, "mm"), long = unit(1.2, "mm")) +
  scale_colour_manual(values = c(oi[["blue"]], "grey60")) + scale_shape_manual(values = c(16, 1)) +
  scale_size_manual(values = c(1.2, 0.6)) +
  labs(x = "Isolation along the divide (km, log scale)", y = "W–E crossing height (m)",
       subtitle = sprintf("Tropical Andes, cols with depth ≥ 100 m (n = %d)", nrow(tr))) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  th + theme(plot.subtitle = element_text(size = 6, colour = "grey30"))
save2(f1a + (f1b / f1c) + plot_layout(widths = c(1, 1.35)) + plot_annotation(tag_levels = "a"), "Fig1_barrier_gateways", 170, 150)

# =====================================================================================================================
# Fig. 2 - Janzen along the Andes (H2) and the hydric dimension (H3)
# =====================================================================================================================
jm <- readRDS(file.path(tb, "M2_modelo_janzen.rds"))$main$model
jd <- jm$model                                                                    # the gateways used in the fit
# fixed-effect curves only: the pointwise SE is dominated by intercept-smooth confounding (CIs of the terms in Table 2)
pr <- function(nd) { nd$d_km <- median(jd$d_km); p <- predict(jm, nd, exclude = "s(d_km)", se.fit = TRUE)
  cbind(nd, fit = p$fit, lo = p$fit - 1.96 * p$se.fit, hi = p$fit + 1.96 * p$se.fit) }
gw3 <- cw[gateway == TRUE][, `:=`(abs_lat = abs(lat), elev_km = elev_cross / 1000)]
# curves only over the observed support (gateways within 0.5 km of the height, or of the zone) and up to 1, the index maximum
nd_a <- rbindlist(lapply(c(1.5, 3.5), function(h) { r <- range(jd$abs_lat[abs(jd$elev_km - h) <= 0.5])
  pr(data.frame(abs_lat = seq(r[1], r[2], 0.25), elev_km = h)) }))[fit <= 1]
f2a <- ggplot() +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.2, linetype = "13") +
  geom_line(data = nd_a, aes(abs_lat, fit, linetype = factor(elev_km)), linewidth = 0.4) +
  geom_point(data = gw3, aes(abs_lat, TO_ann, fill = elev_cross), shape = 21, size = 1.4, stroke = 0.2) +
  scale_elev("fill", "Crossing height (m)") +
  scale_linetype_manual("GAM fit at", values = c("22", "solid"), labels = c("1.5 km", "3.5 km")) +
  labs(x = "|Latitude| (°)", y = "Thermal overlap, crest vs flank") + th +
  theme(legend.title = element_text(size = 6), legend.key.width = unit(5, "mm"), legend.box = "vertical", legend.box.just = "left",
        legend.spacing.y = unit(0.5, "mm"))
zl <- c("Tropical (5°)", "Temperate (30°)")
gw3[, zone := factor(ifelse(tropical, zl[1], zl[2]), zl)]
nd_b <- rbindlist(lapply(c(5, 30), function(a) { r <- range(gw3[tropical == (a == 5), elev_km])
  pr(data.frame(elev_km = seq(r[1], r[2], 0.025), abs_lat = a)) }))[fit <= 1]
nd_b[, zone := factor(ifelse(abs_lat == 5, zl[1], zl[2]), zl)]
f2b <- ggplot() +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.2, linetype = "13") +
  geom_line(data = nd_b, aes(elev_km * 1000, fit, colour = zone, linetype = zone), linewidth = 0.5) +
  geom_point(data = gw3, aes(elev_cross, TO_ann, colour = zone, shape = zone), size = 1) +
  scale_colour_manual(values = c(oi[["blue"]], oi[["sky"]])) +
  scale_shape_manual(values = c(16, 2)) + scale_linetype_manual(values = c("solid", "22")) +
  labs(x = "W–E crossing height (m)", y = "Thermal overlap, crest vs flank") + th
# (c) H3 typology: median cuts fixed before classifying (BITACORA 2026-11-03), types from 17_janzen_M2.R
h3 <- fread(file.path(tb, "M2_H3_tipologia.csv"))
h3 <- cols[, .(col_id, lab)][h3, on = "col_id"][, txt := ifelse(is.na(lab), latlab(lat, 1), lab)]
h3[, txt := paste0(toupper(substr(txt, 1, 1)), substring(txt, 2))]                   # "la Viuda" after dropping "Abra"
tl3 <- c(open = "Open", hydric_only = "Closed in moisture only", thermal_only = "Closed in temperature only", closed = "Closed")   # manuscript terms
h3[, type := factor(tl3[h3_type], tl3)]
med <- h3[, .(TO = median(TO_ann), AI = median(AI_gap))]
h3[, ny := fcase(txt == "Meta-Huila", -0.7, txt == "Portete de Tarqui", 0.6, default = 0)]   # ponytail: hand-nudge the one clash
f2c <- ggplot(h3, aes(TO_ann, AI_gap)) +
  geom_vline(xintercept = med$TO, linetype = "22", linewidth = 0.25, colour = "grey50") +
  geom_hline(yintercept = med$AI, linetype = "22", linewidth = 0.25, colour = "grey50") +
  geom_point(aes(colour = type, shape = type), size = 1.6) +
  geom_text_repel(aes(label = txt), size = 1.8, segment.size = 0.2, min.segment.length = 0, seed = 1, max.overlaps = Inf,
                  box.padding = 0.4, point.padding = 0.15, force = 3, nudge_y = h3$ny,
                  colour = "grey25") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = 2, colour = "grey30", label = "open in temperature, closed in moisture") +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.05, vjust = -1, size = 2, colour = "grey30", label = "closed in temperature, open in moisture") +
  scale_colour_manual(values = unname(oi[c("blue", "green", "yell", "pink")]), drop = FALSE) +
  scale_shape_manual(values = c(16, 17, 15, 18), drop = FALSE) +
  scale_x_continuous(breaks = seq(-0.25, 1, 0.25)) +
  labs(x = "Thermal overlap, crest vs flank (higher = more permeable)", y = "Aridity gap, |log(AI crest / AI flank)|",
       subtitle = sprintf("%d tropical gateways; dashed lines: medians (type cuts)", nrow(h3))) +
  th + theme(plot.subtitle = element_text(size = 6, colour = "grey30"))
save2((f2a | f2b) / f2c + plot_layout(heights = c(1, 1.15)) +
        plot_annotation(tag_levels = "a", caption = "Non-terminal gateways, 1000 m flank, worse side. Lines in a–b: GAM fixed effects (divide smooth excluded),\ndrawn over the observed range and up to the index maximum of 1 (dotted).",
                        theme = theme(plot.caption = element_text(size = 6, colour = "grey30"))),
      "Fig2_janzen_moisture", 170, 150)

# =====================================================================================================================
# Fig. 3 - future permeability (H4)
# =====================================================================================================================
ens <- fread(file.path(tb, "M3_futuro_ensamble.csv"))[col_id %in% tg$col_id]
dts <- fread(file.path(tb, "M3_dTstar.csv"))[col_id %in% tg$col_id]
now <- dts[, .(dTstar = max(dTstar)), by = col_id][tg[, .(col_id, lat, elev_cross)], on = "col_id"]
now <- cols[, .(col_id, lab)][now, on = "col_id"][, txt := ifelse(is.na(lab), latlab(lat, 1), lab)]
plab <- function(p) sub("-", "–", p)                                              # period labels with en dash
ssp_lab <- function(s) factor(paste0("SSP", substr(s, 1, 1), "-", substr(s, 2, 2), ".", substr(s, 3, 3)))
e3 <- tg[, .(col_id, lat)][ens, on = "col_id"][, SSP := ssp_lab(ssp)]
# (a) change vs present, every tropical gateway; filled = robust sign (>= 80 % of GCMs, i.e. 7 of 8)
ch <- melt(e3, id.vars = c("col_id", "lat", "SSP", "period"),
           measure.vars = list(c("dTO.mean", "dPg.mean"), c("dTO.min", "dPg.min"), c("dTO.max", "dPg.max"), c("dTO.agree", "dPg.agree")),
           value.name = c("mean", "min", "max", "agree"))
ch[, var := factor(variable, 1:2, c("Δ thermal overlap", "Δ absolute precipitation gap"))]
ch[, period := plab(period)]
ch[, robust := factor(ifelse(agree >= 0.8, "Sign agreed by ≥ 7 of 8 GCMs", "Weaker agreement"),
                      c("Sign agreed by ≥ 7 of 8 GCMs", "Weaker agreement"))]
dg <- position_dodge(width = 0.9)
f3a <- ggplot(ch, aes(lat, mean, colour = SSP, group = SSP)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.2) +
  geom_linerange(aes(ymin = min, ymax = max), linewidth = 0.2, position = dg, show.legend = FALSE) +
  geom_point(aes(shape = robust), size = 0.9, stroke = 0.3, position = dg) +
  facet_grid(var ~ period, scales = "free_y", switch = "y", labeller = labeller(var = label_wrap_gen(18))) +
  scale_x_reverse(breaks = seq(5, -20, -5), labels = latlab) +
  scale_colour_manual(values = unname(oi[c("sky", "blue", "verm")])) + scale_shape_manual(values = c(16, 1)) +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  labs(x = "Latitude", y = NULL, subtitle = sprintf("%d tropical gateways; point: mean of 8 GCMs; bar: GCM range; change vs present (1970–2000)", nrow(tg))) +
  th + theme(strip.placement = "outside", plot.subtitle = element_text(size = 6, colour = "grey30"))
e90 <- e3[period == "2081-2100"][, noan := ifelse(p_over_T > 0.5, "No analogue (> half of GCMs)", "Analogue in present lowlands")]
# (b) against crossing height, the relation the text reports (r = 0.95); the five lowest gateways labelled
f3b1 <- ggplot(now, aes(elev_cross, dTstar)) + geom_point(size = 1.2, colour = "grey30") +
  geom_text_repel(data = now[order(elev_cross)][1:5], aes(label = txt), size = 1.9, segment.size = 0.2,
                  min.segment.length = 0, seed = 1, nudge_x = 400, colour = "grey25") +
  labs(x = "W–E crossing height (m)", y = "ΔT* (°C)",
       subtitle = sprintf("r = %.2f", cor(now$elev_cross, now$dTstar))) +
  th + theme(plot.subtitle = element_text(size = 6, colour = "grey30"))
f3b2 <- ggplot(e90[!is.na(frac_dTstar.mean)], aes(lat, frac_dTstar.mean, colour = SSP, shape = noan)) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.2, linetype = "22") +
  geom_linerange(aes(ymin = frac_dTstar.min, ymax = frac_dTstar.max), linewidth = 0.2, show.legend = FALSE) +
  geom_point(size = 1.2) +
  annotate("text", x = e90[frac_dTstar.max > 1.6, lat][1] - 0.6, y = max(e90$frac_dTstar.max, na.rm = TRUE), hjust = 0,
           size = 1.9, colour = "grey30", label = "Táchira") +
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5), labels = c("0.1", "0.2", "0.5", "1", "2", "5")) +   # Táchira (dT* = 1 C) spans 1.9-5.4
  scale_x_reverse(breaks = seq(10, -20, -10), labels = latlab) +
  scale_colour_manual(values = unname(oi[c("sky", "blue", "verm")])) + scale_shape_manual(values = c(16, 4)) +
  guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(nrow = 1, order = 2)) +
  labs(x = "Latitude", y = "Fraction of ΔT* reached, 2081–2100 (log scale)") + th +
  theme(legend.box = "vertical", legend.box.just = "left")
save2(f3a / (f3b1 | f3b2) + plot_layout(heights = c(1.25, 1)) + plot_annotation(tag_levels = "a"),
      "Fig3_future_permeability", 170, 170)

# =====================================================================================================================
# Fig. 4 - conservation: montane corridor fragmentation, E-W gateways opening, protection and pressure (plan v2.2)
# =====================================================================================================================
bn <- fread(file.path(int, "M3_corredor_bins_proteccion.csv"))
bn[, st := factor(fifelse(gap_now, "Montane gap today", fifelse(loss_robust_370_2090, "Lost by 2081–2100 (≥ 7/8 GCMs)", "Montane corridor kept")),
                  c("Montane corridor kept", "Lost by 2081–2100 (≥ 7/8 GCMs)", "Montane gap today"))]
pa4 <- sf::st_simplify(sf::st_read(file.path(int, "M3_wdpa_andes.gpkg"), quiet = TRUE), dTolerance = 0.01)
pr4 <- fread(file.path(tb, "M3_prioridad_puertas.csv"), encoding = "UTF-8")
pr4[, lab := fifelse(is.na(name), latlab(lat, 1), sub("^Abra (de |el )?", "", name))][abs(lat - 7.686) < 0.01, lab := "Táchira"]
pr4[, cls := factor(fifelse(critical, fifelse(high_pressure, "Critical, unprotected, high pressure", "Critical, unprotected"), "Other gateway"),
                    c("Critical, unprotected, high pressure", "Critical, unprotected", "Other gateway"))]
f4a <- ggplot() +
  geom_sf(data = sa, fill = "grey92", colour = NA) +
  geom_sf(data = pa4, fill = oi[["green"]], colour = NA, alpha = 0.35) +
  geom_path(data = bn, aes(lon, lat, colour = st, group = 1), linewidth = 0.4) +
  geom_point(data = bn[st != "Montane corridor kept"], aes(lon, lat, colour = st), size = 0.9, show.legend = FALSE) +   # short stretches visible
  geom_point(data = pr4, aes(lon, lat, shape = cls, fill = cls), size = 1.6, stroke = 0.3) +
  geom_text_repel(data = pr4[critical == TRUE & !grepl("Cajanuma", lab)], aes(lon, lat, label = lab), size = 1.8, nudge_x = -3.8, direction = "y", hjust = 1,
                  segment.size = 0.2, min.segment.length = 0, seed = 1, max.overlaps = Inf) +
  geom_text_repel(data = pr4[grepl("Cajanuma", lab)], aes(lon, lat, label = lab), size = 1.8, nudge_x = 3, direction = "y", hjust = 0,   # no room to the west
                  segment.size = 0.2, min.segment.length = 0, seed = 1, max.overlaps = Inf) +
  scale_colour_manual(values = c(oi[["blue"]], oi[["verm"]], "grey60"), name = NULL) +
  scale_shape_manual(values = c(24, 24, 21), name = NULL) + scale_fill_manual(values = c("black", "white", "grey60"), name = NULL) +
  coord_sf(xlim = c(-84, -63), ylim = c(-23.5, 11.5), expand = FALSE) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(linewidth = 1.2)), shape = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) +
  labs(x = NULL, y = NULL, subtitle = "Green: protected areas (WDPA)") + th +
  theme(legend.position = "bottom", legend.box = "vertical", axis.line = element_blank(), plot.subtitle = element_text(size = 6, colour = "grey30"))
ce <- fread(file.path(tb, "M3_corredor_montano_ensamble.csv"))[zone == "tropical" & threshold_m == P$montane_gap_m[[1]]]
c0 <- fread(file.path(tb, "M3_corredor_montano.csv"))[zone == "tropical" & threshold_m == P$montane_gap_m[[1]] & gcm == "present"]
ce <- rbind(ce[, .(ssp, period, m = gap_km.mean, lo = gap_km.min, hi = gap_km.max)],
            rbindlist(lapply(c(126, 370, 585), function(s) c0[, .(ssp = s, period = "1970-2000", m = gap_km, lo = gap_km, hi = gap_km)])))
# x at period midpoints, so the 65 years between baseline and 2041-2060 are not drawn like a 20-year step
pm <- c("1970-2000" = 1985, "2041-2060" = 2050, "2061-2080" = 2070, "2081-2100" = 2090)
ce[, `:=`(SSP = ssp_lab(ssp), t = pm[period])]
f4b <- ggplot(ce, aes(t, m, colour = SSP, fill = SSP)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.5) + geom_point(size = 1) +
  scale_x_continuous(breaks = pm, labels = plab(names(pm))) +
  scale_colour_manual(values = unname(oi[c("sky", "blue", "verm")])) + scale_fill_manual(values = unname(oi[c("sky", "blue", "verm")])) +
  labs(x = NULL, y = "Tropical divide in montane gaps (km)", subtitle = sprintf("Gap: crossing below today's %d m isotherm; line: 8-GCM mean, band: range", P$montane_gap_m[[1]])) +
  th + theme(plot.subtitle = element_text(size = 6, colour = "grey30"))
cs <- fread(file.path(tb, "M3_corredor_proteccion.csv"))[stretch != "montane corridor today"]
cs <- melt(cs, id.vars = "stretch", measure.vars = c("pct_in_PA", "pct_high_pressure"))
cs[, `:=`(stretch = factor(stretch, c("montane corridor kept", "robust loss by 2081-2100 (SSP3-7.0)"), c("Kept", "Lost by 2081–2100")),
          variable = factor(variable, c("pct_in_PA", "pct_high_pressure"), c("Inside protected areas", "High human modification
(gHM ≥ 0.4)")))]
f4c <- ggplot(cs, aes(stretch, value, fill = stretch)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", value)), vjust = -0.4, size = 2) +
  facet_wrap(~ variable) + scale_fill_manual(values = c(oi[["blue"]], oi[["verm"]]), guide = "none") +
  scale_y_continuous(limits = c(0, 40), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "% of montane corridor (tropical divide)") + th
save2(f4a + (f4b / f4c) + plot_layout(widths = c(1, 1.1)) + plot_annotation(tag_levels = "a"), "Fig4_conservation", 170, 170)

message("figures in ", fd)
