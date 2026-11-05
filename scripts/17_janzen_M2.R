# 17_janzen_M2.R - M2.6: Janzen test (H2) and M2.7: continental H3 typology. Huancabamba percentiles: L2 scripts/17b_percentiles_huancabamba.R.
# Col-level permeability = the lower of the two sides (the harder climb), flank threshold 1000 m.
# H2: TO_ann ~ elev_km * abs_lat + s(divide km), Gaussian GAM (TO can be < 0, so no beta family), unit = gateway
# (terminal channel crossings excluded); Moran's I of residuals with weights 1 / divide distance.
# Sensitivity: all cols with depth >= 100 m; TO_mon; W and E sides separately; flank 500 / 1500 m.
# Rscript scripts/17_janzen_M2.R [--force]
suppressPackageStartupMessages({ library(mgcv); library(ggplot2) })
cfg <- yaml::read_yaml("config.yml"); P <- cfg$params; o <- cfg$paths$outputs
force <- "--force" %in% commandArgs(TRUE)
out <- file.path(o, "tables", "M2_modelo_janzen.rds")
if (file.exists(out) && !force) { message("exists, skip"); quit(save = "no") }

m <- read.csv(file.path(o, "tables", "M2_metricas_clima_collados.csv"))
m <- m[!m$terminal, ]
H <- P$huancabamba
m$huanca <- m$lat >= H$lat[1] & m$lat <= H$lat[2] & m$elev_cross < H$max_elev_m & m$depth >= H$min_depth_m

# One row per col and threshold: the worst side for barrier metrics
worst <- function(d) {
  a <- aggregate(cbind(TO_ann, TO_mon) ~ col_id + flank_m, d, min, na.action = na.pass)
  b <- aggregate(maha_pct ~ col_id + flank_m, d, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE), na.action = na.pass)
  # Moisture gaps: the side with the larger ABSOLUTE gap, as overlap is symmetric: a crest much drier than its flank is a
  # barrier too. AI_gap / CLD_gap = magnitude; *_signed and *_side keep direction and side.
  big <- function(v) do.call(rbind, lapply(split(d, list(d$col_id, d$flank_m), drop = TRUE), function(x) {
    i <- which.max(abs(x[[v]]))[1]
    setNames(data.frame(x$col_id[1], x$flank_m[1], abs(x[[v]][i]), x[[v]][i], x$side_lab[i]),
             c("col_id", "flank_m", v, paste0(v, "_signed"), paste0(v, "_side")))
  }))
  b <- merge(merge(b, big("AI_gap")), big("CLD_gap"))
  k <- unique(d[, c("col_id", "name", "lon", "lat", "elev_cross", "depth", "d_km", "tropical", "gateway", "lara", "huanca")])
  merge(k, merge(a, b))
}
cw <- worst(m)
cw$elev_km <- cw$elev_cross / 1000; cw$abs_lat <- abs(cw$lat)

fit <- function(d, y = "TO_ann") {
  d <- d[is.finite(d[[y]]), ]
  g <- gam(as.formula(paste(y, "~ elev_km * abs_lat + s(d_km, k =", min(10, nrow(d) %/% 5), ")")), data = d, method = "REML")
  w <- 1 / as.matrix(dist(d$d_km)); diag(w) <- 0; w[!is.finite(w)] <- 0
  mi <- ape::Moran.I(residuals(g), w)
  ct <- summary(g)$p.table
  list(model = g, n = nrow(d), coef = ct, ci = cbind(lo = ct[, 1] - 1.96 * ct[, 2], hi = ct[, 1] + 1.96 * ct[, 2]),
       moran = unlist(mi[c("observed", "p.value")]))
}
def <- cw[cw$flank_m == P$flank_m[1], ]
side1 <- m[m$flank_m == P$flank_m[1], ]; side1$elev_km <- side1$elev_cross / 1000; side1$abs_lat <- abs(side1$lat)
mods <- list(
  main      = fit(def[def$gateway, ]),
  all_cols  = fit(def[def$depth >= P$depth_m[1], ]),
  monthly   = fit(def[def$gateway, ], "TO_mon"),
  side_W    = fit(side1[side1$gateway & side1$side_lab == "W", ]),
  side_E    = fit(side1[side1$gateway & side1$side_lab == "E", ]),
  flank_500 = fit(cw[cw$gateway & cw$flank_m == 500, ]),
  flank_1500 = fit(cw[cw$gateway & cw$flank_m == 1500, ]))
saveRDS(mods, out)
smry <- do.call(rbind, lapply(names(mods), function(k) {
  x <- mods[[k]]; r <- "elev_km:abs_lat"
  data.frame(model = k, n = x$n, b_elev = x$coef["elev_km", 1], b_int = x$coef[r, 1],
             int_lo = x$ci[r, "lo"], int_hi = x$ci[r, "hi"], int_p = x$coef[r, 4],
             dev_expl = summary(x$model)$dev.expl, moran_I = x$moran["observed"], moran_p = x$moran["p.value"])
}))
write.csv(smry, file.path(o, "tables", "M2_janzen_resumen.csv"), row.names = FALSE)
print(smry, digits = 3)

# QC figures
def$grp <- ifelse(def$name %in% "Abra de Porculla", "Porculla", ifelse(def$huanca, "Huancabamba system", "other"))
p1 <- ggplot(def[def$gateway, ], aes(elev_cross, TO_ann, colour = abs_lat)) +
  geom_point(size = 2) + geom_point(data = def[def$gateway & def$grp == "Porculla", ], shape = 21, size = 5, colour = "red") +
  scale_colour_viridis_c("|lat|") + labs(x = "Gateway crossing height (m)", y = "Thermal overlap, crest vs worst flank") +
  theme_bw()
ggsave(file.path(o, "figs_qc", "M2_janzen_gateways.png"), p1, width = 7, height = 5, dpi = 200)

# ---- H3 continental (M2.7): rules fixed in BITACORA 2026-11-03 before classifying ----
# Tropical gateways, worst side. Thermal open: TO_ann >= median; hydric closed: absolute AI_gap >= median (crest wetter or drier than its flank).
# Sensitivity: terciles (middle third = intermediate), CLD_gap as hydric axis, all gateways.
gw <- def[def$gateway, ]; tg <- gw[gw$tropical, ]
typ <- function(d, h = "AI_gap", q = 0.5) {
  if (q == 0.5) { to <- d$TO_ann >= median(d$TO_ann, na.rm = TRUE); hc <- d[[h]] >= median(d[[h]], na.rm = TRUE) } else {
    cut3 <- function(x, hi_is) { b <- quantile(x, c(1, 2) / 3, na.rm = TRUE); ifelse(x >= b[2], hi_is, ifelse(x <= b[1], !hi_is, NA)) }
    to <- cut3(d$TO_ann, TRUE); hc <- cut3(d[[h]], TRUE)
  }
  ifelse(is.na(to) | is.na(hc), "intermediate",
         ifelse(to, ifelse(hc, "hydric_only", "open"), ifelse(hc, "closed", "thermal_only")))
}
tg$h3_type <- typ(tg); tg$h3_type_terc <- typ(tg, q = 1 / 3); tg$h3_type_cld <- typ(tg, "CLD_gap")
gw$h3_type_all <- typ(gw)
tg$h3_type_all <- gw$h3_type_all[match(tg$col_id, gw$col_id)]
write.csv(tg[order(tg$lat, decreasing = TRUE), c("col_id", "name", "lon", "lat", "elev_cross", "depth", "d_km", "TO_ann", "AI_gap",
             "AI_gap_signed", "AI_gap_side", "CLD_gap", "CLD_gap_signed", "CLD_gap_side", "maha_pct", "h3_type", "h3_type_terc", "h3_type_cld", "h3_type_all")],
          file.path(o, "tables", "M2_H3_tipologia.csv"), row.names = FALSE)
sp <- function(a, b, d = tg) { r <- cor.test(d[[a]], d[[b]], method = "spearman", exact = FALSE)
  data.frame(test = "spearman", x = a, y = b, n = sum(complete.cases(d[, c(a, b)])), est = unname(r$estimate), lo = NA, hi = NA,
             p = r$p.value, dev_expl = NA, moran_I = NA, moran_p = NA, dev_smooth = NA) }
h3fit <- function(y, d = tg) {
  d <- d[complete.cases(d[, c(y, "elev_km", "TO_ann")]), ]
  g <- gam(as.formula(paste(y, "~ elev_km + TO_ann + s(d_km, k =", min(10, nrow(d) %/% 5), ")")), data = d, method = "REML")
  w <- 1 / as.matrix(dist(d$d_km)); diag(w) <- 0; w[!is.finite(w)] <- 0; mi <- ape::Moran.I(residuals(g), w)
  ct <- summary(g)$p.table[-1, , drop = FALSE]
  g0 <- gam(as.formula(paste(y, "~ s(d_km, k =", min(10, nrow(d) %/% 5), ")")), data = d, method = "REML")   # the smooth alone
  data.frame(test = "gam", x = rownames(ct), y = y, n = nrow(d), est = ct[, 1], lo = ct[, 1] - 1.96 * ct[, 2],
             hi = ct[, 1] + 1.96 * ct[, 2], p = ct[, 4], dev_expl = summary(g)$dev.expl, moran_I = mi$observed, moran_p = mi$p.value,
             dev_smooth = summary(g0)$dev.expl)
}
h3 <- rbind(sp("TO_ann", "AI_gap"), sp("TO_ann", "CLD_gap"), sp("AI_gap", "CLD_gap"), sp("elev_cross", "AI_gap"),
            sp("elev_cross", "CLD_gap"), h3fit("AI_gap"), h3fit("CLD_gap"), transform(sp("TO_ann", "AI_gap", gw), test = "spearman_all"))
write.csv(h3, file.path(o, "tables", "M2_H3_resumen.csv"), row.names = FALSE)
print(h3, digits = 3); print(table(tg$h3_type)); print(table(tg$h3_type_terc)); print(table(tg$h3_type_cld))
