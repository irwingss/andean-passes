# 33_suplementos.R - Article 1 Supporting Information (plan v2.4 section 8.3), English, from the M1-M3 outputs only.
# S1 census algorithm and validation; S2 sensitivity of the census and metrics; S3 sensitivity of the H3 typology;
# S4 gateway atlas (one page per tropical gateway); S5 GCM selection and M3 results by GCM.
# Writes outputs/supp/L1_Supplementary_Tables.xlsx (sheet "Contents" holds every title and note), Figs S1-S5 and
# outputs/supp/atlas/*.png. The Word file is built from these by scripts/34_build_supp.py.
# Rscript scripts/33_suplementos.R
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(data.table); library(terra) })
suppressMessages(sf::sf_use_s2(FALSE))
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; raw <- cfg$paths$data_dir; tb <- file.path(cfg$paths$outputs, "tables")
P <- cfg$params; f1000 <- P$flank_m[[1]]
er <- cfg$style$elev_ramp                                         # one terrain ramp for every elevation colour (config.yml)
scale_elev <- function(aesthetics, name) scale_colour_gradientn(name, colours = er$col, values = scales::rescale(er$m),
  limits = range(er$m), oob = scales::squish, aesthetics = aesthetics)
sd <- file.path(cfg$paths$outputs, "supp"); dir.create(file.path(sd, "atlas"), recursive = TRUE, showWarnings = FALSE)
oi <- c(verm = "#D55E00", blue = "#0072B2", sky = "#56B4E9", green = "#009E73", yell = "#E69F00", pink = "#CC79A7", grey = "grey55")
th <- theme_classic(base_size = 7) +
  theme(legend.position = "top", legend.title = element_blank(), legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), plot.tag = element_text(size = 9, face = "bold"),
        strip.background = element_blank(), strip.text = element_text(size = 7, face = "bold"))
save2 <- function(p, file, w, h, dpi = 600) {
  ggsave(paste0(file, ".png"), p, width = w, height = h, units = "mm", dpi = dpi, device = ragg::agg_png, bg = "white")
  if (dpi == 600) ggsave(paste0(file, ".pdf"), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
}
latlab <- function(y, d = 1) ifelse(y > 0, sprintf(paste0("%.", d, "f°N"), y), ifelse(y < 0, sprintf(paste0("%.", d, "f°S"), -y), "0°"))
ssp_lab <- c("126" = "SSP1-2.6", "370" = "SSP3-7.0", "585" = "SSP5-8.5")
type_lab <- c(open = "Open", thermal_only = "Closed in temperature only", hydric_only = "Closed in moisture only",
              closed = "Closed", intermediate = "Intermediate")
r <- function(x, d) round(x, d)

# ---- shared data: the 30 tropical gateways, names and countries as in Table 3 ----
gw <- fread(file.path(tb, "Table3_conservation_priority.csv"), encoding = "UTF-8")[, .(col_id, name, country)]
cols <- fread("outputs/database/APD_cols.csv", encoding = "UTF-8")
m2 <- fread(file.path(tb, "M2_metricas_clima_collados.csv"))
fut <- fread(file.path(tb, "M3_futuro.csv"), colClasses = list(character = "ssp"))[gateway & tropical]
ens <- fread(file.path(tb, "M3_futuro_ensamble.csv"), colClasses = list(character = "ssp"))
typ <- fread(file.path(tb, "M2_H3_tipologia.csv"))
stopifnot(nrow(gw) == 30, all(gw$col_id %in% typ$col_id))

sheets <- list(); meta <- list()
add <- function(id, title, d, note = "") { sheets[[id]] <<- d; meta[[id]] <<- data.table(sheet = id, title = title, note = note) }

# =====================================================================================================================
# S1 - validation of crossing heights
# =====================================================================================================================
v <- fread(file.path(tb, "M1_validacion.csv"), encoding = "UTF-8")
vs <- fread(file.path(tb, "M1_validacion_resumen.csv"))
# Exact source of every reference elevation, each checked against the original table or page
bib <- fread("data/refs/M1_reference_bibliography.csv", encoding = "UTF-8")
srcs <- fread("data/refs/M1_reference_sources.csv", encoding = "UTF-8", colClasses = "character")
stopifnot(setequal(srcs$name, v$name))
lab <- function(id) sapply(strsplit(id, "|", fixed = TRUE), \(k) paste(bib$label[match(k, bib$ref_id)], collapse = "; "))
srcs[, `:=`(ref = lab(ref_id), alt_ref = fifelse(alt_ref_id == "NONE", "No publication traced", fifelse(alt_ref_id == "", "No alternative source", lab(alt_ref_id))))]
stopifnot(!anyNA(srcs$ref), !anyNA(srcs$alt_ref))
v <- merge(v, srcs, by = "name", sort = FALSE)
add("S1.1", "Table S1.1. Validation of crossing heights against 28 documented passes.",
    v[order(-lat), .(Pass = name, Country = country, Latitude = r(lat, 4), Longitude = r(lon, 4),
                     `Reference elevation (m)` = elev_ref_m, `Alternative elevation (m)` = elev_alt_m, Source = ref,
                     `Distance to divide (km)` = r(dist_divide_km, 2), `Distance to nearest col (km)` = r(dist_col_km, 2),
                     `Matched (col within 5 km)` = matched, `Col depth (m)` = r(col_depth, 0),
                     `GLO-30 crossing height (m)` = r(elev_glo30, 0), `Error (m)` = r(err_glo30, 0),
                     `Error vs alternative source (m)` = r(err_glo30_alt, 0))],
    sprintf(paste("Error is the GLO-30 crossing height of the matched col minus the reference elevation; it is computed for every",
                  "pass but used only where a col lies within 5 km (Matched). Acceptance criterion, set before matching: median",
                  "absolute error of 50 m or less and median bias within 30 m. Result: %d of %d passes matched, median absolute",
                  "error %.0f m, median bias %.0f m, mean bias %.0f m; with the alternative source, median absolute error %.0f m."),
            vs$n_matched, vs$n_refs, vs$median_abs_err_glo30, vs$median_bias_glo30, vs$mean_bias_glo30, vs$median_abs_err_glo30_altsrc))
vm <- v[matched == TRUE]
fs1 <- ggplot(vm, aes(elev_ref_m, elev_glo30)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_point(shape = 21, size = 1.6, fill = oi[["blue"]], stroke = 0.2) +
  ggrepel::geom_text_repel(data = vm[abs(err_glo30) > 50], aes(label = name), size = 2, seed = 1, min.segment.length = 0) +
  coord_equal() + labs(x = "Published elevation (m)", y = "GLO-30 crossing height of the matched col (m)") + th
save2(fs1, file.path(sd, "FigS1_validation"), 85, 85)
add("S1.2", "Table S1.2. Published source of the reference elevation of the 28 documented passes, and whether each pass appears in Figure S1.2.",
    v[order(-lat), .(Pass = name, `Shown in Fig. S1.2` = matched, `Distance to nearest col (km)` = r(dist_col_km, 2),
                     `Reference elevation (m)` = elev_ref_m, Reference = ref, `Location in the source` = locator,
                     `Alternative elevation (m)` = elev_alt_m, `Reference of the alternative elevation` = alt_ref,
                     `Location of the alternative elevation` = alt_locator)],
    paste("Full references are in Table S1.3. Passes not shown in Fig. S1.2 have no col within 5 km, so their published elevation was not",
          "compared. Each value was checked against the cited table or page. Pass coordinates come from the same INEI tables, except",
          "Anticona, Antajirca and the Chile-Argentina passes, whose coordinates are OpenStreetMap nodes. The alternative elevation of",
          "Anticona, 4818 m, is a commonly cited value that we could not trace to a publication. Weigend (2002) gives 2145 m for",
          "Porculla. Cristo Redentor and Túnel de Galera were excluded before matching, because their published elevations refer to a",
          "border complex and a tunnel, not to a col."))
add("S1.3", "Table S1.3. Full references of the sources in Table S1.2.",
    bib[, .(Reference = label, `Full reference` = citation, URL = url, `Source credited in the cited table` = credited_in_source)],
    "The pages of the Unidad de Pasos Fronterizos carry no publication date. All URLs were last checked on 3 November 2026.")

# =====================================================================================================================
# S2 - sensitivity of the census and of the metrics
# =====================================================================================================================
s21 <- fread(file.path(tb, "M1_sensibilidad_umbral_ventana.csv"))
add("S2.1", "Table S2.1. Sensitivity of the census to the depth threshold and to the gateway window.",
    s21[order(win, depth), .(`Depth threshold (m)` = depth, `Gateway window (km)` = win, `Cols with depth above threshold` = n_cols,
                             Gateways = n_gateways, `Tropical gateways` = n_gateways_trop,
                             `Lowest non-terminal tropical gateway, latitude` = r(lowest_trop_lat, 3),
                             `Lowest non-terminal tropical gateway, crossing height (m)` = r(lowest_trop_elev, 0))],
    "Default: depth 100 m and window 100 km (both directions along the divide). Terminal cols and gateways are included in the counts.")

jz <- fread(file.path(tb, "M2_janzen_resumen.csv"))
env <- regmatches(readLines(file.path(tb, "M2_janzen_mecanismo_modelos.txt")), regexpr("^TO_env.*", readLines(file.path(tb, "M2_janzen_mecanismo_modelos.txt"))))
en <- as.numeric(regmatches(env, gregexpr("[-+]?[0-9.]+e?-?[0-9]*", env))[[1]])          # est, lo, hi, p, dev
jz <- rbind(jz, data.table(model = "envelope", n = 61, b_elev = NA, b_int = en[1], int_lo = en[2], int_hi = en[3], int_p = en[4],
                           dev_expl = en[5], moran_I = NA, moran_p = NA))
jlab <- c(main = "Main: annual index TOann, worse side, flank 1000 m", side_W = "Western side only", side_E = "Eastern side only",
          flank_500 = "Flank at 500 m", flank_1500 = "Flank at 1500 m", all_cols = "All cols with depth of at least 100 m",
          monthly = "Janzen's monthly index TOmon (0-12 scale)", envelope = "TOenv, monthly crest vs annual flank envelope (0-12 scale)")
jz[, lab := factor(jlab[model], levels = jlab)]
add("S2.2", "Table S2.2. Sensitivity of the test of Janzen's prediction (GAM of thermal overlap on crossing height, absolute latitude and their interaction, with a smooth of position along the divide).",
    jz[order(lab), .(Model = as.character(lab), n, `Slope of height at the equator (per km)` = r(b_elev, 3),
                     `Height x latitude interaction (per km per degree)` = r(b_int, 4), `95% CI low` = r(int_lo, 4),
                     `95% CI high` = r(int_hi, 4), p = signif(int_p, 2), `Deviance explained` = r(dev_expl, 3),
                     `Moran's I of residuals` = r(moran_I, 3), `Moran p` = signif(moran_p, 2))],
    paste("Janzen's prediction is a positive interaction. Indices on the 0-12 scale sum 12 monthly overlaps, so their coefficients",
          "are not comparable in size with those of the annual index. The all-cols model includes non-independent neighbouring",
          "cols and has residual spatial autocorrelation; it is shown only to compare the sign."))
fs2 <- ggplot(jz[!model %in% c("monthly", "envelope")], aes(b_int, lab)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = int_lo, xmax = int_hi), height = 0.2, linewidth = 0.3) +
  geom_point(shape = 21, size = 1.8, fill = oi[["blue"]], stroke = 0.2) +
  scale_y_discrete(limits = rev) + labs(x = "Height x latitude interaction (overlap per km per degree), 95% CI", y = NULL) + th
save2(fs2, file.path(sd, "FigS2_janzen_sensitivity"), 120, 55)

mc <- fread(file.path(tb, "M3_corredor_montano.csv"), colClasses = list(character = "ssp"))[zone == "tropical"]
me <- fread(file.path(tb, "M3_corredor_montano_ensamble.csv"), colClasses = list(character = "ssp"))[zone == "tropical"]
add("S2.3", "Table S2.3. Sensitivity of the montane gaps of the tropical divide to the threshold of montane climate.",
    rbind(mc[gcm == "present", .(`Threshold (m)` = threshold_m, Scenario = "Present", Period = "1970-2000",
                                 `Length in gaps (km), mean` = r(gap_km, 0), min = NA_real_, max = NA_real_,
                                 `Gaps longer than 10 km, mean` = n_gap10, `Longest gap (km), mean` = r(max_gap_km, 0))],
          me[, .(`Threshold (m)` = threshold_m, Scenario = ssp_lab[ssp], Period = period, `Length in gaps (km), mean` = r(gap_km.mean, 0),
                 min = r(gap_km.min, 0), max = r(gap_km.max, 0), `Gaps longer than 10 km, mean` = r(n_gap10.mean, 1),
                 `Longest gap (km), mean` = r(max_gap_km.mean, 0))])[order(`Threshold (m)`, Period != "1970-2000", Scenario, Period)],
    { se <- fread(file.path(tb, "M3_corredor_umbral.csv"))[order(threshold_m)]
      paste0(sprintf("Future values are the mean, minimum and maximum of the eight GCMs. Default threshold: %d m. ", P$montane_gap_m[[1]]),
             "Robust corridor loss under SSP3-7.0 in 2081-2100, critical gateways and share of the loss inside protected areas: ",
             paste(se[, sprintf("%d m, %.1f km, %d gateways, %.0f%%", threshold_m, loss_km, n_critical, pct_lost_in_PA)], collapse = "; "), ".") })

# =====================================================================================================================
# S3 - sensitivity of the H3 typology, and side that sets the moisture gap
# =====================================================================================================================
sides <- m2[flank_m == f1000 & col_id %in% gw$col_id, .(AI_W = AI_gap[side_lab == "W"], AI_E = AI_gap[side_lab == "E"],
                                                        CLD_W = CLD_gap[side_lab == "W"], CLD_E = CLD_gap[side_lab == "E"]), by = col_id]
s3 <- merge(merge(gw, typ, by = "col_id", suffixes = c("", ".t")), sides, by = "col_id")[order(-lat)]
s3[, AI_side := ifelse(AI_gap_side == "W", "Pacific (W)", "Atlantic (E)")]
s3[, AI_dir := ifelse(AI_gap_signed > 0, "Crest wetter", "Crest drier")]
add("S3.1", "Table S3.1. Typology of the 30 tropical gateways under the main rule and its three alternatives, and the side and direction of the aridity gap.",
    s3[, .(`APD id` = col_id, Gateway = name, Country = country, Latitude = r(lat, 3), `Crossing height (m)` = r(elev_cross, 0),
           `Thermal overlap` = r(TO_ann, 3), `Aridity gap, absolute` = r(AI_gap, 3), `Cloud gap (%), absolute` = r(CLD_gap, 1),
           `Aridity gap, W side` = r(AI_W, 3), `Aridity gap, E side` = r(AI_E, 3), `Side of larger absolute aridity gap` = AI_side, `Crest on that side` = AI_dir,
           `Type: medians (main)` = type_lab[h3_type], `Type: terciles` = type_lab[h3_type_terc],
           `Type: cloud gap as moisture axis` = type_lab[h3_type_cld], `Type: medians of all 61 gateways` = type_lab[h3_type_all])],
    paste("Main rule: thermal overlap at or above the median of the 30 tropical gateways is open in",
          "temperature; absolute aridity gap at or above the median is closed in moisture. Terciles: the middle third of either axis is",
          "intermediate. Each metric is that of the worse side (lower overlap, larger absolute gap); per-side gaps keep their sign, positive when the crest is wetter; the eastern side includes endorheic",
          "basins connected to the east."))
lv <- names(type_lab)
cnt <- function(x) as.integer(table(factor(x, levels = lv)))
add("S3.2", "Table S3.2. Number of tropical gateways of each type under the main rule and its alternatives.",
    data.table(Type = type_lab, `Medians (main)` = cnt(s3$h3_type), Terciles = cnt(s3$h3_type_terc),
               `Cloud gap as moisture axis` = cnt(s3$h3_type_cld), `Medians of all 61 gateways` = cnt(s3$h3_type_all),
               `Same type as main rule, terciles` = sapply(lv, \(k) s3[h3_type == k, sum(h3_type_terc == k)]),
               `Same type as main rule, cloud gap` = sapply(lv, \(k) s3[h3_type == k, sum(h3_type_cld == k)])),
    sprintf("The larger absolute aridity gap lies on the Pacific side at %d of the 30 gateways, with the crest wetter than that flank, and on the eastern side at %d, with the crest drier than that flank.", s3[AI_gap_signed > 0, sum(AI_side == "Pacific (W)")], s3[AI_gap_signed < 0, sum(AI_side == "Atlantic (E)")]))

# =====================================================================================================================
# S5 - GCMs: selection and results by model
# =====================================================================================================================
hz <- as.data.table(openxlsx::read.xlsx(file.path(raw, "refs", "Hausfather2022_supporting.xlsx"), sheet = "ECSTCR", startRow = 2))
setnames(hz, 1, "gcm")
inst <- c("ACCESS-ESM1-5" = "CSIRO (Australia)", "BCC-CSM2-MR" = "Beijing Climate Center (China)",
          "CMCC-ESM2" = "CMCC (Italy)", "CNRM-ESM2-1" = "CNRM-CERFACS (France)", "GISS-E2-1-G" = "NASA GISS (USA)",
          "MIROC6" = "MIROC consortium (Japan)", "MPI-ESM1-2-HR" = "Max Planck Institute for Meteorology (Germany)",
          "MRI-ESM2-0" = "Meteorological Research Institute (Japan)")
gcms <- sort(unique(fut$gcm)); stopifnot(setequal(gcms, names(inst)))
add("S5.1", "Table S5.1. The eight CMIP6 models used for future climates.",
    hz[gcm %in% gcms, .(Model = gcm, Institution = inst[gcm], `TCR (°C)` = r(TCR, 2), `ECS (°C)` = r(ECS150, 2))][order(Model)],
    paste("TCR and ECS from Hausfather et al. (2022), supporting data (Zelinka and IPCC AR6 Table 7.SM5; ECS from 150-year",
          "regressions). Filters, set before extracting any data: TCR within the AR6 likely range (1.4-2.2 °C) and ECS within",
          "the very likely range (2-5 °C); SSP1-2.6, SSP3-7.0 and SSP5-8.5 available in WorldClim 2.1; one model per institution",
          "or model family. Excluded by the last two filters: GFDL-ESM4 (no SSP5-8.5), ACCESS-CM2, AWI-CM-1-1-MR, GISS-E2-1-H,",
          "MIROC-ES2L and MPI-ESM1-2-LR."))
bym <- fut[, .(`Median change in thermal overlap` = r(median(dTO), 3), `Gateways with overlap increase` = sum(dTO > 0),
               `Median fraction of ΔT* reached` = r(median(frac_dTstar), 3), `Median change in precipitation gap` = r(median(dPg), 3),
               `Gateways with wider precipitation gap` = sum(dPg > 0), `Gateways with no lowland analogue` = sum(over_T)),
           by = .(Model = gcm, s = ssp, Period = period)][order(s, Period, Model)]
bym[, Scenario := ssp_lab[s]][, s := NULL]; setcolorder(bym, c("Model", "Scenario", "Period"))
add("S5.2", "Table S5.2. Future permeability of the 30 tropical gateways by GCM, scenario and period.", bym,
    paste("Future crest compared with the present flank, worse side; the precipitation gap is the larger absolute value of the two sides. ΔT* is the uniform warming that would centre the annual",
          "range of the crest on that of the present flank. No analogue: crest bio5 or bio6 beyond those of every present",
          "lowland cell below 1000 m within 100 km on the same side."))
add("S5.3", sprintf("Table S5.3. Montane gaps of the tropical divide (threshold %d m) by GCM, scenario and period.", P$montane_gap_m[[1]]),
    mc[threshold_m == P$montane_gap_m[[1]] & gcm != "present",
       .(Model = gcm, Scenario = ssp_lab[ssp], Period = period, `Length in gaps (km)` = r(gap_km, 0),
         `Gaps longer than 10 km` = n_gap10, `Longest gap (km)` = r(max_gap_km, 0),
         `Montane fragments longer than 10 km` = n_frag10, `Longest montane fragment (km)` = r(max_frag_km, 0))][order(Scenario, Period, Model)],
    sprintf("Present (1970-2000): %.0f km in gaps, %d gaps longer than 10 km.", mc[gcm == "present" & threshold_m == P$montane_gap_m[[1]], gap_km],
            mc[gcm == "present" & threshold_m == P$montane_gap_m[[1]], n_gap10]))
fd5 <- fut[, .(v = median(frac_dTstar)), by = .(gcm, ssp, period)]
fs5a <- ggplot(fd5, aes(period, v, group = gcm, colour = gcm)) + geom_line(linewidth = 0.35) + geom_point(size = 0.8) +
  facet_wrap(~ ssp, labeller = as_labeller(ssp_lab)) +
  scale_colour_manual(values = unname(c(oi[1:6], "black", "grey55"))) +
  labs(x = NULL, y = "Median fraction of ΔT* reached\n(30 tropical gateways)", tag = "a") + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + guides(colour = guide_legend(nrow = 2))
fs5b <- ggplot(mc[threshold_m == P$montane_gap_m[[1]] & gcm != "present"], aes(period, gap_km, group = gcm, colour = gcm)) +
  geom_hline(yintercept = mc[gcm == "present" & threshold_m == P$montane_gap_m[[1]], gap_km], colour = "grey60", linetype = "22", linewidth = 0.3) +
  geom_line(linewidth = 0.35) + geom_point(size = 0.8) + facet_wrap(~ ssp, labeller = as_labeller(ssp_lab)) +
  scale_colour_manual(values = unname(c(oi[1:6], "black", "grey55"))) +
  labs(x = NULL, y = "Tropical divide in montane gaps (km)", tag = "b") + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
save2(fs5a / fs5b, file.path(sd, "FigS5_gcm"), 170, 140)

# =====================================================================================================================
# S4 - gateway atlas: map, crossing-height profile and annual thermal ranges, one page per tropical gateway
# =====================================================================================================================
dem <- rast(file.path(raw, "hyd_sa_dem_15s", "hyd_sa_dem_15s.tif"))
prof <- fread(file.path(int, "M1_profile.csv"))
fl <- fread(file.path(int, "M2_flancos.csv"))[flank_m == f1000]
routes <- sf::st_read(file.path(int, "M2_flank_routes.gpkg"), quiet = TRUE)
fps <- fread(file.path(int, "M3_futuro_por_lado.csv"), colClasses = list(character = "ssp"))[period == "2081-2100"]
gaps <- fread(file.path(tb, "M3_dTstar.csv"))
# Level-1 units (department, province, state or region) within 500 m of each col: the divide is often their boundary
iso3 <- c(Venezuela = "VEN", Colombia = "COL", Ecuador = "ECU", Peru = "PER", Chile = "CHL", Bolivia = "BOL")
adm1 <- do.call(rbind, lapply(iso3, \(cc) sf::st_as_sf(unwrap(readRDS(file.path(raw, "gadm", "gadm", sprintf("gadm41_%s_1_pk.rds", cc)))))[, c("COUNTRY", "NAME_1")]))
adm_of <- function(lon, lat, country) {
  aeqd <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f", lat, lon)                  # metres around the col, no lwgeom needed
  a <- sf::st_transform(adm1[adm1$COUNTRY == country, ], aeqd); p <- sf::st_transform(sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326), aeqd)
  nm <- a$NAME_1[as.numeric(sf::st_distance(a, p)) <= 500]
  if (!length(nm)) nm <- a$NAME_1[sf::st_nearest_feature(p, a)]
  paste(sort(unique(nm)), collapse = " / ")
}
dms <- function(x, h) { a <- abs(x); d <- floor(a); m <- floor((a - d) * 60); sprintf("%d°%02d'%02.0f\"%s", d, m, ((a - d) * 60 - m) * 60, h[(x < 0) + 1]) }
atl <- list()
for (i in seq_len(nrow(s3))) {
  g <- s3[i]; cl <- cols[col_id == g$col_id]; id <- g$col_id
  fi <- fl[col_id == id]; ri <- routes[routes$col_id == id, ]
  # map: col, both flank points and their routes, at least +-0.4 degrees
  xs <- c(cl$lon, fi$lon); ys <- c(cl$lat, fi$lat); h <- max(0.4, diff(range(xs)) / 2 + 0.1, diff(range(ys)) / 2 + 0.1)
  e <- ext(mean(range(xs)) - h, mean(range(xs)) + h, mean(range(ys)) - h, mean(range(ys)) + h)
  z <- crop(dem, e); hs <- shade(terrain(z, "slope", unit = "radians"), terrain(z, "aspect", unit = "radians"), 40, 315)
  zd <- as.data.frame(z, xy = TRUE); names(zd)[3] <- "z"; hd <- as.data.frame(hs, xy = TRUE); names(hd)[3] <- "h"
  pd <- prof[lon >= e$xmin & lon <= e$xmax & lat >= e$ymin & lat <= e$ymax]
  ma <- ggplot() + geom_raster(data = zd, aes(x, y, fill = z)) +
    geom_raster(data = hd, aes(x, y, alpha = h), fill = "black", show.legend = FALSE) + scale_alpha(range = c(0.45, 0)) +
    scale_elev("fill", "Elevation (m)") +
    geom_path(data = pd, aes(lon, lat, group = grp), colour = "white", linewidth = 0.45) +
    geom_sf(data = ri, colour = oi[["verm"]], linewidth = 0.4, inherit.aes = FALSE) +
    geom_point(data = fi, aes(lon, lat, shape = side_lab), fill = oi[["sky"]], size = 1.6, stroke = 0.3) +
    geom_point(data = cl, aes(lon, lat), shape = 24, fill = oi[["yell"]], size = 2.2, stroke = 0.3) +
    scale_shape_manual(NULL, values = c(W = 22, E = 21), labels = c(W = "W flank", E = "E flank")) +
    coord_sf(xlim = c(e$xmin, e$xmax), ylim = c(e$ymin, e$ymax), expand = FALSE) + scale_x_continuous(n.breaks = 4) +
    labs(x = NULL, y = NULL, tag = "a") + th +
    theme(legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 6),
          legend.key.width = unit(6, "mm"), axis.line = element_blank(), panel.border = element_rect(fill = NA, linewidth = 0.3),
          plot.margin = margin(5.5, 14, 5.5, 5.5))                                # room for the last longitude label
  # profile +-150 km along the divide
  pp <- prof[abs(d_km - cl$d_km) <= 150]
  pb <- ggplot(pp, aes(d_km - cl$d_km, elev_cross, group = grp)) +
    geom_hline(yintercept = P$montane_gap_m[[1]], colour = oi[["sky"]], linetype = "22", linewidth = 0.3) +
    geom_line(linewidth = 0.3) + geom_point(data = cl, aes(0, elev_cross), inherit.aes = FALSE, shape = 24, fill = oi[["yell"]], size = 1.8) +
    labs(x = "Distance along the divide (km; north < 0 < south)", y = "Crossing height (m)", tag = "b") + th
  # annual thermal ranges (bio6 to bio5): flanks and crest now, crest 2081-2100 (mean of eight GCMs)
  mi <- m2[col_id == id & flank_m == f1000]
  fu <- fps[col_id == id, .(lo = mean(crest_b6), hi = mean(crest_b5)), by = ssp]
  rg <- rbind(mi[, .(k = paste(side_lab, "flank"), lo = flank_bio6, hi = flank_bio5)],
              mi[1, .(k = "Crest, present", lo = crest_bio6, hi = crest_bio5)],
              fu[, .(k = paste("Crest 2081-2100,", ssp_lab[ssp]), lo, hi)])
  rg[, k := factor(k, levels = rev(c("W flank", "E flank", "Crest, present", paste("Crest 2081-2100,", ssp_lab))))]
  rg[, grpc := ifelse(grepl("flank", k), "Flank", ifelse(grepl("present", k), "Present crest", "Future crest"))]
  pc <- ggplot(rg, aes(y = k, yend = k, x = lo, xend = hi, colour = grpc)) + geom_segment(linewidth = 2.2) +
    scale_colour_manual(values = c(Flank = oi[["verm"]], `Present crest` = oi[["blue"]], `Future crest` = oi[["sky"]]), guide = "none") +
    labs(x = "Annual temperature range, bio6 to bio5 (°C)", y = NULL, tag = "c") + th
  adm <- adm_of(cl$lon, cl$lat, g$country)
  pg <- ma + (pb / pc) + plot_layout(widths = c(1, 1.15)) +
    plot_annotation(title = g$name, subtitle = sprintf("%s, %s  ·  %s  %s (%.4f, %.4f)", adm, g$country, dms(cl$lat, c("N", "S")),
                                                       dms(cl$lon, c("E", "W")), cl$lat, cl$lon),
                    theme = theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 7.5, colour = "grey25")))
  fn <- sprintf("S4_%02d_%s", i, id); save2(pg, file.path(sd, "atlas", fn), 170, 120, dpi = 300)
  en370 <- ens[col_id == id & ssp == "370" & period == "2081-2100"]
  dt <- gaps[col_id == id]
  atl[[i]] <- data.table(No = i, File = paste0(fn, ".png"), `APD id` = id, Gateway = g$name, Country = g$country, `Level-1 unit` = adm,
                         Latitude = r(cl$lat, 3), Longitude = r(cl$lon, 3), `Crossing height (m)` = r(cl$elev_cross, 0),
                         `Depth (m)` = r(cl$depth, 0), `Isolation N (km)` = r(cl$iso_N_km, 0), `Isolation S (km)` = r(cl$iso_S_km, 0),
                         `Width (km)` = r(cl$width_km, 1),
                         `W flank route (km)` = r(mi[side_lab == "W", route_km], 1), `E flank route (km)` = r(mi[side_lab == "E", route_km], 1),
                         `Flank by nearest cell` = paste(mi[fallback == TRUE, side_lab], collapse = " "),
                         `Thermal overlap` = r(g$TO_ann, 3), `Aridity gap` = r(g$AI_gap_signed, 3), `Cloud gap (%)` = r(g$CLD_gap_signed, 1),
                         Type = type_lab[g$h3_type], `ΔT* (°C)` = r(max(dt$dTstar), 1),
                         `Fraction of ΔT*, SSP3-7.0 2081-2100` = r(en370$frac_dTstar.mean, 2),
                         `Change in overlap, SSP3-7.0 2081-2100` = sprintf("%.2f (%.2f to %.2f)", en370$dTO.mean, en370$dTO.min, en370$dTO.max),
                         `Models agreeing on sign` = en370$dTO.agree)
}
add("S4.1", "Table S4.1. Attributes of the 30 tropical gateways shown in the atlas.", rbindlist(atl),
    paste("Order from north to south. Flank route is the length of the steepest-descent path to 1000 m; 'Flank by nearest cell'",
          "lists the sides where the path did not reach 1000 m within 150 km. ΔT* is that of the worse side. Change in overlap:",
          "mean of the eight GCMs (range)."))

# =====================================================================================================================
# Tables behind the figures quoted in the Results
# =====================================================================================================================
nm <- cols[, .(col_id, name)]; nm[gw, name := fifelse(is.na(name) | name == "", i.name, name), on = "col_id"]
zone <- function(lat) fifelse(lat >= -23.5, "Tropical", "Temperate")

# S1.4 present montane gaps of the tropical divide, located by their lowest point
mg <- fread(file.path(int, "M1_montane_gaps.csv"))[threshold_m == P$montane_gap_m[[1]] & lat_start >= -23.5][order(-lat_start)]
ce <- names(prof)[names(prof) %in% c("elev_cross", "elev_min", "cross", "e")][1]
s14 <- rbindlist(lapply(seq_len(nrow(mg)), function(i) {
  x <- mg[i]; k <- prof[lat <= x$lat_start & lat >= x$lat_end]; k <- k[which.min(k[[ce]])]
  co <- cols[lat <= x$lat_start + 1e-6 & lat >= x$lat_end - 1e-6][order(elev_cross)][1]
  pt <- sf::st_sfc(sf::st_point(c(k$lon, k$lat)), crs = 4326)
  ct <- adm1$COUNTRY[sf::st_nearest_feature(pt, adm1)]
  data.table(Gap = i, `North end (latitude)` = r(x$lat_start, 3), `South end (latitude)` = r(x$lat_end, 3), `Longitude of lowest point` = r(k$lon, 3),
             Country = ct, `Level-1 unit` = adm_of(k$lon, k$lat, ct), `Length (km)` = r(x$length_km, 1), `Lowest crossing height (m)` = r(x$min_elev, 0),
             `At least 10 km` = x$length_km >= 10, `Lowest col, APD id` = co$col_id, `Lowest col, name` = nm[col_id == co$col_id, name])
}))
add("S1.4", sprintf("Table S1.4. Present montane gaps of the tropical divide (crossing height below %d m), from north to south.", P$montane_gap_m[[1]]), s14,
    sprintf("Gaps are continuous runs of the crossing-height profile below the threshold. The %d gaps add up to %.0f km, and %d of them are at least 10 km long. Level-1 unit: department, state or region within 500 m of the lowest point; the divide is often their boundary.",
            nrow(s14), sum(s14$`Length (km)`), sum(s14$`At least 10 km`)))

# S1.5 most isolated cols of the tropics
s15 <- merge(cols[tropical & !terminal & depth >= 100][order(-iso_km)][1:10], nm, by = "col_id", suffixes = c(".old", ""), sort = FALSE)
add("S1.5", "Table S1.5. The ten most isolated non-terminal tropical cols with depth of at least 100 m.",
    s15[, .(Rank = .I, `APD id` = col_id, Col = name, Latitude = r(lat, 3), Longitude = r(lon, 3), `Crossing height (m)` = r(elev_cross, 0), `Depth (m)` = r(depth, 0),
            `Isolation N (km)` = r(iso_N_km, 0), `Isolation S (km)` = r(iso_S_km, 0), `Isolation, smaller value (km)` = r(iso_km, 0), Gateway = gateway)],
    "Isolation is the distance along the divide to the nearest lower point in each direction. Ranked by the smaller value. Every col is in APD_cols.csv.")

# S2.4 / S2.5 the 61 non-terminal gateways and their medians by zone and latitude band
jm <- fread(file.path(tb, "M2_janzen_mecanismo.csv"))
sd1 <- m2[flank_m == f1000 & col_id %in% jm$col_id]
wk <- merge(jm, sd1[, .(col_id, side_lab, crest_bio5, crest_bio6, flank_bio5, flank_bio6)], by = c("col_id", "side_lab"))
wk[, `:=`(d_C = pmin(crest_bio5, flank_bio5) - pmax(crest_bio6, flank_bio6), R_c = crest_bio5 - crest_bio6, R_f = flank_bio5 - flank_bio6, zona = zone(lat))]
stopifnot(nrow(wk) == 61, wk[, max(abs(TO_ann - d_C / sqrt(R_c * R_f)))] < 1e-6)
wk <- merge(wk, sd1[, .(backup = paste(side_lab[fallback == TRUE], collapse = " ")), by = col_id], by = "col_id")[order(-lat)]
wk <- merge(wk, nm, by = "col_id", sort = FALSE)
add("S2.4", "Table S2.4. Thermal overlap of the 61 non-terminal gateways, worse side, 1000 m flank, from north to south.",
    wk[, .(`APD id` = col_id, Gateway = name, Zone = zona, Latitude = r(lat, 3), `Crossing height (m)` = r(elev_cross, 0), `Worse side` = side_lab,
           `R crest (°C)` = r(R_c, 1), `R flank (°C)` = r(R_f, 1), `d (°C)` = r(d_C, 1), TOann = r(TO_ann, 3), TOmon = r(TO_mon, 3), TOenv = r(TO_env, 3),
           `Mean daily range of flank (°C)` = r(flank_diurnal, 1), `Backup flank` = backup)],
    "R is the annual range (bio5 - bio6) and d the shared length of the annual ranges of crest and flank, negative when they do not overlap, so that TOann = d / sqrt(R crest x R flank). Backup flank lists the sides where the descent path did not reach 1000 m within 150 km.")
bk <- sd1[, .(lat = lat[1], w = sum(side_lab == "W" & fallback), e = sum(side_lab == "E" & fallback)), by = col_id]
sum25 <- function(lbl, x, b) data.table(Group = lbl, n = nrow(x), `Median crossing height (m)` = r(median(x$elev_cross), 0), `Median TOann` = r(median(x$TO_ann), 2),
  `Median d (°C)` = r(median(x$d_C), 1), `Gateways with ranges that do not overlap` = sum(x$TO_ann < 0), `Median annual range of flank (°C)` = r(median(x$flank_annual), 1),
  `Median mean daily range of flank (°C)` = r(median(x$flank_diurnal), 1), `Backup flanks, western side` = sum(b$w), `Backup flanks, eastern side` = sum(b$e))
s25 <- rbind(sum25("Tropical", wk[lat >= -23.5], bk[lat >= -23.5]), sum25("Temperate", wk[lat < -23.5], bk[lat < -23.5]),
             sum25("All non-terminal", wk, bk),
             rbindlist(lapply(list(c(0, 12), c(12, 23.5), c(23.5, 36), c(36, 60)), function(q)
               sum25(sprintf("%g to %g° of absolute latitude", q[1], q[2]), wk[abs(lat) > q[1] & abs(lat) <= q[2]], bk[abs(lat) > q[1] & abs(lat) <= q[2]]))))
add("S2.5", "Table S2.5. Median crossing height, thermal overlap and flank temperature ranges of the 61 non-terminal gateways by zone and latitude band.", s25,
    "Medians of Table S2.4. A gateway has two sides, so 61 gateways have 122 sides.")

# S3.3 correlations and GAM deviance of the moisture axis
h3 <- fread(file.path(tb, "M2_H3_resumen.csv"))
vl <- c(TO_ann = "TOann", AI_gap = "aridity gap", CLD_gap = "cloud gap", elev_cross = "crossing height", elev_km = "crossing height")
add("S3.3", "Table S3.3. Moisture axis: Spearman correlations among the 30 tropical gateways and across the 61 gateways, and deviance of the GAMs.",
    rbind(h3[test != "gam", .(Analysis = fifelse(test == "spearman", "Spearman, tropical", "Spearman, all gateways"), Variables = paste(vl[x], "and", vl[y]), n,
                              rho = r(est, 2), p = signif(p, 2), `Deviance, smooth alone` = NA_real_, `Deviance, full model` = NA_real_)],
          unique(h3[test == "gam", .(y, dev_smooth, dev_expl)])[, .(Analysis = "GAM", Variables = paste(vl[y], "on crossing height, TOann and position along the divide"), n = 30L,
                              rho = NA_real_, p = NA_real_, `Deviance, smooth alone` = r(dev_smooth, 2), `Deviance, full model` = r(dev_expl, 2))]),
    "Moisture gaps are the absolute values of the worse side. Terms of the GAMs are in Table 2a.")

# S5.4 future metrics by gateway, S5.5 by type, S5.6 to S5.8 corridor
e4 <- merge(ens[col_id %in% gw$col_id], typ[, .(col_id, h3_type, lat, elev_cross)], by = "col_id")
e4 <- merge(e4, gaps[, .(dT = max(dTstar)), by = col_id], by = "col_id")
e4 <- merge(e4, nm, by = "col_id")[order(-lat, ssp, period)]
add("S5.4", "Table S5.4. Future metrics of the 30 tropical gateways by scenario and period, mean and range of the eight GCMs.",
    e4[, .(`APD id` = col_id, Gateway = name, Latitude = r(lat, 3), Type = type_lab[h3_type], Scenario = ssp_lab[ssp], Period = period, `TOann, present` = r(TO_now, 3),
           `Change in TOann, mean` = r(dTO.mean, 2), min = r(dTO.min, 2), max = r(dTO.max, 2), `Models sharing the sign of the mean` = round(dTO.agree * 8),
           `ΔT* (°C)` = r(dT, 1), `Fraction of ΔT*, mean` = r(frac_dTstar.mean, 2), `Fraction, min` = r(frac_dTstar.min, 2), `Fraction, max` = r(frac_dTstar.max, 2),
           `Change in precipitation gap, mean` = r(dPg.mean, 2), `Models sharing that sign` = round(dPg.agree * 8), `Share of models with no lowland analogue` = r(p_over_T, 2))],
    "Future crest compared with the present flank, worse side; the precipitation gap is the larger absolute value of the two sides.")
p0 <- sd1[, .(col_id, side_lab, b5c = crest_bio5, b6c = crest_bio6, b5f = flank_bio5, b6f = flank_bio6, TO0 = TO_ann)]
fv <- merge(fps[ssp == "370"], p0, by = c("col_id", "side_lab"))
fv[, `:=`(dcrest = (crest_b5 + crest_b6 - b5c - b6c) / 2, dflank = (flank_b5 + flank_b6 - b5f - b6f) / 2,
          TOff = (pmin(crest_b5, flank_b5) - pmax(crest_b6, flank_b6)) / sqrt((crest_b5 - crest_b6) * (flank_b5 - flank_b6)))]
fw <- fv[col_id %in% gw$col_id, .(dc = mean(dcrest), df = mean(dflank), TOff = mean(TOff), TO0 = TO0[1]), by = .(col_id, side_lab)]
tyv <- e4[ssp == "370" & period == "2081-2100", .(n = .N, m = median(dTO.mean), a = median(dPg.mean), s = sum(dPg.agree >= 7 / 8)), by = h3_type]
add("S5.5", "Table S5.5. Change in thermal overlap and in the precipitation gap by gateway type, SSP3-7.0, 2081-2100 (median of the means of eight GCMs).",
    tyv[match(names(type_lab)[1:4], h3_type), .(Type = type_lab[h3_type], n, `Median change in TOann` = r(m, 2), `Median change in precipitation gap` = r(a, 2),
                                                `Gateways with robust sign of the precipitation change` = s)],
    sprintf("Across the 60 sides of the 30 tropical gateways, the warming of the midpoint of the annual range differed between crest and flank by a median absolute value of %.1f °C (mean of eight GCMs), and TOann between future crest and future flank stayed within %.2f of its present value.",
            fw[, median(abs(dc - df))], fw[, max(abs(TOff - TO0))]))
pc <- fread(file.path(tb, "M3_perdida_pais.csv"))
add("S5.6", "Table S5.6. Montane corridor of the tropical divide robustly lost by country, SSP3-7.0, 2081-2100.",
    rbind(pc[order(-lost_km), .(Country = country, `Corridor today (km)` = r(corridor_km, 1), `Robust loss (km)` = r(lost_km, 1), `Share of corridor lost (%)` = r(pct_lost, 1))],
          pc[, .(Country = "Total", `Corridor today (km)` = r(sum(corridor_km), 1), `Robust loss (km)` = r(sum(lost_km), 1), `Share of corridor lost (%)` = r(100 * sum(lost_km) / sum(corridor_km), 1))]),
    "Robust loss: montane today and below the raised isotherm in at least seven of eight GCMs. Divide segments on a border go to one country.")
cp <- fread(file.path(tb, "M3_corredor_proteccion.csv"))
add("S5.7", "Table S5.7. Protection and human modification of the montane corridor of the tropical divide, SSP3-7.0, 2081-2100.",
    cp[, .(Stretch = stretch, `Length (km)` = r(km, 1), `Inside protected areas (%)` = r(pct_in_PA, 1), `Mean gHM` = r(ghm_mean, 2), `With gHM of 0.4 or more (%)` = r(pct_high_pressure, 1))],
    "Each 0.5 km segment was evaluated at its lowest crossing. Shares with the other thresholds are in the note of Table S2.3.")
bp <- fread(file.path(int, "M3_corredor_bins_proteccion.csv"))[tropical == TRUE]
add("S5.8", "Table S5.8. Local lapse rate and rise of the isotherm that now lies at 2000 m along the tropical divide, SSP3-7.0, 2081-2100, by 1° of latitude.",
    bp[, .(`Lapse rate (°C per km)` = r(median(abs(G_C_per_km)), 2), `Divide segments` = .N, `Median rise (m)` = r(median(rise_m_370_2090), 0),
           `Minimum rise (m)` = r(min(rise_m_370_2090), 0), `Maximum rise (m)` = r(max(rise_m_370_2090), 0)), by = .(`Latitude window` = floor(lat))][order(-`Latitude window`)],
    sprintf("Rise: mean warming of the segment divided by the local lapse rate, mean of eight GCMs. Whole tropical divide: median %.0f m, range %.0f to %.0f m; lapse rates %.1f to %.1f °C per km.",
            median(bp$rise_m_370_2090), min(bp$rise_m_370_2090), max(bp$rise_m_370_2090), min(abs(bp$G_C_per_km)), max(abs(bp$G_C_per_km))))

# ---- one CSV per table plus contents.csv (titles and notes); 34_build_supp.py turns them into the workbook and the Word file.
# ponytail: openxlsx 4.2.9 writes dangling drawing links that openpyxl (and Excel's repair prompt) choke on, so Python writes the xlsx.
td <- file.path(sd, "tables"); dir.create(td, showWarnings = FALSE)
for (k in names(sheets)) fwrite(sheets[[k]], file.path(td, paste0(k, ".csv")), na = "", bom = TRUE)
ct <- rbindlist(meta)[order(sheet)]; fwrite(ct, file.path(td, "contents.csv"), bom = TRUE)
cat("S1-S5 tables and figures written to", sd, "\n"); print(ct[, .(sheet, n = sapply(sheets[sheet], nrow))])
