# 27_tablas_manuscrito.R - Article 1 main tables (plan v2.0 section 8.2), from M1-M3 outputs only.
# Table 1: the 15 lowest tropical gateways (census attributes, present climate at flank 1000 m worst side, moisture gaps with
#          their sign, H3 type,
#          change in 2081-2100 under SSP3-7.0: GCM mean, range and sign agreement).
# Table 2: H2 and H3 model terms (estimate, 95 % CI, p, deviance explained, n, Moran's I) + H4 summary by SSP and period.
# Table 3: conservation priority of the 30 tropical gateways (protection, gHM, E-W opening, montane bottleneck).
# Rscript scripts/27_tablas_manuscrito.R
suppressPackageStartupMessages({ library(data.table); library(mgcv) })
cfg <- yaml::read_yaml("config.yml"); tb <- file.path(cfg$paths$outputs, "tables"); P <- cfg$params

cols <- fread("outputs/database/APD_cols.csv", encoding = "UTF-8")
h3 <- fread(file.path(tb, "M2_H3_tipologia.csv"))
m2 <- fread(file.path(tb, "M2_metricas_clima_collados.csv"))[flank_m == P$flank_m[1]]
ens <- fread(file.path(tb, "M3_futuro_ensamble.csv"))

# ---- Table 1 ----
t1 <- cols[col_id %in% h3$col_id][order(elev_cross)][1:15]
t1[, country := maps::map.where("world", lon, lat)]
t1[is.na(country), country := maps::map.where("world", lon + 0.05, lat)]   # divide points on a border line
t1[is.na(name) & abs(lat - 7.686) < 0.01, name := "Táchira depression"]   # regional name, as in Fig. 1
rk <- dcast(m2[col_id %in% t1$col_id], col_id ~ side_lab, value.var = "route_km")
f <- ens[ssp == 370 & period == "2081-2100"]
t1 <- Reduce(function(a, b) merge(a, b, by = "col_id", sort = FALSE), list(
  t1[, .(col_id, name, country, lat = round(lat, 3), lon = round(lon, 3), elev_cross = round(elev_cross),
         depth = round(depth), iso_N_km = round(iso_N_km), iso_S_km = round(iso_S_km), width_km = round(width_km, 1))],
  rk[, .(col_id, flank_W_km = round(W, 1), flank_E_km = round(E, 1))],
  h3[, .(col_id, TO_ann = round(TO_ann, 2), AI_gap = round(AI_gap_signed, 2), CLD_gap = round(CLD_gap_signed, 1), h3_type)],
  f[, .(col_id, dTO_2090_ssp370 = sprintf("%+.2f [%+.2f, %+.2f]", dTO.mean, dTO.min, dTO.max), dTO_agree = dTO.agree,
        dPg_2090_ssp370 = sprintf("%+.2f [%+.2f, %+.2f]", dPg.mean, dPg.min, dPg.max), dPg_agree = dPg.agree)]))
fwrite(t1[order(elev_cross)], file.path(tb, "Table1_gateways.csv"))

# ---- Table 2 ----
jm <- readRDS(file.path(tb, "M2_modelo_janzen.rds"))$main
r <- c("elev_km", "abs_lat", "elev_km:abs_lat")
h2 <- data.table(hyp = "H2", response = "TO_ann", term = r, n = jm$n, est = jm$coef[r, 1], lo = jm$ci[r, "lo"], hi = jm$ci[r, "hi"],
                 p = jm$coef[r, 4], dev_expl = summary(jm$model)$dev.expl, moran_I = jm$moran[[1]], moran_p = jm$moran[[2]])
h3m <- fread(file.path(tb, "M2_H3_resumen.csv"))
h3m <- h3m[, .(hyp = "H3", response = y, term = ifelse(test == "gam", x, paste0(test, ": ", x)), n, est, lo, hi, p, dev_expl, moran_I, moran_p)]
fwrite(rbind(h2, h3m), file.path(tb, "Table2_models.csv"))
e <- ens[col_id %in% h3$col_id]
h4 <- e[, .(n = .N, dTO_median = median(dTO.mean), TO_up_robust = sum(dTO.mean > 0 & dTO.agree >= 0.8),
            TO_down_robust = sum(dTO.mean < 0 & dTO.agree >= 0.8), dPg_median = median(dPg.mean),
            Pg_up_robust = sum(dPg.mean > 0 & dPg.agree >= 0.8), Pg_down_robust = sum(dPg.mean < 0 & dPg.agree >= 0.8),
            no_analogue = sum(p_over_T > 0.5), frac_dTstar_median = median(frac_dTstar.mean, na.rm = TRUE)), by = .(ssp, period)]
fwrite(h4, file.path(tb, "Table2_H4_summary.csv"))
print(t1[order(elev_cross)]); print(h4)

# ---- Table 3: conservation priority of the tropical gateways (30_proteccion_M3.R) ----
pr <- fread(file.path(tb, "M3_prioridad_puertas.csv"), encoding = "UTF-8")
pr[is.na(name) & abs(lat - 7.686) < 0.01, name := "Táchira depression"]
pr <- merge(pr, t1[, .(col_id, country)], all.x = TRUE)[is.na(country), country := maps::map.where("world", lon, lat)]
pr[is.na(country), country := maps::map.where("world", lon + 0.05, lat)]
fwrite(pr[order(-critical, -lat), .(col_id, name, country, lat = round(lat, 3), elev_cross = round(elev_cross), h3_type,
                                   crest_in_PA, pct_PA_5km = round(pct_PA_5km), routes_pct_PA = round(routes_pct_PA),
                                   ghm_crest = round(ghm_crest, 2), frac_dTstar_370_2090 = round(frac_dTstar_370_2090, 2),
                                   montane_loss_km, opening, montane_bottleneck, critical, critical_unprotected, high_pressure)],
       file.path(tb, "Table3_conservation_priority.csv"))
