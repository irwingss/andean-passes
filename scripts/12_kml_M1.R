# 12_kml_M1.R - KML of the M1 census for Google Earth: tropical gateways (orange), Lara cols (yellow),
# other tropical cols (grey, hidden by default), extratropical gateways (blue, hidden) and the W|E divide.
# Same selection as M1_elev_vs_aislamiento (depth >= 100 m, chain ends excluded).
# Rscript scripts/12_kml_M1.R
cfg <- yaml::read_yaml("config.yml"); int <- cfg$paths$interim; P <- cfg$params
cols <- read.csv("outputs/database/APD_cols.csv")
prof <- read.csv(file.path(int, "M1_profile.csv"))
dir.create("outputs/kml", showWarnings = FALSE)
esc <- function(x) { x <- gsub("&", "&amp;", x); x <- gsub("<", "&lt;", x); gsub(">", "&gt;", x) }
ll <- function(lat, lon) sprintf("%.4f°%s, %.4f°W", abs(lat), ifelse(lat >= 0, "N", "S"), abs(lon))

c0 <- cols[cols$depth >= P$depth_m[[1]] & !cols$terminal, ]
groups <- list(
  list(id = "gw_trop", name = "Puertas tropicales", sel = c0$tropical & c0$gateway & !c0$lara, col = "ff005ed5", icon = "shaded_dot", vis = 1),
  list(id = "lara", name = "Lara (extremo NE de los Andes de Mérida)", sel = c0$tropical & c0$lara, col = "ff009fe6", icon = "triangle", vis = 1),
  list(id = "col_trop", name = "Otros collados tropicales (depth >= 100 m)", sel = c0$tropical & !c0$gateway & !c0$lara, col = "ff999999", icon = "shaded_dot", vis = 0),
  list(id = "gw_extra", name = "Puertas extratropicales", sel = !c0$tropical & c0$gateway, col = "ffb27200", icon = "triangle", vis = 0))

style <- function(g) sprintf(paste0('<Style id="%s"><IconStyle><color>%s</color><scale>%s</scale><Icon><href>',
  'http://maps.google.com/mapfiles/kml/shapes/%s.png</href></Icon></IconStyle><LabelStyle><scale>%s</scale></LabelStyle></Style>'),
  g$id, g$col, if (g$id == "col_trop") "0.5" else "0.9", g$icon, if (g$vis) "0.8" else "0")
pm <- function(r, g) {
  nm <- if (!is.na(r$name)) r$name else sprintf("%s %.0f m", ll(r$lat, r$lon), r$elev_cross)
  desc <- sprintf(paste0("Altura de cruce W–E: %.0f m<br/>Elevación topográfica: %.0f m<br/>Profundidad (depth): %.0f m<br/>",
                         "Aislamiento: %.1f km (N %.1f / S %.1f)<br/>Ancho (&lt; col + 100 m): %.1f km<br/>Puerta: %s<br/>",
                         "Bordea cuenca endorreica: %s<br/>%s<br/>ID: %s"),
                  r$elev_cross, r$elev, r$depth, r$iso_km, r$iso_N_km, r$iso_S_km, r$width_km,
                  ifelse(r$gateway, "sí", "no"), ifelse(r$borders_I, "sí", "no"), ll(r$lat, r$lon), r$col_id)
  sprintf('<Placemark><name>%s</name><description><![CDATA[%s]]></description><styleUrl>#%s</styleUrl><Point><coordinates>%.6f,%.6f,0</coordinates></Point></Placemark>',
          esc(nm), desc, g$id, r$lon, r$lat)
}
folders <- vapply(groups, function(g) {
  d <- c0[g$sel, ]; d <- d[order(d$elev_cross), ]
  sprintf("<Folder><name>%s (%d)</name><visibility>%d</visibility>%s</Folder>", esc(g$name), nrow(d), g$vis,
          paste(vapply(seq_len(nrow(d)), function(i) pm(d[i, ], g), ""), collapse = "\n"))
}, "")
# divide: one LineString per continuous stretch, every 4th 0.5 km bin (~2 km) to keep the file light
lines <- vapply(split(prof, prof$grp), function(z) {
  z <- z[unique(c(seq(1, nrow(z), 4), nrow(z))), ]
  sprintf("<Placemark><styleUrl>#div</styleUrl><LineString><tessellate>1</tessellate><coordinates>%s</coordinates></LineString></Placemark>",
          paste(sprintf("%.5f,%.5f,0", z$lon, z$lat), collapse = " "))
}, "")
kml <- c('<?xml version="1.0" encoding="UTF-8"?>',
         '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
         '<name>APD v1.0 · Collados y puertas de los Andes (GLO-30)</name>',
         '<description>Divisoria W|E de los Andes derivada de Copernicus GLO-30 (inundación desde marcadores). Collados con depth ≥ 100 m. Puerta = mínimo de la divisoria en ±100 km.</description>',
         vapply(groups, style, ""),
         '<Style id="div"><LineStyle><color>ff000000</color><width>2</width></LineStyle></Style>',
         sprintf("<Folder><name>Divisoria W|E</name>%s</Folder>", paste(lines, collapse = "\n")),
         folders, "</Document></kml>")
out <- "outputs/kml/M1_collados_andes.kml"
writeLines(enc2utf8(kml), out, useBytes = TRUE)
cat(out, ":", sapply(groups, function(g) sum(g$sel)), "\n")
