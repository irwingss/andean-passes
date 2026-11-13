"""53_suplemento.py - Supplementary Information of the article: outputs/supp/L1_Supplementary_Information.docx (captions,
short tables and figures) and outputs/supp/L1_Supplementary_Tables.xlsx (every table in full, plus the dictionary of
the Andean Passes Database). Reads outputs/supp/tables (51_sintesis_conectividad.py) and outputs/supp/FigS*.png
(52_figuras.py).
python -P scripts/53_suplemento.py
"""
import pandas as pd
from docx import Document
from docx.shared import Cm, Pt

T, F = "outputs/supp/tables", "outputs/supp"
TABLES = [
    ("S1", "S1_census_by_unit", "Passes of the census by zone, GMBA unit and rank: number of passes with a depth of 30 m or "
     "more, number with a depth of 100 m or more, and median height and depth (m). Tropics are 11°N to 23.5°S."),
    ("S2", "S2_control_passes_census", "The 36 control passes inside the area: the nearest pass (near_) and the nearest "
     "pass with a depth of 100 m or more (deep_), with distance (km), height (m), depth (m), rank, separated area (km²), width (m) and "
     "isolation (km). Measured from relief alone."),
    ("S3", "S3_osm_nodes_identity", "OpenStreetMap pass nodes (mountain_pass=yes or natural=saddle, downloaded 8 "
     "November 2026): share with a pass within 500 m and 1 km by unit and rank (nodes_x), and the identity test for "
     "nodes on the divide network (nodes_y): observed share against the mean and 97.5th percentile of 1000 draws of "
     "random divide cells of the same block core and unit, with the empirical p value."),
    ("S4", "S4_reference_heights", "Official heights of passes (INEI, Peru; UPF, Chile), the nearest census pass and its error "
     "(census pass minus reference, m)."),
    ("S5", "S5_height_error_by_zone", "Error of the height of the nearest pass by zone, against official heights and "
     "against the ele tag of OpenStreetMap (which may itself come from an elevation model)."),
    ("S6", "S6_osm_nodes_without_pass", "OpenStreetMap pass nodes without a census pass within 1 km."),
    ("S7", "S7_envelopes", "Montane limit T2000(lat) (°C), interquartile range of each variable over the domain, and the "
     "5th and 95th percentiles of each variable in the lowland and montane guilds, every 0.5° of latitude. "
     "Precipitation variables as log(1 + P)."),
    ("S8", "S8_tile_conductance_profiles", "Relative conductance of the tiles per configuration, direction, variable and "
     "0.5° band of latitude of the tile centre: number of tiles, median and quartiles."),
    ("S9", "S9_control_huancabamba_tachira", "Control: relative conductance of the tiles that cross the continental "
     "divide for longitudinal connectivity with BIO1, per configuration and 0.5° band, around the Huancabamba zone "
     "(2 to 10°S) and the Táchira depression (74 to 70°W, north of 5°N). in_zone marks the bands of each control "
     "zone. A band next to a 0.5° zone is not an independent comparison when its tiles still cover the zone, which "
     "is the case north of the Táchira depression with tiles of 200 and 400 km."),
    ("S10", "S10_spearman_variables", "Spearman correlation between pairs of variables in the normalized current of the "
     "passes (one pass per 1 km cell) and in the relative conductance of the tiles, main configuration."),
    ("S11", "S11_variable_agreement", "Agreement of each variable with the others: median and minimum Spearman "
     "correlation with the other 18 variables at the level of passes, median at the level of tiles, and correlation "
     "with the pass median over the 19 variables. Main configuration."),
    ("S12", "S12_sensitivity", "Sensitivity of the pass medians against the main configuration: cells with a value in "
     "both, cells that lose their value, Spearman correlation, share of the top decile that stays in the top decile, "
     "and share of cells on the same side of 1."),
    ("S13", "S13_by_unit", "Pass medians by GMBA unit and direction: cells, median and 90th percentile of the median over "
     "the 19 variables, and number and share of cells that concentrate (above 1) or divert (below 1) current in all "
     "five configurations."),
    ("S14", "S14_by_latitude", "Pass medians by band of 1° of latitude and direction: cells, median, 90th percentile and "
     "share of cells above 1 in all five configurations."),
    ("S15", "S15_controls_by_variable", "Mean of the 0.5° band medians of the relative conductance of the tiles that "
     "cross the continental divide inside each control zone of the Methods (Huancabamba 4 to 7.5°S; Táchira 7.5 to "
     "8°N, 74 to 70°W) and in the 2° of latitude next to it (Huancabamba north 2 to 4°S and south 7.5 to 9.5°S; "
     "Táchira south 5 to 7.5°N), per variable and direction, main configuration. North of the Táchira depression a "
     "tile of 200 km still covers it, so that side is left out here and given in Table S9. Lowest in zone is true "
     "when the zone is below every neighbour compared."),
    ("S16", "S16_by_height", "Pass medians by band of 500 m of pass height and direction (one pass per 1 km cell, the "
     "deepest; the last band holds every pass above 5500 m): cells, median, 10th and 90th percentiles of the median "
     "over the 19 variables, and share of cells that concentrate (above 1) or divert (below 1) current in all five "
     "configurations."),
    ("S17", "S17_divide_tile_profiles", "Relative conductance of the tiles that cross the continental divide, per 0.5° "
     "band of latitude of the tile centre, direction, group of variables and configuration. A tile crosses the divide "
     "when continental passes lie in its central half along the direction of the current, that is, between its two "
     "walls. Per tile the group is summarized by its median and by its most and least limiting variable, and the "
     "table gives the median of each over the tiles of the band. Table S8 is the same profile over every tile."),
    ("S18", "S18_limiting_variable", "Share of the tiles that cross the continental divide in which each variable is "
     "the most limiting of its group, by direction, main configuration."),
    ("S19", "S19_profile_bands", "Bands of local extremes of the temperature median in the profiles of Table S17 "
     "(main configuration), numbered as the grey bands of Fig. 3c and 3d and the black stretches of the continental "
     "divide in Fig. 3a and 3b, north to south: 1 to 7 are local maxima of transverse connectivity (openings across "
     "the cordillera) and 8 to 10 local minima of longitudinal connectivity (barriers along it). North and south are "
     "the edges of the band (°, edges of 0.5° bins), with its width and number of bins and the range of tiles per bin. "
     "For each group of variables, the mean over the bins of the band of the median and of the most and least limiting "
     "variable, the extreme bin median (highest for transverse, lowest for longitudinal) and the mean median over the "
     "1° on either side of the band (neighbours). Lowest pass is the lowest continental pass with a depth of 100 m or "
     "more inside the band, with its height and depth (m); control passes are those of Table 2 inside the band."),
]
INLINE_MAX = 60       # longer tables are given only in the workbook
HEAD = {"group": "Group", "most_limiting": "Most limiting", "least_limiting": "Least limiting",
        "share_of_tiles": "Share of tiles", "zone": "Zone", "gmba": "GMBA unit", "rank": "Rank", "passes": "Passes", "depth_100m": "Depth ≥ 100 m",
        "median_z_m": "Median height (m)", "median_depth_m": "Median depth (m)", "control": "Pass", "system": "System",
        "lat": "Lat (°)", "lon": "Lon (°)", "group": "Group", "reading": "Direction", "var": "Variable",
        "configuration": "Configuration", "cells": "Cells", "median": "Median", "p10": "P10", "p90": "P90",
        "z_band": "Height band (m)", "band": "Band (°)", "reference": "Reference", "n": "n",
        "median_error_m": "Median error (m)", "median_abs_error_m": "Median absolute error (m)",
        "share_within_50m": "Within 50 m", "cells_both": "Cells in both", "cells_lost": "Cells lost",
        "rho_with_main": "ρ with main", "top_decile_kept": "Top decile kept", "same_side_of_1": "Same side of 1",
        "n_robust_above_1": "n above 1", "n_robust_below_1": "n below 1", "robust_above_1": "Share above 1",
        "robust_below_1": "Share below 1", "median_rho_other_vars_passes": "Median ρ, passes",
        "min_rho_other_vars_passes": "Min. ρ, passes", "median_rho_other_vars_tiles": "Median ρ, tiles",
        "rho_with_pass_median": "ρ with pass median", "north": "North", "south": "South",
        "lowest_in_zone": "Lowest in zone", "name": "Pass", "country": "Country", "elev_ref_m": "Official (m)",
        "census_km": "Census pass distance (km)", "census_z": "Census pass height (m)", "census_depth": "Census pass depth (m)", "error_m": "Error (m)",
        "near_km": "Nearest pass (km)", "near_z": "Height (m)", "near_depth": "Depth (m)", "near_rank": "Rank",
        "near_area_sep_km2": "Separated area (km²)", "near_isolation_km": "Isolation (km)",
        "deep_km": "Deep pass (km)", "deep_z": "Deep height (m)", "deep_depth": "Deep depth (m)",
        "deep_rank": "Deep rank", "deep_area_sep_km2": "Deep separated area (km²)", "deep_width_m": "Deep width (m)",
        "deep_isolation_km": "Deep isolation (km)", "nodes_x": "Nodes", "within_500m": "Within 500 m",
        "within_1km": "Within 1 km", "nodes_y": "Nodes on network", "obs_500m": "Observed, 500 m",
        "null_500m_mean": "Random, 500 m", "null_500m_q975": "Random P97.5, 500 m", "p_500m": "p, 500 m",
        "obs_1km": "Observed, 1 km", "null_1km_mean": "Random, 1 km", "null_1km_q975": "Random P97.5, 1 km",
        "p_1km": "p, 1 km", "no": "No.", "width_deg": "Width (°)", "temp_median": "Temperature median",
        "temp_extreme": "Temperature extreme", "temp_median_neighbours": "Temperature neighbours",
        "prec_median": "Precipitation median", "lowest_pass": "Lowest pass", "lowest_z_m": "Lowest pass height (m)",
        "control_passes": "Control passes"}
DROP = {"S4": ["lat_dms", "lon_dms", "elev_alt_m", "source", "coord_source", "note"],
        "S19": ["bins", "tiles_min", "tiles_max", "temp_most_limiting", "temp_least_limiting", "prec_extreme",
                "prec_most_limiting", "prec_least_limiting", "prec_median_neighbours", "lowest_lat", "lowest_lon",
                "lowest_depth_m"]}    # full in the workbook
VAR = {f"bio{v:02d}": f"BIO{v}" for v in range(1, 20)}
FIGS = [
    ("S1", "FigS1_validation", "Validation of the census against documented passes. (a) Cumulative share of OpenStreetMap pass nodes by distance "
     "to the nearest pass (black), and the mean share of random divide cells of the same block core and unit within "
     "500 m and 1 km (open circles). (b) Height of the nearest census pass minus the official height of the pass, tropical "
     "(filled circles) and temperate (open circles) references; the grey band is 50 m. Tables S3 to S5."),
    ("S2", "FigS2_variables_transverse", "Normalized current for transverse connectivity for each of the 19 WorldClim 2.1 "
     "variables, main configuration. Colour scale as in Fig. 2."),
    ("S3", "FigS3_variables_longitudinal", "Normalized current for longitudinal connectivity for each of the 19 WorldClim "
     "2.1 variables, main configuration. Colour scale as in Fig. 2."),
    ("S4", "FigS4_sensitivity_maps", "Normalized current with BIO1 in the four sensitivity configurations, both "
     "directions: margins q = 0.05 and 0.2 with tiles of 200 km, and tiles of 100 and 400 km with q = 0.1. Colour "
     "scale as in Fig. 2; uncoloured cells have no value in that configuration. Triangles are the control passes of "
     "Table 2."),
    ("S5", "FigS5_sensitivity_profiles", "Median relative conductance of the tiles with BIO1 per 0.5° band in the five "
     "configurations, both directions. The main configuration is in black; solid grey lines are the margins q = 0.05 "
     "(dark) and 0.2 (light), and dashed grey lines are the tiles of 100 km (dark) and 400 km (light). The dotted "
     "line is uniform resistance. Every tile enters, as in Table S8; Table S9 gives the same profile over the tiles "
     "that cross the divide, in the two control zones."),
    ("S6", "FigS6_robust_passes", "Cells with passes by the median normalized current over the 19 variables: above 1 in "
     "all five configurations (red), below 1 in all five (blue), or depending on the configuration (grey); one pass "
     "per 1 km cell, the deepest. Tables S13 and S16."),
]
DICT = [
    ("apd_id", "identifier of the pass in the Andean Passes Database"),
    ("lat, lon", "position of the pass (degrees, WGS84)"),
    ("z_m", "height of the pass in FABDEM V1-2 (m)"),
    ("depth_m", "depth by prominence along the divide (m)"),
    ("upper_max_m", "height of the higher of the two points that enclose the pass (m)"),
    ("width_m", "length of divide below the pass height plus 100 m (m); passes with depth of 100 m or more"),
    ("ring, flat, flat_km2", "pass on the rim of a real sink; flat pass and area of its flat (km²)"),
    ("depth_lower_bound", "1 when the lower side reaches the edge of the processing window"),
    ("lat_a ... lon_b5", "first cells of each side, and the points where each side drains 5 km²"),
    ("gmba", "GMBA v2 unit"),
    ("sea_a, sea_b, rank", "sea or interior reached by each side, and rank (continental, terminal, tributary)"),
    ("hybas_term_a, hybas_term_b, area_sep_km2", "terminal HydroBASINS basin of each side and separated area (km²)"),
    ("isolation_km", "distance to the nearest lower pass with depth of 100 m or more (km)"),
    ("tropics, cell_1km", "1 for passes between 11°N and 23.5°S; 1 km cell of the connectivity grid"),
    ("trans_med, longi_med", "median normalized current over the 19 variables, transverse and longitudinal, main "
     "configuration"),
    ("trans_min, trans_max, longi_min, longi_max", "range over the 19 variables, main configuration"),
    ("*_med_q005, *_med_q02, *_med_100km, *_med_400km", "median in each sensitivity configuration"),
    ("trans_robust, longi_robust", "above 1 or below 1 in all five configurations, or depends"),
]

doc = Document()
st = doc.styles["Normal"]
st.font.name, st.font.size = "Times New Roman", Pt(10)
doc.add_paragraph().add_run("Supplementary Information").bold = True
doc.add_paragraph("A census of Andean passes and their climatic connectivity across and along the tropical Andes. Tables S1 to S16 are given in "
                  "full in L1_Supplementary_Tables.xlsx; tables longer than 60 rows appear there only.")
wb = pd.ExcelWriter(f"{F}/L1_Supplementary_Tables.xlsx", engine="openpyxl")
pd.DataFrame([(n, c) for n, _, c in TABLES], columns=["table", "caption"]).to_excel(wb, sheet_name="contents", index=False)
for n, f, cap in TABLES:
    d = pd.read_csv(f"{T}/{f}.csv")
    d.to_excel(wb, sheet_name=f"Table {n}", index=False)
    p = doc.add_paragraph()
    p.add_run(f"Table {n}. ").bold = True
    p.add_run(cap)
    d = d.drop(columns=DROP.get(n, []))
    if "var" in d:
        d["var"] = d["var"].map(VAR)
    if len(d) > INLINE_MAX:
        doc.add_paragraph(f"{len(d)} rows; in the workbook, sheet Table {n}.")
        continue
    t = doc.add_table(rows=1, cols=len(d.columns))
    t.style = "Table Grid"
    for c, h in zip(t.rows[0].cells, d.columns):
        c.text = HEAD.get(h, h)
    for _, r in d.iterrows():
        for c, v in zip(t.add_row().cells, r):
            c.text = ("" if pd.isna(v) else f"{v:,.0f}".replace(",", " ") if isinstance(v, float) and abs(v) >= 1000
                      else f"{v:.3g}" if isinstance(v, float) else str(v).replace("_", " "))
    for row in t.rows:
        for c in row.cells:
            for q in c.paragraphs:
                for x in q.runs:
                    x.font.size = Pt(6.5)
    doc.add_paragraph()
pd.DataFrame(DICT, columns=["column", "meaning"]).to_excel(wb, sheet_name="APD dictionary", index=False)
wb.close()
for n, f, cap in FIGS:
    doc.add_picture(f"{F}/{f}.png", width=Cm(16))
    p = doc.add_paragraph()
    p.add_run(f"Figure {n}. ").bold = True
    p.add_run(cap)
p = doc.add_paragraph()
p.add_run("Andean Passes Database (outputs/database/APD_passes.csv.gz). ").bold = True
p.add_run("One row per pass of the census, with the columns below.")
t = doc.add_table(rows=0, cols=2)
t.style = "Table Grid"
for k, v in DICT:
    a, b = t.add_row().cells
    a.text, b.text = k, v
doc.save(f"{F}/L1_Supplementary_Information.docx")
print("->", F)
