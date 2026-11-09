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
    ("S1", "S1_census_by_unit", "Cols of the census by zone, GMBA unit and rank: number of cols with a depth of 30 m or "
     "more, number with a depth of 100 m or more, and median height and depth (m). Tropics are 11°N to 23.5°S."),
    ("S2", "S2_control_passes_census", "The 36 control passes inside the area: the nearest col (near_) and the nearest "
     "col with a depth of 100 m or more (deep_), with distance (km), height (m), depth (m), rank, separated area (km²), width (m) and "
     "isolation (km). Measured from relief alone."),
    ("S3", "S3_osm_nodes_identity", "OpenStreetMap pass nodes (mountain_pass=yes or natural=saddle, downloaded 8 "
     "November 2026): share with a col within 500 m and 1 km by unit and rank (nodes_x), and the identity test for "
     "nodes on the divide network (nodes_y): observed share against the mean and 97.5th percentile of 1000 draws of "
     "random divide cells of the same block core and unit, with the empirical p value."),
    ("S4", "S4_reference_heights", "Official heights of passes (INEI, Peru; UPF, Chile), the nearest col and its error "
     "(col minus reference, m)."),
    ("S5", "S5_height_error_by_zone", "Error of the height of the nearest col by zone, against official heights and "
     "against the ele tag of OpenStreetMap (which may itself come from an elevation model)."),
    ("S6", "S6_osm_nodes_without_col", "OpenStreetMap pass nodes without a col within 1 km."),
    ("S7", "S7_envelopes", "Montane limit T2000(lat) (°C), interquartile range of each variable over the domain, and the "
     "5th and 95th percentiles of each variable in the lowland and montane guilds, every 0.5° of latitude. "
     "Precipitation variables as log(1 + P)."),
    ("S8", "S8_tile_conductance_profiles", "Relative conductance of the tiles per configuration, reading, variable and "
     "0.5° band of latitude of the tile centre: number of tiles, median and quartiles."),
    ("S9", "S9_control_huancabamba_tachira", "Control: relative conductance of the tiles in the longitudinal reading "
     "with BIO1, per configuration and 0.5° band, around the Huancabamba zone (2 to 10°S) and the Táchira depression "
     "(74 to 70°W, north of 5°N). in_zone marks the bands of each control zone."),
    ("S10", "S10_spearman_variables", "Spearman correlation between pairs of variables in the normalized current of the "
     "cols (one col per 1 km cell) and in the relative conductance of the tiles, main configuration."),
    ("S11", "S11_variable_agreement", "Agreement of each variable with the others: median and minimum Spearman "
     "correlation with the other 18 variables at the level of cols, median at the level of tiles, and correlation "
     "with the col median over the 19 variables. Main configuration."),
    ("S12", "S12_sensitivity", "Sensitivity of the col medians against the main configuration: cells with a value in "
     "both, cells that lose their value, Spearman correlation, share of the top decile that stays in the top decile, "
     "and share of cells on the same side of 1."),
    ("S13", "S13_by_unit", "Col medians by GMBA unit and reading: cells, median and 90th percentile of the median over "
     "the 19 variables, and number and share of cells that concentrate (above 1) or divert (below 1) current in all "
     "five configurations."),
    ("S14", "S14_by_latitude", "Col medians by band of 1° of latitude and reading: cells, median, 90th percentile and "
     "share of cells above 1 in all five configurations."),
    ("S15", "S15_controls_by_variable", "Mean of the band medians of tile relative conductance inside each control zone "
     "(Huancabamba 4 to 6°S; Táchira 7.5 to 8°N) and in the latitudes next to it (north 2 to 4°S, south 6.5 to "
     "8.5°S; Táchira south 5 to 7.5°N), per variable and reading, main configuration."),
]
INLINE_MAX = 60       # longer tables are given only in the workbook
FIGS = [
    ("S1", "FigS1_validation", "Validation of the census. (a) Cumulative share of OpenStreetMap pass nodes by distance "
     "to the nearest col (black), and the mean share of random divide cells of the same block core and unit within "
     "500 m and 1 km (blue). (b) Height of the nearest col minus the official height of the pass, tropical and "
     "temperate references; the grey band is 50 m. Tables S3 to S5."),
    ("S2", "FigS2_variables_transverse", "Normalized current in the transverse reading for each of the 19 WorldClim 2.1 "
     "variables, main configuration. Colour scale as in Fig. 2."),
    ("S3", "FigS3_variables_longitudinal", "Normalized current in the longitudinal reading for each of the 19 WorldClim "
     "2.1 variables, main configuration. Colour scale as in Fig. 2."),
    ("S4", "FigS4_sensitivity_maps", "Normalized current with BIO1 in the four sensitivity configurations, both "
     "readings: margins q = 0.05 and 0.2 with tiles of 200 km, and tiles of 100 and 400 km with q = 0.1."),
    ("S5", "FigS5_sensitivity_profiles", "Median relative conductance of the tiles with BIO1 per 0.5° band in the five "
     "configurations, both readings. Blue bands mark the Huancabamba zone and the Táchira depression. Table S9."),
    ("S6", "FigS6_robust_cols", "Cells with cols by the median normalized current over the 19 variables: above 1 in "
     "all five configurations (red), below 1 in all five (blue), or depending on the configuration (grey). Table S13."),
]
DICT = [
    ("apd_id", "identifier of the col in the Andean Passes Database"),
    ("lat, lon", "position of the col (degrees, WGS84)"),
    ("z_m", "height of the col in FABDEM V1-2 (m)"),
    ("depth_m", "depth by prominence along the divide (m)"),
    ("upper_max_m", "height of the higher of the two points that enclose the col (m)"),
    ("width_m", "length of divide below the col height plus 100 m (m); cols with depth of 100 m or more"),
    ("ring, flat, flat_km2", "col on the rim of a real sink; flat col and area of its flat (km²)"),
    ("depth_lower_bound", "1 when the lower side reaches the edge of the processing window"),
    ("lat_a ... lon_b5", "first cells of each side, and the points where each side drains 5 km²"),
    ("gmba", "GMBA v2 unit"),
    ("sea_a, sea_b, rank", "sea or interior reached by each side, and rank (continental, terminal, tributary)"),
    ("hybas_term_a, hybas_term_b, area_sep_km2", "terminal HydroBASINS basin of each side and separated area (km²)"),
    ("isolation_km", "distance to the nearest lower col with depth of 100 m or more (km)"),
    ("tropics, cell_1km", "1 for cols between 11°N and 23.5°S; 1 km cell of the connectivity grid"),
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
doc.add_paragraph("Census and climatic connectivity of the cols of the tropical Andes. Tables S1 to S15 are given in "
                  "full in L1_Supplementary_Tables.xlsx; tables longer than 60 rows appear there only.")
wb = pd.ExcelWriter(f"{F}/L1_Supplementary_Tables.xlsx", engine="openpyxl")
pd.DataFrame([(n, c) for n, _, c in TABLES], columns=["table", "caption"]).to_excel(wb, sheet_name="contents", index=False)
for n, f, cap in TABLES:
    d = pd.read_csv(f"{T}/{f}.csv")
    d.to_excel(wb, sheet_name=f"Table {n}", index=False)
    p = doc.add_paragraph()
    p.add_run(f"Table {n}. ").bold = True
    p.add_run(cap)
    if len(d) > INLINE_MAX:
        doc.add_paragraph(f"{len(d)} rows; in the workbook, sheet Table {n}.")
        continue
    t = doc.add_table(rows=1, cols=len(d.columns))
    t.style = "Table Grid"
    for c, h in zip(t.rows[0].cells, d.columns):
        c.text = str(h)
    for _, r in d.iterrows():
        for c, v in zip(t.add_row().cells, r):
            c.text = "" if pd.isna(v) else (f"{v:.3g}" if isinstance(v, float) else str(v))
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
p.add_run("Andean Passes Database (outputs/database/APD_cols.csv.gz). ").bold = True
p.add_run("One row per col of the census, with the columns below.")
t = doc.add_table(rows=0, cols=2)
t.style = "Table Grid"
for k, v in DICT:
    a, b = t.add_row().cells
    a.text, b.text = k, v
doc.save(f"{F}/L1_Supplementary_Information.docx")
print("->", F)
