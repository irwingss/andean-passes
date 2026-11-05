"""34_build_supp.py - Supporting Information of Article 1 from the outputs of 33_suplementos.R.
Writes outputs/supp/L1_Supplementary_Tables.xlsx (Contents sheet plus one sheet per table, title and note above it) and
outputs/supp/L1_Supplementary_Information.docx (text of S1-S5, Figs S1.1, S1.2, S2, S5, the S4 atlas and the tables that fit a page;
S2.4, S4.1, S5.2, S5.3 and S5.4 are in the workbook only). Same look as manuscript/fuente/build.py.
python scripts/34_build_supp.py <project root>"""
import csv, os, sys, docx, openpyxl
from docx.enum.section import WD_ORIENT
from docx.enum.text import WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Mm, Pt
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
SUPP = ROOT + '/outputs/supp'
rd = lambda f: list(csv.reader(open(f, encoding='utf-8-sig')))
contents = {r['sheet']: r for r in csv.DictReader(open(SUPP + '/tables/contents.csv', encoding='utf-8-sig'))}
tables = {k: rd(f"{SUPP}/tables/{k}.csv") for k in contents}

def num(v):
    try: return int(v) if v.lstrip('-').isdigit() else float(v)
    except ValueError: return {'TRUE': 'yes', 'FALSE': 'no', 'NA': ''}.get(v, v)

# ---- Workbook ----
wb = openpyxl.Workbook(); ws = wb.active; ws.title = 'Contents'
ws.append(['Sheet', 'Title', 'Note'])
for k, c in contents.items(): ws.append([k, c['title'], c['note']])
for col, w in zip('ABC', (8, 70, 120)): ws.column_dimensions[col].width = w
for row in ws.iter_rows(min_row=2):
    for c in row: c.alignment = Alignment(wrap_text=True, vertical='top')
for c in ws[1]: c.font = Font(bold=True)
for k, rows in tables.items():
    s = wb.create_sheet(k)
    s.append([contents[k]['title']]); s['A1'].font = Font(bold=True)
    s.append([contents[k]['note']]); s.append([])
    s.append(rows[0])
    for c in s[4]: c.font = Font(bold=True); c.alignment = Alignment(wrap_text=True, vertical='bottom')
    for r in rows[1:]: s.append([num(v) for v in r])
    for j, h in enumerate(rows[0], 1):
        width = max([len(h) / 2] + [len(str(r[j - 1])) for r in rows[1:]])
        s.column_dimensions[get_column_letter(j)].width = min(max(8, width + 2), 45)
    s.freeze_panes = 'A5'
wb.save(SUPP + '/L1_Supplementary_Tables.xlsx')

# ---- Word ----
d = docx.Document()
st = d.styles['Normal']; st.font.name = 'Times New Roman'; st.font.size = Pt(12)
st.paragraph_format.space_after = Pt(6); st.paragraph_format.line_spacing = 1.5

def P(text='', lead=None, bold=False, size=None):
    p = d.add_paragraph()
    if lead: p.add_run(lead).bold = True
    r = p.add_run(text); r.bold = bold
    if size:
        for r in p.runs: r.font.size = Pt(size)
    return p
H = lambda t: P(t, bold=True)

def rule(cell, edge):
    tcPr = cell._tc.get_or_add_tcPr()
    b = tcPr.find(qn('w:tcBorders'))
    if b is None: b = OxmlElement('w:tcBorders'); tcPr.append(b)
    e = OxmlElement('w:' + edge)
    for k, v in (('val', 'single'), ('sz', '6'), ('space', '0'), ('color', '000000')): e.set(qn('w:' + k), v)
    b.append(e)

def table(k, size=7):
    c = contents[k]; rows = tables[k]
    P(c['title'])
    t = d.add_table(rows=len(rows), cols=len(rows[0]))
    for i, r in enumerate(rows):
        for j, v in enumerate(r):
            cell = t.rows[i].cells[j]; cell.text = str(num(v)) if i else v
            if i == 0: rule(cell, 'top'); rule(cell, 'bottom')
            if i == len(rows) - 1: rule(cell, 'bottom')
    for i, row in enumerate(t.rows):
        for cell in row.cells:
            for p in cell.paragraphs:
                p.paragraph_format.space_after = Pt(0); p.paragraph_format.line_spacing = 1.0
                for r in p.runs: r.font.size = Pt(size); r.bold = i == 0
    if c['note']: P(c['note'], size=9)

def section(landscape):
    s = d.add_section()
    w, h = (Mm(297), Mm(210)) if landscape else (Mm(210), Mm(297))
    s.orientation = WD_ORIENT.LANDSCAPE if landscape else WD_ORIENT.PORTRAIT
    s.page_width, s.page_height = w, h
    for m in ('left_margin', 'right_margin', 'top_margin', 'bottom_margin'): setattr(s, m, Mm(20))

def figure(path, caption, lead, width=165):
    d.add_picture(path, width=Mm(width)); P(caption, lead=lead)

s0 = d.sections[0]; s0.page_width, s0.page_height = Mm(210), Mm(297)
P('Supporting Information', bold=True, size=14)
P('Climatic Permeability of Andean Mountain passes: a continental census under present and future climates')
P('Contents: S1, census algorithm and validation. S2, sensitivity of the census and of the metrics. S3, sensitivity of the '
  'typology of tropical gateways. S4, atlas of the 30 tropical gateways. S5, selection of climate models and results by model. '
  'All tables are also in L1_Supplementary_Tables.xlsx; Tables S2.4, S4.1, S5.2, S5.3 and S5.4 are in the workbook only.')

# S1
H('S1. Census algorithm and validation')
for i, s in enumerate([
    'Label each level-12 basin of HydroBASINS as western (Pacific and the Caribbean basins west of the Eastern Cordillera of '
    'Colombia and of the Mérida Andes), eastern (Orinoco, Amazon, Río de la Plata and the Atlantic slopes of Patagonia) or '
    'interior (endorheic basins connected to the eastern domain; coastal endorheic basins count as western) from the position '
    'of its outlet (Figure S1.1).',
    'Place western and eastern seeds on Copernicus GLO-30 (1 arc-second) more than 12 km from any contact between sides; '
    'interior basins carry no seed.',
    'Flood the elevation model from both sets of seeds by priority flooding (a marker-controlled watershed), in blocks of 1° of '
    'latitude except over the Altiplano–Puna (17–27°S), which is flooded in a single window. The divide is the '
    'boundary where the western and eastern floods meet, and each of its edges carries the crossing height, the lowest '
    'west to east route forced through that edge.',
    'Clip the divide to the Andes polygon of GMBA v2 with a 20 km margin, from Lara (Venezuela) to Tierra del Fuego, and '
    'keep the minimum crossing height of each 0.5 km segment of distance along the divide.',
    'Mark as a col every local minimum of this profile. Depth is the smaller of the highest points reached northwards and '
    'southwards before the profile falls below the col, minus the height of the col. Isolation is the distance along the '
    'divide to the nearest lower point in each direction. Width is the continuous length of divide below the height of the '
    'col plus 100 m.',
    'A gateway is a col with depth of at least 100 m that is the lowest point within 100 km along the divide in both '
    'directions. Cols within 50 km of the ends of the chain, and sea channels below 50 m south of 52°S, are terminal. The '
    'census counts 65 gateways, 31 of them tropical; the four terminal gateways stay out of the analyses, which leaves 61 '
    'gateways, 30 of them tropical.',
    'Validate crossing heights against passes with published elevations, matching each to the nearest col within 5 km.'], 1):
    P(s, lead=f'Step {i}. ')
section(True); figure(ROOT + '/outputs/figs_qc/M1_dominios.png', 'Drainage domains and divide of the census. Fill: level-12 basins of HydroBASINS '
       'classified by outlet as Pacific and Caribbean (W), Atlantic (E) or endorheic (I). Coastal endorheic basins such as Sechura and '
       'Pisco are shown as endorheic but count as western in the flood. Solid line: W|E divide computed on Copernicus GLO-30. Dashed '
       'line: I|E boundary of HydroBASINS. Thin lines: Andes of the GMBA Mountain Inventory v2. Points: the 30 non-terminal tropical '
       'gateways, with their crossing heights. The horizontal dotted line marks 23.5°S. (A) Whole cordillera. (B) Northern Andes. '
       '(C) Central Andes of Peru. (D) Altiplano–Puna plateau.', 'Figure S1.1. ', width=165)
section(True); table('S1.1', size=6)
section(False); figure(SUPP + '/FigS1_validation.png', 'Crossing height of the matched col (Copernicus GLO-30) against the published elevation '
       'for the 18 documented passes with a col within 5 km. The dashed line marks equal values; labelled passes have an '
       'absolute error above 50 m. The source of each published elevation is in Tables S1.2 and S1.3.', 'Figure S1.2. ', width=110)
section(True); table('S1.2', size=6); table('S1.3', size=8); table('S1.4', size=7); table('S1.5', size=7)

# S2
section(True); H('S2. Sensitivity of the census and of the metrics')
table('S2.1'); table('S2.2')
figure(SUPP + '/FigS2_janzen_sensitivity.png', 'Interaction between crossing height and absolute latitude in the GAM of '
       'thermal overlap, with its 95% confidence interval, for the main model and the sensitivity analyses of Table S2.2. '
       "A positive interaction is Janzen's prediction. The two indices on the scale of 12 months are in Table S2.2 only.",
       'Figure S2. ')
table('S2.3'); table('S2.5')

# S3 (wide table)
section(True); H('S3. Sensitivity of the typology of tropical gateways')
table('S3.1', size=6); table('S3.2'); table('S3.3')

# S4 atlas
section(False); H('S4. Atlas of the 30 tropical gateways')
P('Each page shows one tropical gateway, from north to south; their attributes are in Table S4.1 (workbook). (a) Relief '
  '(HydroSHEDS 15 arc-seconds) with the drainage divide (white line), the col (triangle), the flank points at 1000 m on each '
  'side and the steepest-descent routes that reached them (red lines); a flank without a route is the nearest cell below '
  '1000 m on that side. (b) Crossing height along the divide within 150 km of the col; the dashed line marks 2000 m, the '
  'lower limit of montane climate. (c) Annual temperature range, from the minimum temperature of the coldest month to the '
  'maximum of the warmest month, of each flank and of the crest at present, and of the crest in 2081–2100 under three '
  'scenarios (mean of eight CMIP6 models).')
atlas = rd(SUPP + '/tables/S4.1.csv'); hd = atlas[0]
for r in atlas[1:]:
    a = dict(zip(hd, r))
    d.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
    figure(f"{SUPP}/atlas/{a['File']}", f"{a['Gateway']}, {a['Country']} ({a['Latitude']}°, {a['Crossing height (m)']} m). "
           f"Type: {a['Type'].lower()}.", f"Figure S4.{a['No']}. ")

# S5
section(False); H('S5. Selection of climate models and results by model')
table('S5.1', size=8)
GAP0 = next(round(float(r['gap_km'])) for r in csv.DictReader(open('outputs/tables/M3_corredor_montano.csv'))
            if (r['threshold_m'], r['zone'], r['gcm']) == ('2000', 'tropical', 'present'))
figure(SUPP + '/FigS5_gcm.png', '(a) Median fraction of ΔT* reached by the 30 tropical gateways and (b) length of tropical '
       'divide in montane gaps (threshold 2000 m), for each of the eight CMIP6 models by scenario and period. The dashed '
       f'line in b is the present length ({GAP0} km). Values by model are in Tables S5.2 and S5.3.', 'Figure S5. ')

out = SUPP + '/L1_Supplementary_Information.docx'
table('S5.5', size=8); table('S5.6', size=8); table('S5.7', size=8); table('S5.8', size=8)
d.save(out); print('ok', out, SUPP + '/L1_Supplementary_Tables.xlsx')
