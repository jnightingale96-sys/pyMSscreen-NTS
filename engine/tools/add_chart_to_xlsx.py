"""Add the rebuilt PEC-vs-MEC chart and its supporting tables to work/grain_PEC_MEC_semiquant_correction.xlsx."""
import os
import pandas as pd
from openpyxl import load_workbook
from openpyxl.drawing.image import Image
from openpyxl.styles import Font

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); W = os.path.join(ROOT, "work")
X = os.path.join(W, "grain_PEC_MEC_semiquant_correction.xlsx")
KEEP = ["Caffeine", "Carbamazepine", "DIOH-CBZ", "Venlafaxine"]
O = pd.read_csv(os.path.join(W, "grain_pec_mec_by_field.csv"), dtype={"field": str})
D = O[O.chemical.isin(KEEP) & (O.MEC > 0)][["chemical", "field", "MEC", "targeted_ng_L", "corrected_ng_L", "SQ_ng_L", "PEC_targ", "PEC_corr", "PEC_SQ", "ratio_targ", "ratio_corr", "ratio_SQ"]]
D = D.rename(columns={"MEC": "MEC_ng_g_dw", "PEC_targ": "PEC_targeted_water", "PEC_corr": "PEC_NTS_corrected", "PEC_SQ": "PEC_NTS_uncorrected",
                      "ratio_targ": "PEC/MEC targeted", "ratio_corr": "PEC/MEC NTS corrected", "ratio_SQ": "PEC/MEC NTS uncorrected"})
with pd.ExcelWriter(X, engine="openpyxl", mode="a", if_sheet_exists="replace") as xw:
    D.to_excel(xw, sheet_name="Chart_data", index=False)
    pd.read_csv(os.path.join(W, "nts_only_ww.csv")).to_excel(xw, sheet_name="NPN_WW_areas", index=False)
    pd.read_csv(os.path.join(W, "nts_only_npn_bound.csv")).to_excel(xw, sheet_name="NPN_PEC_bound", index=False)
wb = load_workbook(X)
if "Chart_PEC_vs_MEC" in wb.sheetnames: del wb["Chart_PEC_vs_MEC"]
ws = wb.create_sheet("Chart_PEC_vs_MEC", 0)
ws["A1"] = "Maize grain: PEC vs MEC - chemicals with an assessable irrigation-water signal"; ws["A1"].font = Font(bold=True, size=13)
img = Image(os.path.join(W, "grain_PEC_vs_MEC_rebuilt.png")); img.width, img.height = 1150, 640; ws.add_image(img, "A3")
notes = ["Factor difference = geometric mean of PEC/MEC over fields (x = PEC above MEC, / = PEC below). Data: sheet Chart_data.",
         "Excluded: gabapentin (zwitterion), thiabendazole (MEC not water-driven), crotamiton (unconfirmed), alprazolam/diltiazem (no water signal),",
         "DEHP / octyl hydrogen phthalate / N,N-diethyldodecanamide (blank level), 17b-estradiol (8.2 min peak is not estradiol by MS2), dodecyl sulphate (cannot map to fields).",
         "NPN: PEC is an upper bound over logKow 4-5, Kaw 1e-8..1e-4, DT50 30-1000 d (sheet NPN_PEC_bound); grain MEC not settled (CD 19-28 vs re-integrated 2-3 ng/g).",
         "Venlafaxine rests on n=1-2 fields."]
for i, t in enumerate(notes): ws.cell(row=36 + i, column=1, value=t)
wb.save(X); print(wb.sheetnames)
