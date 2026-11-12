"""Fill-in templates for the project workflow (Colab box 3 / 5 / 6) and a checker for filled-in lists.

Templates (CSV, one header row, example rows to overwrite):
  sample_sheet_template.csv      file, type, group, replicate, cal_level, notes
  calibration_list_template.csv  name, formula, monoisotopic_mass, smiles, label, role, is_name, adduct, rt_expected, units,
                                 L1 ... L8 (concentration of THIS compound at each calibration level)
  is_list_template.csv           name, formula, monoisotopic_mass, smiles, label, role, conc, units, adduct, polarity,
                                 rt_expected, notes
Sample types: sample | blank | qc | cal | is_mix | solvent_std | matrix_spike | exclude
  is_mix        IS (and/or target) mix in solvent - response factors for the IS-surrogate semi-quant
  solvent_std   standard mix in neat solvent (e.g. 10 % MeOH / water)         } matrix effect = matrix_spike / solvent_std
  matrix_spike  the same mix spiked into a control (blank) matrix EXTRACT      } (calculated after the run, box 6b)
"""
from __future__ import annotations
import os
import re
import numpy as np
import pandas as pd

SAMPLE_TYPES = ("sample", "blank", "qc", "cal", "is_mix", "solvent_std", "matrix_spike", "exclude")
LEVELS = [f"L{i}" for i in range(1, 9)]

SAMPLE_SHEET = pd.DataFrame([
    dict(file="Blank_01", type="blank", group="procedural_blank", replicate="1", cal_level="", notes="procedural / solvent blank"),
    dict(file="Cal_L1", type="cal", group="cal", replicate="1", cal_level="L1", notes="lowest calibration level"),
    dict(file="Cal_L2", type="cal", group="cal", replicate="1", cal_level="L2", notes=""),
    dict(file="Cal_L3", type="cal", group="cal", replicate="1", cal_level="L3", notes=""),
    dict(file="Cal_L4", type="cal", group="cal", replicate="1", cal_level="L4", notes=""),
    dict(file="Cal_L5", type="cal", group="cal", replicate="1", cal_level="L5", notes="highest calibration level"),
    dict(file="ISmix_solvent_01", type="is_mix", group="is_mix", replicate="1", cal_level="", notes="IS mix in solvent: response factors"),
    dict(file="Std_10pcMeOH_01", type="solvent_std", group="matrix_test", replicate="1", cal_level="", notes="mix in 10 % MeOH / water"),
    dict(file="Std_spikedExtract_01", type="matrix_spike", group="matrix_test", replicate="1", cal_level="", notes="same mix in control extract"),
    dict(file="QC_pool_01", type="qc", group="QC", replicate="1", cal_level="", notes="pooled QC, inject every ~10 samples"),
    dict(file="SiteA_rep1", type="sample", group="SiteA", replicate="1", cal_level="", notes=""),
    dict(file="SiteA_rep2", type="sample", group="SiteA", replicate="2", cal_level="", notes=""),
    dict(file="SiteA_rep3", type="sample", group="SiteA", replicate="3", cal_level="", notes=""),
])

CAL_LIST = pd.DataFrame([
    dict(name="Carbamazepine", formula="C15H12N2O", monoisotopic_mass=236.09496, smiles="NC(=O)N1c2ccccc2C=Cc2ccccc12", label="",
         role="target", is_name="Carbamazepine-D10", adduct="[M+H]+", rt_expected="", units="ng/mL",
         L1=0.5, L2=1, L3=5, L4=10, L5=50, L6=100, L7="", L8=""),
    dict(name="Venlafaxine", formula="C17H27NO2", monoisotopic_mass=277.20418, smiles="COc1ccc(cc1)C(CN(C)C)C1(O)CCCCC1", label="",
         role="target", is_name="Venlafaxine-D6", adduct="[M+H]+", rt_expected="", units="ng/mL",
         L1=0.5, L2=1, L3=5, L4=10, L5=50, L6=100, L7="", L8=""),
    dict(name="Carbamazepine-D10", formula="C15H2[2H]10N2O", monoisotopic_mass=246.15773, smiles="NC(=O)N1c2ccccc2C=Cc2ccccc12",
         label="D10", role="IS", is_name="", adduct="[M+H]+", rt_expected="", units="ng/mL",
         L1=10, L2=10, L3=10, L4=10, L5=10, L6=10, L7="", L8=""),
    dict(name="Venlafaxine-D6", formula="C17H21[2H]6NO2", monoisotopic_mass=283.24184, smiles="COc1ccc(cc1)C(CN(C)C)C1(O)CCCCC1",
         label="D6", role="IS", is_name="", adduct="[M+H]+", rt_expected="", units="ng/mL",
         L1=10, L2=10, L3=10, L4=10, L5=10, L6=10, L7="", L8=""),
])

IS_LIST = pd.DataFrame([
    dict(name="Carbamazepine-D10", formula="C15H2[2H]10N2O", monoisotopic_mass=246.15773, smiles="NC(=O)N1c2ccccc2C=Cc2ccccc12",
         label="D10", role="IS", conc=0.1, units="ug/mL", adduct="[M+H]+", polarity="positive", rt_expected="", notes=""),
    dict(name="Diazepam-D5", formula="C16H8[2H]5ClN2O", monoisotopic_mass=289.10302, smiles="CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
         label="D5", role="IS", conc=0.1, units="ug/mL", adduct="[M+H]+", polarity="positive", rt_expected="", notes=""),
    dict(name="Diclofenac-D4", formula="C14H7[2H]4Cl2NO2", monoisotopic_mass=299.04179, smiles="OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
         label="D4", role="IS", conc=0.1, units="ug/mL", adduct="[M-H]-", polarity="negative", rt_expected="", notes="ESI-"),
])

README = """HOW TO FILL THESE TEMPLATES (delete the example rows, keep the header row; save as CSV)

sample_sheet_template.csv - one row per injection
  file       mzML/raw file name without extension (must match the data folder)
  type       sample | blank | qc | cal | is_mix | solvent_std | matrix_spike | exclude
               is_mix       = internal-standard mix in solvent (response factors for semi-quantification)
               solvent_std  = standard / IS mix in 10 % MeOH / water
               matrix_spike = the SAME mix spiked into a control (blank) matrix extract
               (matrix effect = matrix_spike / solvent_std, calculated after the run)
  group      replicate injections of one sample share a group (the replicate rule, e.g. 2/3, is applied per group)
  replicate  1, 2, 3 ...
  cal_level  for type cal only: L1 ... L8, matching the level columns of the calibration list

calibration_list_template.csv - one row per calibrated compound and per internal standard
  formula            neutral molecular formula; isotope labels written as [2H] / [13C] / [15N]
  monoisotopic_mass  neutral monoisotopic mass, 5 decimal places (checked against the formula)
  role               target | IS
  is_name            for targets: the IS the response is normalised to (must be a row with role IS)
  L1 ... L8          concentration of THIS compound at each calibration level (same units for the row)

is_list_template.csv - one row per internal standard (or target) to screen in every injection
  conc / units       concentration in the IS mix (used for response factors / semi-quantification)
  adduct, polarity   expected ion, e.g. [M+H]+ / positive or [M-H]- / negative
"""


def write_templates(out_dir: str) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    for name, df in (("sample_sheet_template.csv", SAMPLE_SHEET), ("calibration_list_template.csv", CAL_LIST),
                     ("is_list_template.csv", IS_LIST)):
        p = os.path.join(out_dir, name)
        df.to_csv(p, index=False)
        paths[name] = p
    p = os.path.join(out_dir, "README_templates.txt")
    open(p, "w", encoding="utf-8").write(README)
    paths["README_templates.txt"] = p
    return paths


def _read(path_or_df):
    if isinstance(path_or_df, pd.DataFrame):
        return path_or_df.copy()
    p = str(path_or_df)
    return pd.read_excel(p) if p.lower().endswith((".xlsx", ".xls")) else pd.read_csv(p)


def check_compound_list(path_or_df, kind="is") -> pd.DataFrame:
    """Row-by-row check of a filled-in IS / calibration list. Returns the table with 'check' and 'mass_from_formula'
    columns: formula parses, monoisotopic mass within 5 ppm of the formula (5 decimal places expected), role valid,
    is_name points to an IS row, levels numeric (calibration), conc present (IS list)."""
    from .targeted import parse_formula, mono
    df = _read(path_or_df)
    df.columns = [str(c).strip() for c in df.columns]
    low = {c.lower(): c for c in df.columns}
    out = []
    names_is = set(df.loc[df[low["role"]].astype(str).str.upper() == "IS", low["name"]]) if "role" in low and "name" in low else set()
    lv = [c for c in df.columns if re.fullmatch(r"L\d+", c)]
    for _, r in df.iterrows():
        msg = []
        f = str(r.get(low.get("formula", ""), "") or "").strip()
        m_given = pd.to_numeric(r.get(low.get("monoisotopic_mass", ""), np.nan), errors="coerce")
        m_calc = np.nan
        if f and f.lower() != "nan":
            try:
                m_calc = mono(parse_formula(f))
            except Exception:
                msg.append(f"formula '{f}' not understood")
        elif not str(r.get(low.get("smiles", ""), "") or "").strip():
            msg.append("needs a formula or a SMILES")
        if np.isfinite(m_given) and np.isfinite(m_calc):
            ppm = (m_given - m_calc) / m_calc * 1e6
            if abs(ppm) > 5:
                msg.append(f"monoisotopic mass {m_given:.5f} differs from formula ({m_calc:.5f}) by {ppm:+.1f} ppm")
            s = f"{m_given}"
            if "." in s and len(s.split(".")[1]) < 4:
                msg.append("give the monoisotopic mass to 5 decimal places")
        role = str(r.get(low.get("role", ""), "") or "").strip()
        if role and role.upper() not in ("IS", "TARGET"):
            msg.append(f"role '{role}' should be IS or target")
        if kind == "cal":
            isn = str(r.get(low.get("is_name", ""), "") or "").strip()
            if role.lower() == "target" and isn and isn.lower() != "nan" and isn not in names_is:
                msg.append(f"is_name '{isn}' is not a row with role IS")
            vals = pd.to_numeric(r[lv], errors="coerce") if lv else pd.Series(dtype=float)
            if lv and vals.notna().sum() < 3 and role.lower() == "target":
                msg.append(f"only {int(vals.notna().sum())} calibration levels (>= 3 needed, >= 5 recommended)")
        if kind == "is" and "conc" in low and not np.isfinite(pd.to_numeric(r.get(low["conc"]), errors="coerce")):
            msg.append("no concentration (needed for response factors)")
        out.append(dict(mass_from_formula=round(m_calc, 5) if np.isfinite(m_calc) else np.nan, check="; ".join(msg) or "OK"))
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(out)], axis=1)


def level_concentrations(cal_list: pd.DataFrame) -> pd.DataFrame:
    """Long table name, cal_level, conc from the L1..L8 columns of a calibration list."""
    lv = [c for c in cal_list.columns if re.fullmatch(r"L\d+", str(c))]
    if not lv:
        return pd.DataFrame(columns=["name", "cal_level", "conc"])
    long = cal_list.melt(id_vars=["name"], value_vars=lv, var_name="cal_level", value_name="conc")
    long["conc"] = pd.to_numeric(long["conc"], errors="coerce")
    return long.dropna(subset=["conc"])
