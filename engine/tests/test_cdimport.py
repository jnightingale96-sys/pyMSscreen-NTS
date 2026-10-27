"""cd-screen: Compound Discoverer sheet parsing, structure resolution (aligned rows only), water broadcast, RT identity check."""
import copy
import numpy as np, pandas as pd, pytest
from conftest import ROOT
from hrms_nts import cdimport as CD
from hrms_nts import pecreport as PR


def _ww(pol, name, rt, area, ngl, smiles="C", ik="AAAAAAAAAAAAAA-UHFFFAOYSA-N", inp=None):
    r = [None] * 30
    r[0], r[1], r[2], r[4], r[9], r[10], r[11] = pol, name, "C1H4", rt, area, ngl / 120000, ngl
    r[13], r[17], r[18], r[20], r[23] = inp or name, "1-1-1", ik, smiles, smiles
    return r


def _pl(pol, name, rt, ngg, field, sample="1-3", inp=None, smiles="C", ik="BBBBBBBBBBBBBB-UHFFFAOYSA-N"):
    r = [None] * 30
    r[0], r[1], r[7], r[8], r[9], r[10], r[13], r[14] = pol, name, rt, 1e8, ngg / 166.6667, ngg, sample, field
    r[15], r[19], r[20], r[22], r[25] = inp or name, "2-2-2", ik, smiles, smiles
    return r


RAW = {CD.WW: [[None] * 30, ["hdr"] * 30,
               _ww("+VE", "Chemx", 6.15, 1e7, 1000.0, "CCO", "XXXXXXXXXXXXXX-UHFFFAOYSA-N"),
               _ww("+VE", "Chemx", 6.14, 5e6, 500.0, "CCO", "XXXXXXXXXXXXXX-UHFFFAOYSA-N"),       # duplicate row: max wins
               _ww("+VE", "Estra", 8.29, 1e7, 1000.0, "CCC", "EEEEEEEEEEEEEE-UHFFFAOYSA-N"),
               _ww("-ve", "Wateronly", 3.0, 1e7, 700.0, "CCCC", "WWWWWWWWWWWWWW-UHFFFAOYSA-N"),
               _ww("+VE", "beta-Estra", 8.30, 1e7, 900.0, "CCC", "EEEEEEEEEEEEEE-UHFFFAOYSA-N")],
       CD.PL: [["hdr"] * 30,
               _pl("+VE", "Chemx", 6.13, 20.0, "FIELD 1", "1-3", smiles="CCO", ik="XXXXXXXXXXXXXX-UHFFFAOYSA-N"),
               _pl("+VE", "Chemx", 6.13, 40.0, "4", "10-12", inp="Something else", smiles="CCN", ik="ZZZZZZZZZZZZZZ-UHFFFAOYSA-N"),   # misaligned CompTox block: ignored
               _pl("+VE", "17B-Estra", 6.18, 1.5, "FIELD 2", "4-6", smiles="CCC", ik="EEEEEEEEEEEEEE-UHFFFAOYSA-N")],
       "Semi quant modelling_Soil": [[None] * 30, ["hdr"] * 30]}


def test_parse_structures_and_candidates(cfg):
    water, grain, st = CD.build_cd_tables(RAW, cfg)
    assert len(water) == 5 and len(grain) == 3 and set(grain.field) == {"1", "4", "2"}
    assert st.set_index("name_l").loc["chemx", "smiles"] == "CCO"                    # taken from the aligned row, not the misaligned one
    C, F, D = CD.cd_candidates(water, grain, st, cfg)
    c = C.set_index("name")
    assert c.loc["Chemx", "in_water"] and c.loc["Chemx", "in_grain"] and not c.loc["Wateronly", "in_grain"] and c.loc["Chemx", "ik1"] == "XXXXXXXXXXXXXX"
    w = F[(F.cand_id == "XXXXXXXXXXXXXX") & (F.matrix == "water")]
    assert set(w.field) == {"1", "2", "4"} and (w.conc == 1000.0).all()               # max over duplicate rows, broadcast to every grain field
    g = F[(F.cand_id == "XXXXXXXXXXXXXX") & (F.matrix == "grain")].set_index("field").conc
    assert g["1"] == 20.0 and g["4"] == 40.0
    e = C[C.cand_id == "EEEEEEEEEEEEEE"].iloc[0]                                      # 'Estra' / 'beta-Estra' (water) and '17B-Estra' (grain) = one structure
    assert len(C[C.cand_id == "EEEEEEEEEEEEEE"]) == 1 and e.rt_water == pytest.approx(8.295) and e.rt_grain == pytest.approx(6.18) and e.in_water and e.in_grain


def test_cd_screen_end_to_end_with_fake_fate(cfg, tmp_path, monkeypatch):
    monkeypatch.setattr(CD, "read_cd_sheets", lambda path, cache=None, max_rows=600: RAW)
    cfg2 = copy.deepcopy(cfg); cfg2["paths"]["cache_dir"] = str(tmp_path)
    cfg2["pec_screen"]["fate"]["cache"] = str(tmp_path / "k.csv")
    backend = lambda mw, logkow, pka, base, kaw, dt50_soil_d=None, biowin4=None: 0.01
    res = PR.run_cd_screen(cfg2, "fake.xlsx", str(tmp_path / "out"), props_path=None, backend=backend)
    S = res["summary"].set_index("name")
    est = res["summary"][res["summary"].cand_id == "EEEEEEEEEEEEEE"].iloc[0]
    assert est["class"].startswith("identity conflict")                               # RT 8.29 vs 6.18 although the names differ between the sheets
    assert S.loc["Chemx", "n_pairs"] == 2 and S.loc["Wateronly", "n_pairs"] == 0
    assert (tmp_path / "out" / "pec_screen_results.xlsx").exists()
