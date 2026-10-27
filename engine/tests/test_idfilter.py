import copy
import numpy as np, pandas as pd, pytest
from hrms_nts.massbank_eu import split_names, load_massbank_eu, fragment_support
from hrms_nts import idfilter as IDF
from hrms_nts.ms2 import peaks_to_string


def test_split_names_drops_code_tokens_and_blanks():
    assert split_names("LID_235.1805_10.1|Lidocaine| 2-(diethylamino)-N-(2,6-dimethylphenyl)acetamide|nan|Lidocaine") == ["Lidocaine", "2-(diethylamino)-N-(2,6-dimethylphenyl)acetamide"]
    assert split_names("VEN_278.2114_12.5|Venlafaxine") == ["Venlafaxine"] and split_names(float("nan")) == []


def test_load_merges_rows_per_structure(tmp_path):
    cols = ["NORMAN_ID", "NAME", "SMILES", "NEUTRAL_EXACT_MASS", "FORMULA", "CAS_RN", "InChIKey", "DTXSID", "MB_Records Positive Ionization", "Simplified_fragments_Positive",
            "Institute Positive Ionization", "MB_Records Negative Ionization", "Simplified_fragments_Negative", "Institute Negative Ionization"]
    rows = [["a", "LID_235.18_10.1|Lidocaine|Lignocaine", "C", "234.17", "C14H22N2O", "CAS_RN: 137-58-6", "NNJVILVZKWQKPM-UHFFFAOYSA-N", "DTX1", "R1", "86.0964;58.0651", "UFZ", None, None, None],
            ["b", "Lidocaine", "C", "234.17", "C14H22N2O", "CAS_RN: 137-58-6", "NNJVILVZKWQKPM-UHFFFAOYSA-N", "DTX1", "R2", "86.0964;164.1070", "Eawag", None, None, None]]
    p = tmp_path / "mb.csv"; pd.DataFrame(rows, columns=cols).to_csv(p, index=False)
    m = load_massbank_eu(str(p)); r = m.loc["NNJVILVZKWQKPM"]
    assert len(m) == 1 and r["name"] == "Lidocaine" and r.n_rows == 2 and list(r.pos_frags) == [58.0651, 86.0964, 164.107] and r.has_pos and not r.has_neg and r.cas == "137-58-6"


def test_fragment_support_excludes_precursor():
    n, nref = fragment_support([86.0964, 58.0651, 120.0], [235.18, 86.0964, 58.0651, 164.107], precursor_mz=235.1805)
    assert (n, nref) == (2, 3)


def _cands(n_anchor=10):
    rows = []
    for i in range(n_anchor):                                        # anchors on RT = 1 + 0.01 x RTI
        rows.append(dict(cand_id=f"a{i}", name=f"anchor{i}", ik1=f"A{i}", level="2a", lib_score=0.9, polarity="positive", rt=1 + 0.01 * (300 + 40 * i) + (0.1 if i % 2 else -0.1), mz=300.0 + i, ms2_peaks=""))
    rows += [dict(cand_id="good", name="good2b", ik1="G1", level="2b", lib_score=0, polarity="positive", rt=1 + 0.01 * 500, mz=400.0, ms2_peaks=""),
             dict(cand_id="bad", name="bad2b", ik1="B1", level="2b", lib_score=0, polarity="positive", rt=1 + 0.01 * 500 + 6.0, mz=410.0, ms2_peaks=""),
             dict(cand_id="bad_sup", name="bad2b_supported", ik1="S1", level="2b", lib_score=0, polarity="positive", rt=1 + 0.01 * 500 + 6.0, mz=420.0,
                  ms2_peaks=peaks_to_string([100.0, 120.0, 140.0, 160.0], [1e4] * 4)),
             dict(cand_id="nortI", name="nortI2b", ik1="N1", level="2b", lib_score=0, polarity="positive", rt=5.0, mz=430.0, ms2_peaks="")]
    return pd.DataFrame(rows)


def test_id_filter_keeps_2a_and_plausible_2b_and_mbeu_rescue(cfg):
    c = _cands(); ids = list(c.ik1)
    rti = pd.DataFrame({"ik1": ids, "pred_rti_pos": [300 + 40 * i for i in range(10)] + [500, 500, 500, np.nan], "pred_rti_neg": np.nan})
    mb = pd.DataFrame({"ik1": ["S1"], "pos_frags": [np.array([100.0, 120.0, 140.0, 160.0, 400.0])], "neg_frags": [np.array([])]}).set_index("ik1", drop=False)
    cfg2 = copy.deepcopy(cfg); cfg2["pec_screen"]["id_filter"]["massbank_eu"] = None
    kept, removed, models = IDF.apply_id_filter(c, cfg2, None, rti=rti, mbeu=mb)
    st = kept.set_index("cand_id").id_filter
    assert all(f"a{i}" in st.index for i in range(10))                       # every 2a stays
    assert st["good"].startswith("2b + RT plausible")
    assert st["bad_sup"].startswith("2b + MassBank-EU")                      # RT-implausible but supported by measured fragments
    assert list(removed.cand_id.sort_values()) == ["bad", "nortI"]           # RT implausible / no RTI and no support
    assert models["positive"]["n_anchors"] == 10 and models["positive"]["rmse"] < 0.3 and models["positive"]["tol"] >= 1.0


def test_id_filter_cannot_veto_without_enough_anchors(cfg):
    c = _cands(3); rti = pd.DataFrame({"ik1": list(c.ik1), "pred_rti_pos": 400.0, "pred_rti_neg": np.nan})
    cfg2 = copy.deepcopy(cfg); cfg2["pec_screen"]["id_filter"]["massbank_eu"] = None
    kept, removed, models = IDF.apply_id_filter(c, cfg2, None, rti=rti, mbeu=None)
    assert models == {} and len(removed) == 0 and kept[kept.level == "2b"].id_filter.str.contains("cannot be judged").all()


def test_applicability_domain_limits_rt_veto(cfg):
    c = _cands(); ids = list(c.ik1)
    rti = pd.DataFrame({"ik1": ids, "pred_rti_pos": [300 + 40 * i for i in range(10)] + [500, 500, 500, 500], "pred_rti_neg": np.nan,
                        "ad_pos": [True] * 10 + [True, False, True, True], "ad_neg": True})        # 'bad' is outside the RTI domain -> RT cannot judge it
    cfg2 = copy.deepcopy(cfg); cfg2["pec_screen"]["id_filter"]["massbank_eu"] = None
    kept, removed, models = IDF.apply_id_filter(c, cfg2, None, rti=rti, mbeu=None)
    r = removed.set_index("cand_id")
    assert "outside AD" in r.at["bad", "id_filter"] and "IMPLAUSIBLE" not in r.at["bad", "id_filter"] and models["positive"]["n_anchors"] == 10
