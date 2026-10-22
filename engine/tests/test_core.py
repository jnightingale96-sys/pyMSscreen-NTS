import numpy as np, pandas as pd, pytest
from hrms_nts import features as F, level4 as L4, levels as LV, semiquant as SQ, insilico as IS, library as LIB
from hrms_nts.common import isotope_pattern, parse_formula, formula_mass, normalize_adduct, read_table

CBZ_SMI = "NC(=O)N1c2ccccc2C=Cc2ccccc12"
PK = ([194.0964, 192.0808, 180.0808, 167.073], [1e5, 3e4, 2e4, 1e4])


def test_formula_and_isotopes():
    assert parse_formula("C15H12N2O") == {"C": 15, "H": 12, "N": 2, "O": 1}
    assert formula_mass("C15H12N2O") == pytest.approx(236.09496, abs=1e-4)
    assert isotope_pattern("C15H12N2O")[1] == pytest.approx(0.171, abs=0.01)
    assert isotope_pattern("C6H4Cl2")[2] == pytest.approx(0.64, abs=0.03)
    assert normalize_adduct("M+H") == "[M+H]+" and normalize_adduct("M-H") == "[M-H]-"


def test_feature_detection(fake_run):
    f = F.detect_features(fake_run, 10, 0.3, 3, 1e5)
    cb = f[abs(f.mz - 237.1022) < 0.01]
    assert len(cb) == 1 and abs(cb.rt_apex.iloc[0] - 8.0) < 0.1
    assert len(f[abs(f.mz - 150.0) < 0.01]) >= 2, "continuous background must be valley-split"
    assert bool(f[abs(f.mz - 150.0) < 0.01].background_flag.any()) is True or f[abs(f.mz - 150.0) < 0.01].mz_occupancy.max() > 0.3
    assert not bool(cb.background_flag.iloc[0])


def test_level4_matching_and_isotopes(fake_run, tiny_db):
    f = F.annotate_relations(F.detect_features(fake_run, 10, 0.3, 3, 1e5))
    hits = L4.match_features(f, tiny_db, {"[M+H]+": 1.007276, "[M+Na]+": 22.989218}, 5, fake_run, 1e4, True)
    assert set(hits.name) == {"Carbamazepine"}, "decoy 0.01 Da away must be rejected"
    assert abs(hits.ppm_error.iloc[0]) < 2 and hits.iso_score.iloc[0] > 0.7


def test_insilico_discriminates():
    good = IS.score_candidate(CBZ_SMI, "[M+H]+", *PK, 237.1022)
    bad = IS.score_candidate("CCCCCCCCCCCCCCCC(=O)O", "[M+H]+", *PK, 237.1022)
    assert good["bi_score"] > bad["bi_score"] and good["bi_expl_adj_n"] >= 2


def test_predict_product_ions_contains_known_losses():
    p = IS.predict_product_ions(CBZ_SMI, "[M+H]+", top=200)
    assert p["precursor_mz"].iloc[0] == pytest.approx(237.1022, abs=1e-3)
    assert (abs(p["mz"] - 194.0964) < 0.003).any(), "loss of HNCO (194.0964) is the textbook carbamazepine fragment"


def test_cosine():
    a = np.array([100., 150., 200.]); b = np.array([1000., 500., 100.])
    s, n = LIB.cosine(a, b, a + 0.001, b)
    assert s == pytest.approx(1.0) and n == 3
    assert LIB.cosine(a, b, a + 0.5, b) == (0.0, 0)


def test_semiquant_single_reference(cfg):
    top = pd.DataFrame({"area": [564162031.8, 5641620.318], "max_int": [1, 1], "rt_apex": [5, 6], "smiles": ["CC", "CC"], "name": ["a", "b"]})
    out = SQ.semiquantify(top, cfg, "grain")
    assert out["conc_grain"].iloc[0] == pytest.approx(83.3333, abs=0.01)      # 0.5 ug/mL * 166.667
    assert out["conc_grain"].iloc[1] == pytest.approx(0.8333, abs=0.001)
    assert out["conc_grain_low"].iloc[0] == pytest.approx(out["conc_grain"].iloc[0] / 10)


BASE = dict(ppm_error=1.0, iso_score=np.nan, n_adducts_seen=1, n_candidates_for_feature=1, broad_flag=False, ms2_n_peaks=5, ms2_good_int=True,
            lib_score=0.0, lib_matched=0, rt_available=True, rt_agree=False, rt_conflict=False, std_rt_diff=np.nan, bi_expl_n=0, bi_expl_adj_n=0.0,
            bi_score=0.0, mf_score=np.nan, clear_top=False)


def _levels(cfg, rows):
    h = pd.DataFrame([{**BASE, **kw, "feature_id": fid, "ik1": fid, "inchikey": fid + "-1"} for fid, kw in rows.items()])
    h["ins_score"] = h["bi_score"]
    return LV.assign_levels(h, cfg).set_index("feature_id")["level"].to_dict()


@pytest.mark.parametrize("name,kw,expected", [
    ("lib+RT -> 2a", dict(lib_score=0.9, lib_matched=4, rt_agree=True), {"2a"}),
    ("strong in-silico but RT conflict -> not 2b", dict(bi_expl_adj_n=4.0, bi_score=0.6, clear_top=True, rt_conflict=True), {"3", "4"}),
    ("standard RT only (no MS2) -> 4", dict(ms2_n_peaks=0, ms2_good_int=False, rt_available=False, std_rt_diff=0.05), {"4"}),
    ("standard RT + library MS2 -> 1", dict(lib_score=0.95, lib_matched=6, rt_agree=True, std_rt_diff=0.05), {"1"}),
    ("clear-top in-silico + RT -> 2b", dict(bi_expl_adj_n=4.0, bi_score=0.6, clear_top=True, rt_agree=True), {"2b"}),
    ("2 ions, RT unknown -> 3", dict(bi_expl_adj_n=2.0, bi_score=0.3, rt_available=False), {"3"}),
    ("2b blocked by background flag", dict(bi_expl_adj_n=4.0, bi_score=0.6, clear_top=True, rt_agree=True, background_flag=True), {"4"}),
    ("2b blocked for tiny molecules", dict(bi_expl_adj_n=4.0, bi_score=0.6, clear_top=True, rt_agree=True, exact_mass=90.0), {"4"}),
    ("library match with RT conflict still stays 2a-capable only if not conflicting", dict(lib_score=0.9, lib_matched=4, rt_conflict=True), {"4", "3"}),
])
def test_level_logic(cfg, name, kw, expected):
    assert _levels(cfg, {"f": kw})["f"] in expected, name


def test_read_table_ignores_comment_lines_but_keeps_hash_in_smiles(tmp_path):
    p = tmp_path / "s.csv"
    p.write_text("name,smiles\n# a comment\nacetonitrile,CC#N\n")
    t = read_table(str(p))
    assert len(t) == 1 and t.smiles.iloc[0] == "CC#N"
