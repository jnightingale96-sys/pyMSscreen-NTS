import os, numpy as np, pandas as pd, pytest
from hrms_nts import insilico as IS, rt_prediction as RT
from hrms_nts.aggregate import blank_filter, align_features
from hrms_nts.common import load_config


def test_rt_models_recover_descriptor_driven_rt():
    from rdkit import Chem
    from rdkit.Chem import Crippen
    smi = ["CCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "c1ccc2ccccc2c1", "CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
           "NC(=O)N1c2ccccc2C=Cc2ccccc12", "CCCCCCCC", "CCCCCCCCCCO", "Clc1ccccc1", "Brc1ccccc1", "OC(=O)c1ccccc1", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
           "COc1ccc(CCN)cc1", "CCN(CC)CC", "c1ccc(cc1)c1ccccc1", "CC(C)(C)c1ccccc1", "O=C(O)CCCCCCCC(=O)O"]
    rng = np.random.default_rng(1)
    rt = [2 + 1.1 * Crippen.MolLogP(Chem.MolFromSmiles(s)) + rng.normal(0, 0.1) for s in smi]
    rti_db = pd.DataFrame({"ik1": [RT._smiles_ik1(s) for s in smi], "pred_rti_pos": [200 + 60 * r for r in rt], "pred_rti_neg": np.nan})
    m = RT.RTModel("positive", tol=0.5).fit(pd.DataFrame({"name": smi, "smiles": smi, "rt": rt}), rti_db)
    assert m.qsar is not None and m.qsar_rmse < 1.5
    assert m.rti is not None and m.rti_rmse < 0.3
    right = Crippen.MolLogP(Chem.MolFromSmiles("CCCCCCCCO")) * 1.1 + 2.0
    a = m.annotate(pd.DataFrame({"smiles": ["CCCCCCCCO", "CCCCCCCCO"], "rt_apex": [right, 15.0], "pred_rti_pos": [np.nan, np.nan]}))
    assert bool(a.rt_agree.iloc[0]) and bool(a.rt_conflict.iloc[1])


def test_blank_filter_and_alignment():
    f = pd.DataFrame({"feature_id": list("abcd"), "sample": ["s1", "s2", "blank", "s1"], "polarity": ["positive"] * 4,
                      "mz": [200.0, 200.0005, 200.0002, 300.0], "rt_apex": [5.0, 5.05, 5.02, 8.0], "max_int": [1e7, 2e7, 1e6, 5e6], "area": [1e7, 2e7, 5e6, 5e6]})
    bf = blank_filter(f, ["blank"], 3.0)
    assert bool(bf.loc[0, "blank_flag"]) and not bool(bf.loc[1, "blank_flag"])
    al = align_features(f)
    assert al[al.mz < 201].consensus_id.nunique() == 1 and al[al.mz > 201].consensus_id.nunique() == 1


def test_config_layering_and_overrides(tmp_path):
    base = tmp_path / "config.yaml"; base.write_text("paths: {susdat: a.csv, java: java}\nlevel4: {ppm: 5}\n")
    (tmp_path / "config.local.yaml").write_text("paths: {susdat: b.csv}\n")
    c = load_config(str(base), {"level4": {"ppm": 3}})
    assert c["paths"]["susdat"] == "b.csv" and c["paths"]["java"] == "java" and c["level4"]["ppm"] == 3


@pytest.mark.skipif(not IS.metfrag_available(load_config())[0], reason="Java / MetFrag jar not configured")
def test_metfrag_wrapper(tmp_path):
    cfg = load_config(); ok, jar, java = IS.metfrag_available(cfg)
    c = pd.DataFrame({"db_id": ["A"], "name": ["Carbamazepine"], "formula": ["C15H12N2O"], "exact_mass": [236.094963],
                      "smiles": ["NC(=O)N1c2ccccc2C=Cc2ccccc12"], "inchikey": ["FFGPTBGBLSHEPO-UHFFFAOYSA-N"]})
    r = IS.run_metfrag(c, np.array([167.0730, 180.0808, 192.0808, 194.0964]), np.array([1e4, 2e4, 3e4, 1e5]), 237.1022, "positive", "[M+H]+",
                       str(tmp_path), jar, java, tag="cbz")
    assert not r.empty and float(r["mf_expl"].iloc[0]) >= 3


def test_mb_format_and_duplicate_removal(tmp_path):
    from hrms_nts.library import build_library
    import sqlite3
    mb = tmp_path / "mona.mb"
    mb.write_text("# SampleName = Mecoprop\n# InChIKey = WNTGYJSOUMFZEP-UHFFFAOYSA-N\n# IsPositiveIonMode = False\n# IonizedPrecursorMass = 213.03\n# NumPeaks = 3\n"
                  "141 100\n105 45\n71 9\n\n")
    msp = tmp_path / "x.msp"
    msp.write_text("NAME: Mecoprop\nPRECURSORMZ: 213.03\nPRECURSORTYPE: [M-H]-\nIONMODE: Negative\nINCHIKEY: WNTGYJSOUMFZEP-UHFFFAOYSA-N\nNum Peaks: 3\n141\t1000\n105\t450\n71\t90\n\n"
                   "NAME: Other\nPRECURSORMZ: 300.1\nPRECURSORTYPE: [M+H]+\nIONMODE: Positive\nINCHIKEY: AAAAAAAAAAAAAA-UHFFFAOYSA-N\nNum Peaks: 3\n100\t1000\n120\t500\n130\t200\n\n")
    db = tmp_path / "lib.sqlite"; st = []
    n = build_library([str(mb), str(msp)], str(db), ["mona", "msdial"], dedupe=True, stats=st)
    con = sqlite3.connect(db); rows = con.execute("select name, source, polarity, precursor_mz from spectra order by id").fetchall(); con.close()
    assert n == 2 and len(rows) == 2 and rows[0][1] == "mona" and rows[0][2] == "negative"       # the identical MSP copy of Mecoprop is dropped, 'Other' kept
    assert st[1]["duplicates"] == 1


def test_orbitrap_kept_gc_ei_dropped(tmp_path):
    from hrms_nts.library import build_library
    msp = tmp_path / "x.msp"
    peaks = "100\t1000\n120\t500\n130\t200\n\n"
    rec = lambda n, inst: ("NAME: %s\nPRECURSORMZ: 300.1\nPRECURSORTYPE: [M+H]+\nIONMODE: Positive\nINSTRUMENTTYPE: %s\nINCHIKEY: AAAAAAAAAAAAAA-UHFFFAOYSA-N\nNum Peaks: 3\n" % (n, inst)) + peaks
    msp.write_text(rec("a", "Orbitrap") + rec("b", "GC-EI-TOF") + rec("c", "LC-ESI-QTOF") + rec("d", "MALDI-TOF"))
    assert build_library([str(msp)], str(tmp_path / "l.sqlite"), ["t"]) == 2
