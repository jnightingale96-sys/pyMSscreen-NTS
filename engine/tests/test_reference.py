import pandas as pd

from hrms_nts import reference as R
from hrms_nts.database import load_sle


def test_annotate_labels_and_hint(tmp_path):
    idx = pd.DataFrame([
        dict(ik1="FFGPTBGBLSHEP", name="Carbamazepine", norman_lists="S010;S057", n_norman_lists=2, known_tp=False, np_sources="",
             lotus_n_organisms=0, npatlas_origin="", chebi_roles="", lipid_category="", lipid_class=""),
        dict(ik1="AAAAAAAAAAAAAA", name="Sitosterol", norman_lists="", n_norman_lists=0, known_tp=False,
             np_sources="COCONUT;LOTUS;ChEBI", lotus_n_organisms=900, npatlas_origin="", chebi_roles="plant metabolite",
             lipid_category="Sterol Lipids [ST]", lipid_class="Sterols"),
        dict(ik1="BBBBBBBBBBBBBB", name="TP1", norman_lists="S066", n_norman_lists=1, known_tp=True, np_sources="",
             lotus_n_organisms=0, npatlas_origin="", chebi_roles="", lipid_category="", lipid_class=""),
    ])
    idx["ik1"] = idx.ik1.str[:14]
    idx.to_csv(tmp_path / "ri.csv.gz", index=False)
    pd.DataFrame([dict(ik1_tp="BBBBBBBBBBBBBB", tp_name="TP1", parent_name="Parent X", parent_ik1="", transformation="hydroxylation",
                       source="S066")]).to_csv(tmp_path / "tp.csv.gz", index=False)
    cfg = {"paths": {"reference_index": str(tmp_path / "ri.csv.gz"), "tp_pairs": str(tmp_path / "tp.csv.gz")}, "_root": str(tmp_path)}
    R._CACHE.clear()
    T = pd.DataFrame(dict(name=["CBZ", "sterol", "tp", "unknown"],
                          ik1=["FFGPTBGBLSHEP"[:14], "AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB", "ZZZZZZZZZZZZZZ"]))
    out = R.annotate(T, cfg).set_index("name")
    assert out.loc["CBZ", "origin_hint"] == "anthropogenic (NORMAN lists)"
    assert out.loc["sterol", "origin_hint"] == "natural / endogenous"
    assert bool(out.loc["tp", "known_tp"]) and out.loc["tp", "origin_hint"] == "anthropogenic (NORMAN lists)"
    assert not bool(out.loc["CBZ", "known_tp"])
    assert "Parent X (hydroxylation) [S066]" in out.loc["tp", "tp_of"]
    assert out.loc["unknown", "origin_hint"] == "" and out.loc["unknown", "n_norman_lists"] == 0


def test_annotate_without_index_is_a_no_op():
    R._CACHE.clear()
    T = pd.DataFrame(dict(name=["x"], ik1=["AAAAAAAAAAAAAA"]))
    assert R.annotate(T, {"paths": {}, "_root": "."}).equals(T)


def test_load_sle(tmp_path):
    p = tmp_path / "sle.csv.gz"
    pd.DataFrame([dict(name="Carbamazepine", formula="C15H12N2O", exact_mass=236.09496, smiles="NC(=O)N1c2ccccc2C=Cc2ccccc12",
                       inchikey="FFGPTBGBLSHEPO-UHFFFAOYSA-N", lists="S010;S057")]).to_csv(p, index=False)
    df = load_sle(str(p))
    assert df.source.iloc[0] == "NORMAN-SLE" and df.ik1.iloc[0] == "FFGPTBGBLSHEPO" and df.db_id.iloc[0] == "SLE:S010"
