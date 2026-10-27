import pytest
from hrms_nts.properties import structure_flags


@pytest.mark.parametrize("name,smiles,expect", [
    ("glycine", "NCC(=O)O", True),
    ("guanidine + carboxylic acid (peramivir-like)", "CCC(CC)[C@H](NC(C)=O)[C@@H]1[C@H](O)[C@@H](C(=O)O)C[C@H]1NC(=N)N", True),
    ("amine N-oxide", "CCCCCCCCCCCCCC[N+](C)(C)[O-]", True),
    ("betaine", "C[N+](C)(C)CC([O-])=O", True),
    ("ciprofloxacin", "OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O", True),
    ("valsartan (acid + tetrazole, no base)", "CCCCC(=O)N(Cc1ccc(cc1)-c1ccccc1-c1nn[nH]n1)C(C(C)C)C(O)=O", False),
    ("tramadol (base only)", "COc1cccc(c1)C1(O)CCCCC1CN(C)C", False),
    ("ibuprofen (acid only)", "CC(C)Cc1ccc(cc1)C(C)C(O)=O", False),
    ("4-aminobenzoic acid (aniline is not basic at pH 7)", "Nc1ccc(cc1)C(O)=O", False),
    ("nitrobenzene", "c1ccccc1[N+](=O)[O-]", False),
])
def test_zwitterion_flag(name, smiles, expect):
    assert structure_flags(smiles)["zwitterion"] == expect, name
