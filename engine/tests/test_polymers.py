import pytest
from hrms_nts.polymers import classify_polymer


@pytest.mark.parametrize("name,smiles,expect", [
    ("Hexaethylene glycol", "OCCOCCOCCOCCOCCOCCO", True),
    ("PPG n7", "", True),
    ("3,6,9,12,15-Pentaoxanonadecan-1-ol", "CCCCOCCOCCOCCOCCOCCO", True),
    ("Tetraethylene glycol monoethyl ether", "CCOCCOCCOCCOCCO", True),
    ("Polyoxyethylene (9) lauryl ether", "", True),
    ("Diethylene glycol", "OCCOCCO", False),                                        # solvent, only two units
    ("Carbamazepine", "NC(=O)N1c2ccccc2C=Cc2ccccc12", False),
    ("2-(2-(2-Hydroxyethoxy)ethoxy)ethyl octanoate", "CCCCCCCC(=O)OCCOCCOCCO", False),     # ester: not a pure polyether
    ("Tripropylene glycol diacrylate", "O=C(OCC(OCC(OCC(OC(=O)C=C)C)C)C)C=C", False),      # monomer
    ("15-Crown-5", "C1COCCOCCOCCOCCO1", False),                                      # ring
])
def test_classify_polymer(name, smiles, expect):
    assert classify_polymer(name, smiles)[0] == expect
