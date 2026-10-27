import pandas as pd, pytest
from hrms_nts.endogenous import classify, flag_table


@pytest.mark.parametrize("name,smiles,expect", [
    ("Citric acid", "OC(=O)CC(O)(CC(O)=O)C(O)=O", True), ("Leucylproline", "CC(C)CC(N)C(=O)N1CCCC1C(O)=O", True),
    ("some fatty acid", "CCCCCCCCCCCCCCCC(=O)O", True), ("unnamed amino acid", "NC(Cc1ccccc1)C(=O)O", True),
    ("Carbamazepine", "NC(=O)N1c2ccccc2C=Cc2ccccc12", False), ("Venlafaxine", "COc1ccc(cc1)C(CN(C)C)C1(O)CCCCC1", False),
    ("PEG n7", "OCCOCCOCCOCCOCCOCCOCCO", False), ("Di(2-ethylhexyl) phthalate", "CCCCC(CC)COC(=O)c1ccccc1C(=O)OCC(CC)CCCC", False),
    ("N-Phenyl-1-naphthylamine", "c1ccc(Nc2cccc3ccccc23)cc1", False), ("Oleamide", "CCCCCCCCC=CCCCCCCCC(N)=O", False),      # slip agent, not a free acid
    ("PEG carboxymethylated acid", "CCCCCCCCCCCOCCOCCOCCOCCOCC(=O)O", False),
    ("Trans-Zeatin", "CC(=CCNc1ncnc2[nH]cnc12)CO", False)])
def test_classify(name, smiles, expect):
    # trans-Zeatin is caught by NAME (curated), the SMILES-only rows test the structure rules
    flag = classify(name, smiles)[0]
    assert flag == (expect or name.lower().startswith("trans-zeatin")), (name, classify(name, smiles))


def test_cd_class_text_and_extra(tmp_path):
    assert classify("Whatever", "", "Metabolic intermediate in the Krebs cycle")[0]
    x = tmp_path / "e.csv"; x.write_text("name,ik1,category\nmyplantthing,,own list\n")
    from hrms_nts.endogenous import load_extra
    assert classify("MyPlantThing", "", "", "", load_extra(str(x)))[:2] == (True, "own list")
    df = pd.DataFrame([dict(name="Succinic acid", smiles="OC(=O)CCC(O)=O", ik1="KDYFGRWQOYBRFD"), dict(name="Lidocaine", smiles="CCN(CC)CC(=O)Nc1c(C)cccc1C", ik1="NNJVILVZKWQKPM")])
    f = flag_table(df); assert list(f.endogenous) == [True, False]


@pytest.mark.parametrize("name,smiles,expect", [("Butanedioic acid", "OC(=O)CCC(O)=O", True), ("Phosphoric acid", "OP(O)(O)=O", True), ("D-Glucitol", "OCC(O)C(O)C(O)C(O)CO", True),
                                                  ("some sugar phosphate", "OCC(O)C(O)C(O)C(=O)COP(O)(O)=O", True), ("Nonanedioic acid", "OC(=O)CCCCCCCC(O)=O", True),
                                                  ("4-Nitrophenol", "Oc1ccc(cc1)[N+]([O-])=O", False), ("Diethyl phthalate", "CCOC(=O)c1ccccc1C(=O)OCC", False),
                                                  ("PEG (hexaethylene glycol)", "OCCOCCOCCOCCOCCOCCO", False), ("Lauryl sulfate", "CCCCCCCCCCCCOS(O)(=O)=O", False)])
def test_more_primary_metabolites(name, smiles, expect):
    assert classify(name, smiles)[0] == expect, (name, classify(name, smiles))
