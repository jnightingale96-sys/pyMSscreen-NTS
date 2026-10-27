"""Build user_data/properties_project.csv (toolkit property-table format) from the project's grain_props.csv (the 'Model data' block used for
the grain PEC), adding the InChIKey block-1 from the SMILES so NTS candidates match by structure, not by name."""
import os, sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SMILES = {"Alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(=NC2)c1ccccc1", "Caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Carbamazepine": "NC(=O)N1c2ccccc2C=Cc2ccccc12",
          "Crotamiton": "CC=CC(=O)N(CC)c1ccccc1C", "Diltiazem": "COc1ccc(cc1)C1Sc2ccccc2N(CCN(C)C)C(=O)C1OC(C)=O", "DIOH-CBZ": "OC1C(O)c2ccccc2N(C(N)=O)c2ccccc12",
          "Gabapentin": "NCC1(CC(O)=O)CCCCC1", "Thiabendazole": "c1ccc2[nH]c(nc2c1)-c1cscn1", "Venlafaxine": "COc1ccc(cc1)C(CN(C)C)C1(O)CCCCC1"}
g = pd.read_csv(os.path.join(ROOT, "work", "grain_props.csv"))
g["smiles"] = g.chemical.map(SMILES)
g["ik1"] = [MolToInchiKey(Chem.MolFromSmiles(s)).split("-")[0] for s in g.smiles]
out = g.rename(columns={"chemical": "name", "solubility_mg_L": "solubility_mg_L"})[["name", "ik1", "smiles", "mw", "logkow", "pka_acid", "pka_base", "type", "kaw", "vp_pa", "solubility_mg_L", "biowin4"]]
out["kaw"] = out["kaw"].where(out["kaw"] > 0, 1e-12)      # as in the earlier grain PEC: missing / zero Kaw -> non-volatile
out["source"] = "project Model-data block (ACD/EPI Suite/OPERA as used in the manuscript model)"
os.makedirs(os.path.join(ROOT, "user_data"), exist_ok=True)
dest = os.path.join(ROOT, "user_data", "properties_project.csv"); out.to_csv(dest, index=False); print(out.round(5).to_string(index=False)); print("->", dest)
