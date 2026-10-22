"""Create a tiny synthetic mzML (carbamazepine peak at 8.0 min + background + MS2) so the whole pipeline can be tried without instrument data:

    python examples/make_demo_mzml.py            -> examples/demo_POS_01.mzML
    python -m hrms_nts --set paths.susdat=null --set paths.pubchemlite=null --set paths.user_compounds=examples/compounds_example.csv \\
           run examples/demo_POS_01.mzML --no-library --out work/demo
"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(here), "tests"))
from mzml_writer import write_mzml, synthetic_scans

p = write_mzml(os.path.join(here, "demo_POS_01.mzML"), synthetic_scans())
print("wrote", p)
