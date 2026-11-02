# HRMS NTS in Google Colab

`HRMS_NTS_Colab.ipynb` runs the whole workflow in Colab, one box per step:
setup → folders → .raw conversion → sample sheet → QA/QC settings → calibration → IS / targets → non-target screening →
QA (replicates, QC, blanks, peak rating, presence confidence) → semi-quantification (elucidation distance + CompTox) →
chromatogram viewer (Q1 / Q3) → report.

## First time
1. Put the package zip (`HRMS_NTS_Ultimate_*.zip`, 1.2 GB: engine + NORMAN SusDat, PubChemLite, 970k-spectrum MS/MS library,
   NORMAN RTI table) in your Google Drive, e.g. `MyDrive/HRMS_NTS/`.
2. Put your data in a Drive folder (`.raw` and/or `.mzML`; one folder per batch). SCIEX `.wiff/.wiff2` must be converted
   on Windows with ProteoWizard msconvert first (vendor DLLs are Windows-only); profile mzML is fine.
3. Open `HRMS_NTS_Colab.ipynb` in Colab (File → Upload notebook, or open it from Drive), set the paths in boxes 0 and 1, run the boxes in order.
4. Optional, for semi-quantification: request a free CompTox CTX API key from the EPA and add it as a Colab secret named `CTX_API_KEY`
   (key icon in the left bar). Without a key, a Dashboard batch-search export or RDKit logP estimates are used.

Every step saves its results in the results folder (pickles + CSV + Excel), so you can reconnect and continue at any step.

## The same steps without Colab
```python
import sys; sys.path.insert(0, r"...\engine")
from hrms_nts.workflow import Project
P = Project(r"D:\batch", r"D:\batch_results")
P.set_sample_sheet("sample_sheet.csv")          # file,type(sample|blank|qc|cal),group,replicate,conc
P.calibrate("calibration_list.csv")            # name,smiles|formula,label,role(target|IS),is_name
P.targeted("IS_list.csv"); P.screen(); P.qa("2/3"); P.semiquant(api_key="..."); P.report()
P.chromatogram("Carbamazepine", "WW1_rep1")
```

## Tested / not tested
* Tested on Windows with the engine's real data (SCIEX X500R SWATH, 4 injections) and unit tests: sample sheet, calibration
  (linear/quadratic choice), targeted IS screen, NTS screening, QA, semi-quant, report, chromatograms.
* **Not yet tested inside Colab itself** (Drive mount, widget rendering, apt/mono install of ThermoRawFileParser for .raw
  conversion). Run box 0-2 once with a small batch and report problems.
* Regenerate the notebook after editing the cells: `python colab/make_notebook.py`.
