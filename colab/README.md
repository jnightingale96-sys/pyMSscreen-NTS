# HRMS NTS in Google Colab

`HRMS_NTS_Colab.ipynb` runs the whole workflow in Colab, one box per step:
setup → folders → .raw conversion → sample sheet → QA/QC settings → calibration → IS / targets → non-target screening →
QA (replicates, QC, blanks, peak rating, presence confidence) → semi-quantification (elucidation distance + CompTox) →
chromatogram viewer (Q1 / Q3) → report.

## First time - from GitHub (recommended)
1. The code lives in the GitHub repository; the databases (`hrms_nts_data.zip`, 0.33 GB: NORMAN SusDat, PubChemLite,
   970k-spectrum MS/MS library, NORMAN RTI table, MassBank-EU) are attached to a Release tagged `data-v1`
   (build it with `python colab/make_data_zip.py <engine/data> hrms_nts_data.zip`).
2. Open the notebook straight from GitHub: in Colab, *File → Open notebook → GitHub*, tick *Include private repos* if the repository
   is private, choose `colab/HRMS_NTS_Colab.ipynb`. For a public repository the link is
   `https://colab.research.google.com/github/<user>/<repo>/blob/main/colab/HRMS_NTS_Colab.ipynb`.
3. Private repository: create a fine-grained GitHub token with read-only *Contents* access to that repository and add it as a Colab
   secret named `GH_TOKEN` (key icon in the left bar). Box 0 then clones the code and downloads the Release zip with it.
4. Put your data in a Google Drive folder (`.raw` and/or `.mzML`; one folder per batch). SCIEX `.wiff/.wiff2` must be converted
   on Windows with ProteoWizard msconvert first (vendor DLLs are Windows-only); profile mzML is fine.
5. Set the repository / folder names in boxes 0 and 1 and run the boxes in order.
Licences: the databases have their own licences (`engine/DATA_SOURCES.md`; some library records are non-commercial) - keep the
repository / Release **private** unless you have checked that redistribution is allowed.

## Alternative - Google Drive zip
Put the full package zip (`HRMS_NTS_Ultimate_*.zip`, 1.2 GB) in your Google Drive and choose `SOURCE = Google Drive zip` in box 0.
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
