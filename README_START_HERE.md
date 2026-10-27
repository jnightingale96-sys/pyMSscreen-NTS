# HRMS non-target screening package — start here

One folder with everything needed to screen LC-HRMS data (mzML, DDA): a point-and-click web app (pyMSscreen-NTS) and a
command-line engine (`engine/`, `hrms_nts`), with all databases bundled — NORMAN SusDat, PubChemLite, a 970,218-spectrum MS/MS
library, MetFrag with its own Java runtime, and ThermoRawFileParser for Thermo `.raw` files.

## Use it
1. Unzip to a short path without spaces, e.g. `C:\NTS\HRMS_NTS_Ultimate` (about 5 GB unpacked).
2. Install Python 3.10+ from python.org if needed (tick "Add python.exe to PATH").
3. Double-click **`START_pyMSscreen-NTS.bat`**. The first start installs the Python packages (internet, 5–10 min); then the app opens
   in your browser at http://127.0.0.1:8765. Close the black window to stop it.
   - **NTS levels**: choose a folder of mzML files → Start screening → results by confidence level, structures, MS2 mirror plots, Excel.
   - **Semi-quant**, **Product ions**, and the original pyMSscreen **Inputs / Extraction / Plotting** prescreen pages.
4. Before trusting identifications, put retention times of authentic standards measured on your method in `engine/user_data/standards.csv`
   (≥ 8 per polarity); semi-quantification needs calibrant areas in `engine/user_data/calibrants.csv` (Semi-quant page → Measure areas).

## With Claude Code
Open this folder in Claude Code and ask for the analysis (for example "screen the mzML files in D:\data\run1, blanks are blank1 and blank2,
matrix water"). `CLAUDE.md` tells Claude how to install, check, run, read and report the results.

## More
- Engine method, commands and validation status: `engine/README.md`, `engine/README_START_HERE.md`, `engine/VALIDATION.md`
- Data sources and licences: `engine/DATA_SOURCES.md` (cite them). Web app: modified pyMSscreen (University of Luxembourg, Artistic License 2.0), see `README.md`.
- Research prototype — not validated for regulatory use; semi-quantitative values are order-of-magnitude estimates.
