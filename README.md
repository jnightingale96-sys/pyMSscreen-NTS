# pyMSscreen-NTS

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jnightingale96-sys/pyMSscreen-NTS/blob/main/colab/HRMS_NTS_Colab.ipynb)
[![Tests](https://github.com/jnightingale96-sys/pyMSscreen-NTS/actions/workflows/tests.yml/badge.svg)](https://github.com/jnightingale96-sys/pyMSscreen-NTS/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24573714.svg)](https://doi.org/10.5281/zenodo.24573714)
[![Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24573798.svg)](https://doi.org/10.5281/zenodo.24573798)

**Fastest way to try it: open the notebook in Google Colab** (badge above), run box 0, put your mzML files on Google Drive and
work through the boxes. Nothing is installed on your computer. Step-by-step instructions: [`colab/README.md`](colab/README.md).

**pyMSscreen-NTS** is a modified version of [pyMSscreen](https://gitlab.com/uniluxembourg/lcsb/eci/pymsscreen)
(University of Luxembourg, LCSB Environmental Cheminformatics; Artistic License 2.0). It keeps the complete pyMSscreen
targeted-prescreening workflow (inputs → extraction → plotting/QA) and adds the **HRMS NTS engine** (`engine/`, package `hrms_nts`):

| Page | What it does |
|---|---|
| Inputs / Extraction / Plotting | original pyMSscreen: compound list → EIC + MS2 extraction → QA flags → interactive review, PDF, MSP library |
| **NTS levels** | non-target screening of DDA mzML files: Level 5 → 4 → 3 → 2b → 2a → 1 (NORMAN SusDat, PubChemLite, MassBank/GNPS libraries, RT prediction, product-ion evidence), results tables, MS2 mirror plots, Excel export, *Send to pyMSscreen prescreen* |
| **Semi-quant** | concentrations without standards: single reference, logIE + elucidation distance, or own calibration curve; matrix conversion (water / grain / soil / extract); calibrant table editing and area measurement |
| **Product ions** | predicted fragments / neutral losses and MRM transitions for a SMILES; score an experimental MS2 spectrum against the prediction |

### Start on Windows (no Docker)
Double-click **`START_pyMSscreen-NTS.bat`**. The first start creates `.venv` and installs everything (needs Python ≥ 3.10, and
Node.js once to build the interface; 5–10 min). Later starts take seconds and open `http://127.0.0.1:8765` in your browser.
Close the black window to stop. The app listens on this computer only.

Before the first non-target run, point the engine at your databases: copy `engine/config.local.example.yaml` to
`engine/config.local.yaml` and edit the paths (SusDat, PubChemLite, MSP libraries); the *NTS levels* page shows what is found.
Build the spectral library once: `.venv\Scripts\python -m hrms_nts build-library`. Engine details, validation status and
the command line: [`engine/README.md`](engine/README.md), [`engine/VALIDATION.md`](engine/VALIDATION.md).

Tests: `.venv\Scripts\python -m pytest backend/tests` and `cd engine && ..\.venv\Scripts\python -m pytest tests`.

Changes to the original pyMSscreen files are limited to: navigation bar and three new pages (`frontend/src/components/`),
same-origin API base for the native build (`api.js`, `npm run build:windows`), Windows working-directory defaults on the
Inputs page, run-number-tolerant polarity detection from file names (`backend/app/utils.py`), and Flask serving the built
frontend (`web_app.py`, `run_windows.py`). The Docker route below still runs the original pyMSscreen pages; the backend
image does not yet install `engine/`, so the three NTS pages currently need the native Windows start.

---

# pyMSscreen (original documentation)

## Overview

pyMSscreen provides a web-based interface for initial prescreening of mass spectrometry data. It combines a Flask backend with a React frontend, and is designed to run via Docker for reproducibility. The application allows users to upload compound information and mzML files, extract chromatographic and MS/MS data using pyOpenMS, and perform automated prescreening.

## Project Structure

- **Backend**: Flask + Gunicorn server with `pyOpenMS` for mzML parsing and data processing
- **Frontend**: React (Vite) served via Nginx
- **Deployment**: Docker + Docker Compose with bind-mounted working directory

## Requirements

- Docker
- Docker Compose
- Local working directory with input files:
  - `compounds.csv`
  - One or more `.mzML` files

## Setup Instructions

### 1. Clone the repository

```bash
git clone https://gitlab.com/uniluxembourg/lcsb/eci/pymsscreen.git
cd your-repo
```

### 2. Create a `.env` file

After cloning the repository, create a `.env` file in the root of the project. This file defines the path to your local data directory, which will be mounted into the backend container at `/app/uploads`.

```bash
touch .env
```

Then, open `.env` in a text editor and add the following line:

```
UPLOADS_PATH=/absolute/path/to/your/data/folder
```

On Windows with Docker Desktop + WSL2, the path should look like(important to include "\\"):

```
UPLOADS_PATH=C:\\Users\\name\\path-to-input-data
```

Ensure that this path exists and contains:

- A `compounds.csv` file
- One or more `.mzML` files. 
- Only `.mzML` files are supported. Vendor formats (e.g., Thermo `.raw`) must be converted beforehand (e.g., using MSConvert / ProteoWizard).


The `compounds.csv` file must include the following columns:

- **`ID`** (required)  
  Unique identifier for each compound.
- **`Name`** (optional)  
  Human-readable compound name.
- **`SMILES`**  
  SMILES representation of the compound structure.
- **`mz`** 
  Neutral mass of the compound.
- **`precurso_mz`**
  precursor mass where no adduct information is required. 

At least one of **`SMILES`** or **`mz`** or **`precursor_mz`** must be provided for each compound.

### Priority Rules

- If **both `SMILES` and `mz` are provided**, the application **prioritizes `SMILES`** for compound processing.
- If only one of these fields is provided, that value will be used.

### 3. Build the Docker containers

Run the following command to build the backend and frontend images, including system libraries required by `pyOpenMS`:

```bash
docker-compose build --no-cache
```

### 4. Start the application

Start both containers using:

```bash
docker-compose up
```

This will launch:

- Frontend: http://localhost:3000
- Backend: http://localhost:8080

The frontend connects to the backend using the environment variable `VITE_API_BASE`, configured in frontend/.env. By default this is set to `http://localhost:8080`, matching the host port mapped in `docker-compose.yml`. If you need to change the backend URL, copy `frontend/.env.example` to `frontend/.env` and update the value:

```bash
VITE_API_BASE=http://your-backend-url
```
Then rebuild the frontend container:

```bash
docker-compose build frontend
docker-compose up
```

If no `.env` file is present, the frontend falls back to `http://localhost:8080` automatically.

## Example Workflow (Use Case 1 from Manuscript)

1. Navigate to `http://localhost:3000`
2. Download the [mzML files](https://zenodo.org/records/14950291) and input [compound list](https://gitlab.com/uniluxembourg/lcsb/eci/pymsscreen/-/blob/main/example_inputs/AAs.csv?ref_type=heads). Place all the input files into a local directory (e.g., /path/to/ example-data)
3. Set this directory in your .env file (as explained above).
4. In the application UI, the Working Directory: /app/uploads/ is set by default. 
5. Upload **Input Files**:
   - The compound CSV ( e.g. AA CSV)
   - One or more `.mzML` files
6. **Tags** (Manual Input): Each mzML file requires a tag to be entered manually in the interface.
   - Tags are user-defined.
   - They are used for identifying samples, grouping and interpretation, output naming
   - Tags must be non-empty and must not contain spaces, as they are used directly in output filenames and sample matching. 

      **Note:**
       - Tags are not automatically extracted from filenames
       - Users should use consistent and meaningful tags, especially for larger datasets.

    Assign a tag to each mzML file as follows:

    | mzML file | Tag |
    |-----------|-----|
    | ISB401_AA_KO_NEG_DLD3_1.mzML | KO-|
    | ISB401_AA_KO_POS_DLD3_1.mzML | KO+ |
    | ISB401_AA_STD_NEG_METmix_1.mzML | STD- |
    | ISB401_AA_STD_POS_METmix_1.mzML | STD+ |
    | ISB401_AA_WT_NEG_FY4_1.mzML | WT- |
    | ISB401_AA_WT_POS_FY4_1.mzML | WT+ |

    
7. **Adduct Selection:** Selected adduct(s): [e.g., [M+H]+, [M-H]- depending on your example]
  - If you are using SMILES or neutral mass (mz) as input, you must select an adduct from the dropdown menu. The selected adduct is used to calculate the expected m/z for extraction.
  - If you are screening directly using precursor mass, selecting an adduct is not required. In this case, extraction is performed around the provided precursor m/z value.

  **Adduct Selection:**

   | Polarity | Adducts |
   |----------|---------|
   | Positive | [M+H]+ |
   | Negative | [M-H]- |


8. Click **"Save State"**
   - This generates a `state.json` file that captures the current application state, including selected inputs and settings.
9. Click **"Load State"**
   The application reads from the previously generated **`state.json`** file.
    - The user **must select the same compound CSV and `.mzML` files** that were used to create the `state.json`.
    - Once **Load State** is clicked, the application resumes from the point where it was last saved, allowing the workflow to continue without repeating previous steps.
10. Click **"Generate Table"**
11. Proceed to the **Extraction** and **Prescreening** steps as needed
12. **Extraction** Settings
    - Mass tolerance (MS1 Fine error): [e.g., 5 ppm or specify your value]
    -  Retention time window: [e.g., ±0.5 min]

    **Extraction Parameters (default settings):**

    | Parameter | Value |
    |-----------|-------|
    | MS1 coarse error | 0.01 Da |
    | MS1 fine error | 5 ppm |
    | MS1 EIC window | 0.01 Da |
    | Retention time window | 0.5 min |

13. **Prescreening / QC Criteria**
    - Signal-to-noise threshold: 3
    - MS1 intensity threshold: 100000
    - Retention time delay: 0.5 min

     **Prescreening / QC Parameters (default settings):**

    | Parameter | Value |
    |-----------|-------|
    | MS1 intensity threshold | 100000 |
    | MS2 intensity threshold | 2500 |
    | Signal-to-noise ratio | 3 |
    | Retention time delay | 0.5 min |

      **Notes**
       - These parameters correspond to the example dataset provided via Zenodo.
       - Users should adjust these values depending on instrument type, acquisition method, and data quality.
       - The parameters are set through the GUI during runtime and stored internally in `extract_config.yaml`.

     Due to environment-specific paths (e.g., Docker-mounted working directories), these configuration files are not directly portable and are therefore not included as static examples in the repository.

14. Finally the *manual QA steps and plots are visualized in the **pyMSscreen Visualization Dashboard**

All outputs (e.g., EIC/fragment CSVs, summary tables) will be saved in the same mounted working directory.

## Troubleshooting

- `FileNotFoundError`: Ensure the path in `.env` exists and contains all required files
- `500 Internal Server Error`: Check backend logs via `docker logs backend`
- `ImportError: libgthread-2.0.so.0`: Make sure your backend Dockerfile installs `libglib2.0-0`

## Platform Compatibility Notes

This application is containerized and intended to run across Linux, macOS, and Windows systems.

- **Linux (native)**: Fully supported
- **macOS**: Fully supported with Docker Desktop
- **Windows**: Requires Docker Desktop with WSL2 backend enabled. Use Linux-style paths in `.env` (e.g., `/mnt/c/Users/...`).

Ensure Docker Desktop has access to the drive via:
- Docker Desktop → Settings → Resources → File Sharing

## Notes

- This application is currently optimized for Linux environments (due to `pyOpenMS` requirements).
- Nginx is configured to support client-side routing (e.g., `/extraction`, `/plotting`).
- All data is accessed through the bind-mounted directory specified via `.env

### How to cite
Nightingale, J. (2026). *pyMSscreen-NTS: HRMS non-target screening, identification confidence and semi-quantification*
(version 0.6.0) [software]. Zenodo. https://doi.org/10.5281/zenodo.24573715

Cite the version you used (each release has its own DOI; https://doi.org/10.5281/zenodo.24573714 always points to the newest).
The reference databases are a separate record (CC BY-NC-SA 4.0; the Colab notebook downloads it automatically):
Nightingale, J. (2026). *pyMSscreen-NTS reference databases (data package v2)* [data set]. Zenodo.
https://doi.org/10.5281/zenodo.24573798

Please also cite pyMSscreen and the databases you rely on ([`engine/DATA_SOURCES.md`](engine/DATA_SOURCES.md)).
