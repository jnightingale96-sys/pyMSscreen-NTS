# Data sources and tools bundled in `data/`

## Licence summary for public release

* **Code** (this repository): Artistic License 2.0, as the pyMSscreen code it modifies.
* **Database package** (`hrms_nts_data.zip`): a compilation of public databases, each under its own licence (table below).
  Most are CC0 or CC BY 4.0. The most restrictive components are **non-commercial**: individual MassBank / MoNA records
  (per-record CC BY-NC / CC BY-NC-SA, inherited by the MS-DIAL compilations) and NP Atlas (CC BY-NC 4.0, used only as an
  origin label). The package is therefore released under **CC BY-NC-SA 4.0** (non-commercial research use with attribution); anyone
  needing commercial use can rebuild it from the original sources with `hrms-nts build-library`
  (`library_sources_dir`), `tools/build_reference_index.py` and `tools/train_global_rt.py`.
* Thermo's RawFileReader libraries (inside ThermoRawFileParser) are **not** in the database package (Colab installs the
  converter itself); the Windows package bundles them under the Thermo licence.

The files below come from public sources and are bundled for convenience. **Licences and terms can change - check each source before passing the package on or publishing,
and cite the sources you rely on.** Nothing here is our own measured data except the in-silico settings.

| Item | Source | Notes |
|---|---|---|
| `susdat.csv` | NORMAN Suspect List Exchange / SusDat (norman-network.com/nds/susdat), export dated 2026-03-29 | candidate chemicals, predicted RT index (RTI), toxicity predictions; cite the NORMAN network and the SusDat Zenodo record |
| `Updated_NORMAN_with_TPSA.xlsx` | NORMAN SusDat table extended with RTI applicability-domain flags and calculated TPSA | used for the RT plausibility check; RTI method: Aalizadeh et al., *Anal. Chem.* 2021 |
| `PubChemLite_exposomics.csv` | PubChemLite for Exposomics; data-v1: release 2024-12-27, data-v2: CCSbase release 2026-09-25 (Zenodo 22953182; 469,602 compounds, predicted CCS for 11 adducts) (pubchemlite.lcsb.uni.lu) | candidate chemicals with annotation counts; CC BY 4.0; cite Schymanski/Bolton et al. PubChemLite |
| `sle_candidates.csv.gz` (data-v2) | NORMAN Suspect List Exchange, all collections S1-S137 (Zenodo community norman-sle, 2026); 447,453 structures from 128 collections with structure columns | extra candidates + list membership; CC BY 4.0; cite Mohammed Taha et al., *Environ. Sci. Eur.* 2022 and the individual lists |
| `reference_index.csv.gz` (data-v2) | built by `tools/build_reference_index.py` from NORMAN SLE, **COCONUT** 2.0 (2026-10, CC0), **LOTUS** (frozen 2023-01-06, Zenodo 7534071, CC BY 4.0), **NP Atlas** v2024_09 (CC BY-NC 4.0, academic use), **LIPID MAPS** LMSD (CC BY 4.0), **ChEBI** 2026-10 with ontology roles (CC BY 4.0); 957,782 structures | origin labels (anthropogenic / natural / lipid / metabolite roles); cite Sorokina et al. 2021 (COCONUT), Rutz et al. 2022 (LOTUS), van Santen et al. 2022 (NP Atlas), Fahy et al. 2009 (LIPID MAPS), Hastings et al. 2016 (ChEBI) |
| `tp_pairs.csv.gz` (data-v2) | literature parent -> transformation-product pairs from NORMAN SLE S66 (Eawag TPs), S68 (HSDB TPs), S73 (METXBIODB), S74 (REFTPS), S113, S121 (Eawag-BBD); 7,898 pairs | "TP of" column; cite the individual lists (e.g. Schymanski et al., PubChem transformations) |
| `global_rt_model.pkl` (data-v2) | gradient-boosting model trained by `tools/train_global_rt.py` on the **METLIN SMRT** dataset (80,038 compounds, RP C18; Domingo-Almenara et al., *Nat. Commun.* 2019, CC BY 4.0) | method-independent retention index, transferred to each method with its standards / internal standards as anchors |
| `MassBankEU_compounds.csv` | NORMAN / MassBank-EU compound list (2017, with measured MS/MS fragment annotations) | extra fragment support for in-silico identifications |
| `spectral_library.sqlite` | built by `hrms-nts build-library` on 2026-10-23 from **GNPS** (ALL_GNPS_NO_PROPAGATED), **MassBank** release 2026.03, **MoNA** LC-MS/MS (2024-10-14) and the **MS-DIAL** public libraries (VS19 experimental spectra, Fiehn natural products, PFAS, RESPECT, BMDMS-NP, PlaSMA, MassBank-EU, HILIC ...). 970,218 spectra, 75,561 structures; duplicates removed. `spectral_library_manifest.json` lists every source file with checksum | respect the individual library licences (mostly CC0 / CC BY); cite GNPS (Wang et al. 2016), MassBank (Horai et al. 2010), MoNA, MS-DIAL (Tsugawa et al. 2015) |
| `MetFragCommandLine-2.6.11.jar`, `mb_converted/` | MetFrag command line (IPB Halle, github.com/ipb-halle/MetFragRelaunched, LGPL) and MassBank / GNPS / MoNA spectra converted to MetFrag's offline format | in-silico fragmentation; cite Ruttkies et al., *J. Cheminform.* 2016 |
| `jre/` | Eclipse Temurin 21 Java runtime (adoptium.net, GPLv2 with Classpath Exception) | only used to run MetFrag |
| `ThermoRawFileParser/` | ThermoRawFileParser (CompOmics, github.com/CompOmics/ThermoRawFileParser, Apache-2.0) incl. Thermo's RawFileReader libraries under the Thermo licence (see `THERMO_LICENSE`) | `.raw` to mzML conversion; Hulstaert et al., *J. Proteome Res.* 2020 |
| Python packages (installed by `install_nts.bat`) | numpy, pandas, scipy, scikit-learn, RDKit, pyteomics, lxml, matplotlib, openpyxl, pyyaml | installed from PyPI under their own licences |
