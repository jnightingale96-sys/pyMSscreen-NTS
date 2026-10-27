# Data sources and tools bundled in `data/`

The files below come from public sources and are bundled for convenience. **Licences and terms can change - check each source before passing the package on or publishing,
and cite the sources you rely on.** Nothing here is our own measured data except the in-silico settings.

| Item | Source | Notes |
|---|---|---|
| `susdat.csv` | NORMAN Suspect List Exchange / SusDat (norman-network.com/nds/susdat), export dated 2026-03-29 | candidate chemicals, predicted RT index (RTI), toxicity predictions; cite the NORMAN network and the SusDat Zenodo record |
| `Updated_NORMAN_with_TPSA.xlsx` | NORMAN SusDat table extended with RTI applicability-domain flags and calculated TPSA | used for the RT plausibility check; RTI method: Aalizadeh et al., *Anal. Chem.* 2021 |
| `PubChemLite_exposomics.csv` | PubChemLite for Exposomics, release 2024-12-27 (pubchemlite.lcsb.uni.lu, Zenodo) | candidate chemicals with annotation counts; cite Schymanski/Bolton et al. PubChemLite |
| `MassBankEU_compounds.csv` | NORMAN / MassBank-EU compound list (2017, with measured MS/MS fragment annotations) | extra fragment support for in-silico identifications |
| `spectral_library.sqlite` | built by `hrms-nts build-library` on 2026-10-23 from **GNPS** (ALL_GNPS_NO_PROPAGATED), **MassBank** release 2026.03, **MoNA** LC-MS/MS (2024-10-14) and the **MS-DIAL** public libraries (VS19 experimental spectra, Fiehn natural products, PFAS, RESPECT, BMDMS-NP, PlaSMA, MassBank-EU, HILIC ...). 970,218 spectra, 75,561 structures; duplicates removed. `spectral_library_manifest.json` lists every source file with checksum | respect the individual library licences (mostly CC0 / CC BY); cite GNPS (Wang et al. 2016), MassBank (Horai et al. 2010), MoNA, MS-DIAL (Tsugawa et al. 2015) |
| `MetFragCommandLine-2.6.11.jar`, `mb_converted/` | MetFrag command line (IPB Halle, github.com/ipb-halle/MetFragRelaunched, LGPL) and MassBank / GNPS / MoNA spectra converted to MetFrag's offline format | in-silico fragmentation; cite Ruttkies et al., *J. Cheminform.* 2016 |
| `jre/` | Eclipse Temurin 21 Java runtime (adoptium.net, GPLv2 with Classpath Exception) | only used to run MetFrag |
| `ThermoRawFileParser/` | ThermoRawFileParser (CompOmics, github.com/CompOmics/ThermoRawFileParser, Apache-2.0) incl. Thermo's RawFileReader libraries under the Thermo licence (see `THERMO_LICENSE`) | `.raw` to mzML conversion; Hulstaert et al., *J. Proteome Res.* 2020 |
| Python packages (installed by `install_nts.bat`) | numpy, pandas, scipy, scikit-learn, RDKit, pyteomics, lxml, matplotlib, openpyxl, pyyaml | installed from PyPI under their own licences |
