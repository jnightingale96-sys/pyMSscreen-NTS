"""Rebuild the MS/MS spectral library from every installed open library (same as `python -m hrms_nts build-library` when
paths.library_sources_dir is set) and print a summary.

    python tools/rebuild_library.py [--sources-dir D:/databases] [--no-lipids]

Main library  : work/spectral_library.sqlite          (MassBank, MS-DIAL public libraries, MoNA, GNPS; previous file kept as spectral_library_prev.sqlite)
Lipid library : work/spectral_library_lipids.sqlite    (MS-DIAL Tandem Mass Spectral Atlas VS69 pos/neg, if present)
Then `python tools/package_library.py` zips it for sharing."""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hrms_nts.common import setup_logging, load_config, resolve
from hrms_nts.library_build import rebuild


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources-dir", help="folder with massbank/ msdial/ mona/ gnps/ lipidblast/ (default: paths.library_sources_dir)")
    ap.add_argument("--no-lipids", action="store_true")
    a = ap.parse_args(); setup_logging()
    cfg = load_config()
    src = a.sources_dir or resolve(cfg, cfg["paths"].get("library_sources_dir"))
    if not src or not os.path.isdir(src):
        sys.exit(f"library source folder not found ({src}); set paths.library_sources_dir in config.local.yaml or pass --sources-dir")
    m = rebuild(resolve(cfg, cfg["paths"]["library_db"]), src, resolve(cfg, cfg["paths"].get("inhouse_msp")), lipids=not a.no_lipids)
    info, prev = m["main"], m["previous"] or {}
    print("MAIN:", info["n_spectra"], "spectra,", info["unique_structures"], "structures", info["polarity"],
          "\nPREVIOUS:", prev.get("n_spectra"), "spectra,", prev.get("unique_structures"), "structures")
    for s in info["by_source"][:40]:
        print("  %-62s %8d spectra %7d structures" % (s["source"], s["spectra"], s["structures"]))
    print("SANITY:"); [print("  ", c) for c in m["sanity_searches"]]


if __name__ == "__main__":
    main()
