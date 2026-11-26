"""Build the database zips for the GitHub Release (no Java runtime, no ThermoRawFileParser - Colab installs Java itself).

    python colab/make_data_zip.py <engine/data> hrms_nts_data.zip [metfrag_spectra.zip] [--extra <folder>] [--pubchemlite <csv>]

hrms_nts_data.zip   candidate databases, spectral library, RT lookup, MassBank-EU, MetFrag jar, plus (data-v2) the free-database
                    reference index, literature TP pairs, NORMAN SLE candidates and the SMRT global RT model from --extra
                    (output folder of tools/build_reference_index.py and tools/train_global_rt.py)
metfrag_spectra.zip MetFrag offline spectra (mb_converted/), optional
--pubchemlite       a newer PubChemLite csv to ship instead of data/PubChemLite_exposomics.csv (stored under that name)
"""
import os
import sys
import zipfile

need = ["susdat.csv", "PubChemLite_exposomics.csv", "spectral_library.sqlite"]
EXTRA = ["reference_index.csv.gz", "tp_pairs.csv.gz", "sle_candidates.csv.gz", "global_rt_model.pkl", "global_rt_model_cv.txt"]


def build(src, out, include, extra=None, pubchemlite=None):
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as z:
        for root, _, files in os.walk(src):
            for f in files:
                p = os.path.join(root, f)
                rel = os.path.relpath(p, src).replace("\\", "/")
                if include(rel):
                    if pubchemlite and rel == "PubChemLite_exposomics.csv":
                        continue
                    z.write(p, rel)
        if pubchemlite:
            z.write(pubchemlite, "PubChemLite_exposomics.csv")
            print("PubChemLite from", pubchemlite)
        for f in EXTRA if extra else []:
            p = os.path.join(extra, f)
            if os.path.exists(p):
                z.write(p, f)
                print("added", f)
            else:
                print("missing (skipped):", p)
    import hashlib
    hs = hashlib.sha256()
    with open(out, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            hs.update(b)
    size = os.path.getsize(out)
    print(f"{out}: {size / 1e9:.2f} GB, sha256 {hs.hexdigest()}")
    if size > 2 * 1024 ** 3:
        print("WARNING: larger than the 2 GB GitHub Release asset limit")


def main(a):
    src, out = a[0], a[1]
    ms = a[2] if len(a) > 2 and not a[2].startswith("--") else None
    extra = a[a.index("--extra") + 1] if "--extra" in a else None
    pcl = a[a.index("--pubchemlite") + 1] if "--pubchemlite" in a else None
    miss = [n for n in need if not os.path.exists(os.path.join(src, n))]
    if miss:
        sys.exit(f"missing in {src}: {miss}")
    skip_dirs = ("jre/", "ThermoRawFileParser/", "mb_converted/")
    build(src, out, lambda r: not r.startswith(skip_dirs) and not r.endswith((".zip", ".log")), extra, pcl)
    if ms:
        build(src, ms, lambda r: r.startswith("mb_converted/"))


if __name__ == "__main__":
    main(sys.argv[1:])
