"""Builds the GitHub Release assets used by the notebook (GitHub source):

    python colab/make_data_zip.py <engine data folder> [out folder]

  hrms_nts_data.zip        databases + MetFrag jar: SusDat, PubChemLite, MS/MS library, NORMAN RTI, MassBank-EU, MetFrag
  metfrag_spectra.zip      MetFrag's offline MassBank-format spectra (mb_converted, optional - only for MetFrag's spectral score)

Windows-only tools (bundled Java runtime, ThermoRawFileParser) are left out: Colab installs Java / the converter itself.
GitHub Release assets may be up to 2 GB each. Check the licences in engine/DATA_SOURCES.md before attaching them to a PUBLIC release."""
import os, sys, zipfile, hashlib

src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "engine", "data")
outdir = sys.argv[2] if len(sys.argv) > 2 else "."
need = ["susdat.csv", "PubChemLite_exposomics.csv", "spectral_library.sqlite"]
missing = [n for n in need if not os.path.exists(os.path.join(src, n))]
if missing:
    sys.exit(f"{src}: missing {missing}")


def build(out, include):
    n = 0
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as z:
        for root, dirs, files in os.walk(src):
            dirs[:] = [d for d in dirs if d not in ("jre", "ThermoRawFileParser", "__pycache__")]
            for f in files:
                p = os.path.join(root, f); rel = os.path.relpath(p, src)
                if include(rel):
                    z.write(p, rel); n += 1
    hs = hashlib.sha256()
    with open(out, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            hs.update(b)
    size = os.path.getsize(out)
    print(f"{out}: {n} files, {size / 1e9:.2f} GB, sha256 {hs.hexdigest()}")
    if size > 2 * 1024 ** 3:
        print("WARNING: larger than the 2 GB GitHub Release asset limit")


is_mb = lambda rel: rel.replace("\\", "/").startswith("mb_converted/")
build(os.path.join(outdir, "hrms_nts_data.zip"), lambda rel: not is_mb(rel))
if os.path.isdir(os.path.join(src, "mb_converted")):
    build(os.path.join(outdir, "metfrag_spectra.zip"), is_mb)
