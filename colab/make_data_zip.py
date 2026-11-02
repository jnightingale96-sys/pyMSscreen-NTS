"""Builds hrms_nts_data.zip (the databases only) for a GitHub Release asset used by the Colab notebook (GitHub source).

    python colab/make_data_zip.py <engine data folder> [out.zip]

Leaves out Windows-only tools (bundled Java runtime, ThermoRawFileParser) and the MetFrag jar + offline spectra (MetFrag is off
in Colab; ~2 GB unpacked). GitHub Releases accept
files up to 2 GB; the zip is ~1.1 GB. Check the licences in engine/DATA_SOURCES.md before attaching it to a PUBLIC release."""
import os, sys, zipfile, hashlib

src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "engine", "data")
out = sys.argv[2] if len(sys.argv) > 2 else "hrms_nts_data.zip"
SKIP = ("jre", "ThermoRawFileParser", "mb_converted")          # Windows tools; MetFrag offline spectra (MetFrag is off in Colab)
SKIP_FILES = ("MetFragCommandLine-2.6.11.jar",)
need = ["susdat.csv", "PubChemLite_exposomics.csv", "spectral_library.sqlite"]
missing = [n for n in need if not os.path.exists(os.path.join(src, n))]
if missing:
    sys.exit(f"{src}: missing {missing}")
n = 0
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as z:
    for root, dirs, files in os.walk(src):
        dirs[:] = [d for d in dirs if d not in SKIP and d != "__pycache__"]
        for f in files:
            if f in SKIP_FILES:
                continue
            p = os.path.join(root, f)
            z.write(p, os.path.relpath(p, src)); n += 1
hs = hashlib.sha256()
with open(out, "rb") as fh:
    for b in iter(lambda: fh.read(1 << 20), b""):
        hs.update(b)
h = hs.hexdigest()
size = os.path.getsize(out)
print(f"{out}: {n} files, {size / 1e9:.2f} GB, sha256 {h}")
if size > 2 * 1024 ** 3:
    print("WARNING: larger than the 2 GB GitHub Release asset limit")
