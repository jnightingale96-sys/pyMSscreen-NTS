"""Zip the rebuilt spectral library for sharing: dist/spectral_library_<date>.zip = spectral_library.sqlite + manifest (JSON + README with sources, counts, checks, usage).
The lipid atlas library (if built) goes in its own zip so the main package stays small.    python tools/package_library.py"""
import hashlib, json, os, sys, time, zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); WORK = os.path.join(ROOT, "work"); DIST = os.path.join(ROOT, "dist")
os.makedirs(DIST, exist_ok=True)
M = json.load(open(os.path.join(WORK, "spectral_library_manifest.json"), encoding="utf-8"))
day = time.strftime("%Y%m%d")


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()


def readme(kind, info, srcs, prev=None, sanity=None):
    L = [f"# HRMS NTS toolkit - {kind} spectral library ({time.strftime('%Y-%m-%d')})", "",
         f"* Spectra: **{info['n_spectra']:,}**  ({', '.join(f'{k}: {v:,}' for k, v in info['polarity'].items())})", f"* Unique structures (InChIKey block 1): **{info['unique_structures']:,}**"]
    if prev: L += [f"* Previous library: {prev['n_spectra']:,} spectra, {prev['unique_structures']:,} structures (kept in the toolkit as `work/spectral_library_prev.sqlite`)"]
    L += ["", "## Use", "", "1. Unzip `spectral_library.sqlite` anywhere.", "2. In `config.local.yaml` set `paths: library_db: \"<path>/spectral_library.sqlite\"` (the NTS toolkit / pyMSscreen-NTS reads it directly, SQLite table `spectra`).",
          "3. Level 2a/1 matching then uses it (cosine >= 0.70 and >= 2 matched fragments).  Identifications already cached in a previous run must be re-run (`--no-reuse`) to use the new library.", "",
          "## Contents (per source, after removing duplicates)", "", "| Source | Spectra | Structures |", "|---|---|---|"]
    L += [f"| {s['source']} | {s['spectra']:,} | {s['structures']:,} |" for s in info["by_source"]]
    L += ["", "## How it was built", "", "`tools/rebuild_library.py`: MSP / MetFrag-format files -> SQLite; MS/MS only, ESI/LC-type instruments (GC-EI, MALDI, DART spectra dropped), >= 3 peaks, top 100 peaks, intensities scaled to 1000.",
          "Cross-source duplicates (same structure, polarity, precursor to 0.01 and the same 8 most intense fragments) are removed - the first source in the list wins, so curated libraries are kept in preference to the large GNPS aggregate.", "",
          "| Source file | kept | skipped | duplicates | MD5 |", "|---|---|---|---|---|"]
    L += [f"| `{os.path.basename(s['file'])}` ({s['source']}) | {s['kept']:,} | {s['skipped']:,} | {s['duplicates']:,} | {s.get('md5', '')} |" for s in srcs]
    if sanity:
        L += ["", "## Self-check", "", "Each compound's own library spectrum, perturbed (m/z noise 3 mDa, intensities +-30 %), is searched against the library:", "", "| Compound | Source of query | Top-1 correct | Cosine |", "|---|---|---|---|"]
        L += [f"| {c['compound']} | {c.get('source', '-')} | {c.get('top1_ok', 'not in library')} | {c.get('top1_score', '')} |" for c in sanity]
    L += ["", "Sources / licences: MassBank (CC BY 4.0 / per record), MS-DIAL public libraries (CC BY 4.0, Zenodo 10951913 & 13261737), MoNA (per record), GNPS (CC0 / per library).  Cite the original libraries."]
    return "\n".join(L)


def pack(db, name, info, srcs, prev=None, sanity=None):
    out = os.path.join(DIST, f"{name}_{day}.zip")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        z.write(db, "spectral_library.sqlite")
        z.writestr("README.md", readme(name.replace("_", " "), info, srcs, prev, sanity))
        z.writestr("manifest.json", json.dumps(dict(M, sqlite_md5=md5(db)), indent=1, default=str))
    print(f"{out}: {os.path.getsize(out) / 1e6:.0f} MB  (database {os.path.getsize(db) / 1e6:.0f} MB, md5 {md5(db)})")
    return out


pack(os.path.join(WORK, "spectral_library.sqlite"), "spectral_library", M["main"], M["sources"], M.get("previous"), M.get("sanity_searches"))
ldb = os.path.join(WORK, "spectral_library_lipids.sqlite")
if "lipids" in M and os.path.exists(ldb):
    pack(ldb, "spectral_library_lipids_optional", M["lipids"], M["lipids"]["sources"])
