"""Full MS/MS spectral-library rebuild from every installed open library (MassBank, MS-DIAL public libraries, MoNA, GNPS):
cross-source duplicate removal, InChIKeys from SMILES, manifest with counts/MD5s and self-check searches.

Source files are looked up under `paths.library_sources_dir` in the layout massbank/ msdial/ mona/ gnps/ lipidblast/; missing files are
skipped and listed in the manifest.  Source order matters for de-duplication (first copy of a spectrum wins): curated/experimental
libraries first, the large aggregate (GNPS) last."""
from __future__ import annotations
import glob, hashlib, json, logging, os, sqlite3, time
import numpy as np

from .library import build_library, SpectralLibrary, fill_inchikeys_from_smiles

log = logging.getLogger("hrms_nts")

MSDIAL_ORDER = ["MSMS-Pos-MassBankEU", "MSMS-Neg-MassBankEU", "MSMS-Pos-Respect", "MSMS-Neg-Respect", "MSMS-Pos-PlaSMA", "MSMS-Neg-PlaSMA",
                "BioMSMS-Pos-PlaSMA", "BioMSMS-Neg-PlaSMA", "MSMS-Neg-PFAS_20200806", "MSMS-Pos-Vaniya-Fiehn_Natural_Products_Library_20200109",
                "MSMS-Neg-Vaniya-Fiehn_Natural_Products_Library_20200109", "MSMS-Pos-FiehnHILIC", "MSMS-Neg-FiehnHILIC", "MSMS-Pos-CASMI2016",
                "MSMS-Neg-CASMI2016", "MSMS-Pos-MetaboBASE", "MSMS-Neg-MetaboBASE", "MSMS-Pos-Pathogen_Box_20200109", "KI-GIAR_zic-HILIC_Pos_v0.90",
                "MSMS-Neg-RikenOxPLs", "MSMS-Pos-bmdms-np_20200811", "MSMS-Public_experimentspectra-neg-VS19", "MSMS-Public_experimentspectra-pos-VS19"]
SANITY_COMPOUNDS = ("carbamazepine", "venlafaxine", "lidocaine", "diazinon", "tramadol", "caffeine", "atenolol", "perfluorooctanoic acid",
                    "perfluorooctanesulfonic acid")


def sources(db_dir, inhouse_msp=None):
    """(name, path) of the source files that exist, and the names of those that do not."""
    S = [("MassBank_2026.03", os.path.join(db_dir, "massbank", "MassBank_NISTformat.msp"))]
    S += [("msdial:" + n, os.path.join(db_dir, "msdial", n + ".msp")) for n in MSDIAL_ORDER]
    S.append(("MoNA_LCMSMS_20241014", os.path.join(db_dir, "mona", "MoNA-export-LC-MS-MS_Spectra-20241014-0.005.mb")))
    S.append(("GNPS_pesticides_pos", os.path.join(db_dir, "gnps", "GNPS-COLLECTIONS-PESTICIDES-POSITIVE.msp")))
    S.append(("GNPS_pesticides_neg", os.path.join(db_dir, "gnps", "GNPS-COLLECTIONS-PESTICIDES-NEGATIVE.msp")))
    S.append(("GNPS_ALL_NO_PROPAGATED", os.path.join(db_dir, "gnps", "ALL_GNPS_NO_PROPOGATED.msp")))
    if inhouse_msp:
        S.insert(0, ("inhouse", inhouse_msp))           # own standards first: they win de-duplication and are tagged 'inhouse' (Level 1)
    return [(n, p) for n, p in S if os.path.exists(p)], [n for n, p in S if not os.path.exists(p)]


def lipid_sources(db_dir):
    d = os.path.join(db_dir, "lipidblast")
    return [(f"msdial_lipid_atlas:{os.path.basename(p)[:-4]}", p) for p in sorted(glob.glob(os.path.join(d, "MSDIAL-TandemMassSpectralAtlas-VS69-*.msp")))]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def summarize(dbp):
    con = sqlite3.connect(dbp)
    n = con.execute("select count(*) from spectra").fetchone()[0]
    by = con.execute("select polarity, count(*) from spectra group by polarity").fetchall()
    u = con.execute("select count(distinct ik1) from spectra where ik1<>''").fetchone()[0]
    src = con.execute("select source, count(*), count(distinct ik1) from spectra group by source order by 2 desc").fetchall()
    con.close()
    return dict(n_spectra=n, polarity=dict(by), unique_structures=u, by_source=[dict(source=a, spectra=b, structures=c) for a, b, c in src])


def sanity(dbp, compounds=SANITY_COMPOUNDS):
    """Search the library with a perturbed copy (m/z noise 3 mDa, intensity +-30 %) of each compound's own spectrum: it must come back top-1."""
    lib = SpectralLibrary(dbp); con = sqlite3.connect(dbp); out = []
    for nm in compounds:
        r = con.execute("select id,name,inchikey,precursor_mz,polarity,source,peaks from spectra where lower(name)=? order by n_peaks desc limit 1", (nm,)).fetchone()
        if not r:
            out.append(dict(compound=nm, found_in_library=False)); continue
        v = np.frombuffer(r[6], dtype=np.float32)
        pk = (v if v.size % 2 == 0 else v[:-1]).reshape(-1, 2)
        rng = np.random.default_rng(0)
        q_mz = pk[:, 0].astype(float) + rng.normal(0, 0.003, len(pk)); q_int = pk[:, 1].astype(float) * rng.uniform(0.7, 1.3, len(pk))
        o = np.argsort(q_mz); res = lib.search(float(r[3]), r[4], q_mz[o], q_int[o], 20.0, 0.01, 20.0, 0.3, True)
        top = res.iloc[0] if len(res) else None
        exp = r[2].split("-")[0]
        out.append(dict(compound=nm, found_in_library=True, source=r[5], top1_ik=(top["ik1"] if top is not None else None), expected_ik=exp,
                        top1_ok=bool(top is not None and top["ik1"] == exp), top1_score=round(float(top["score"]), 3) if top is not None else None))
    con.close(); lib.con.close()
    return out


def rebuild(main_db, db_dir, inhouse_msp=None, lipids=True, limit=None):
    """Build <main_db> from all sources under db_dir.  The previous library is kept as *_prev.sqlite; a manifest JSON is written next to it."""
    work = os.path.dirname(os.path.abspath(main_db))
    new_db = main_db.replace(".sqlite", "_new.sqlite"); prev = main_db.replace(".sqlite", "_prev.sqlite")
    prev_info = summarize(main_db) if os.path.exists(main_db) else None
    srcs, missing = sources(db_dir, inhouse_msp)
    if not srcs:
        raise FileNotFoundError(f"no library source files found under {db_dir} (expected massbank/ msdial/ mona/ gnps/)")
    if missing:
        log.warning("library sources not found (skipped): %s", missing)
    st = []; t0 = time.time()
    n = build_library([p for _, p in srcs], new_db, [s for s, _ in srcs], dedupe=True, stats=st, limit=limit)
    log.info("main library: %d spectra in %.0f min", n, (time.time() - t0) / 60)
    log.info("InChIKeys derived from SMILES for %d spectra", fill_inchikeys_from_smiles(new_db))
    info = summarize(new_db); chk = sanity(new_db)
    if os.path.exists(main_db):
        if os.path.exists(prev):
            os.replace(prev, prev + ".old")
        os.replace(main_db, prev)
    os.replace(new_db, main_db)
    manifest = dict(built=time.strftime("%Y-%m-%d %H:%M:%S"), main=info, previous=prev_info,
                    sources=[dict(s, md5=md5(s["file"]) if os.path.getsize(s["file"]) < 3e9 else "") for s in st],
                    missing_sources=missing, sanity_searches=chk)
    if lipids and lipid_sources(db_dir):
        ls = lipid_sources(db_dir); st2 = []; t1 = time.time(); ldb = os.path.join(work, "spectral_library_lipids.sqlite")
        n2 = build_library([p for _, p in ls], ldb, [s for s, _ in ls], dedupe=True, stats=st2)
        manifest["lipids"] = dict(summarize(ldb), sources=st2)
        log.info("lipid library: %d spectra in %.0f min", n2, (time.time() - t1) / 60)
    with open(os.path.join(work, "spectral_library_manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=1, default=str)
    return manifest
