"""Full library rebuild from a database folder (massbank/ msdial/ mona/ gnps/): de-duplication across sources, previous library kept, manifest."""
import json
import os

from hrms_nts.library import build_library
from hrms_nts.library_build import rebuild

REC = ("NAME: {n}\nPRECURSORMZ: {mz}\nPRECURSORTYPE: [M+H]+\nIONMODE: Positive\nINCHIKEY: {ik}-UHFFFAOYSA-N\nNum Peaks: 3\n"
       "100\t1000\n120\t500\n{p3}\t200\n\n")


def test_rebuild_from_sources_dir(tmp_path):
    src = tmp_path / "databases"
    (src / "massbank").mkdir(parents=True); (src / "gnps").mkdir()
    (src / "massbank" / "MassBank_NISTformat.msp").write_text(
        REC.format(n="carbamazepine", mz=237.1022, ik="FFGPTBGBLSHEPO", p3=194) + REC.format(n="other", mz=300.1, ik="AAAAAAAAAAAAAA", p3=130))
    (src / "gnps" / "ALL_GNPS_NO_PROPOGATED.msp").write_text(
        REC.format(n="carbamazepine", mz=237.1022, ik="FFGPTBGBLSHEPO", p3=194)      # same spectrum as MassBank -> dropped
        + REC.format(n="third", mz=310.2, ik="BBBBBBBBBBBBBB", p3=150))
    db = tmp_path / "work" / "spectral_library.sqlite"
    old = tmp_path / "old.msp"
    old.write_text(REC.format(n="old", mz=200.1, ik="CCCCCCCCCCCCCC", p3=140))
    build_library([str(old)], str(db), ["old"])                   # an older library that must be kept as *_prev

    m = rebuild(str(db), str(src), lipids=False)

    assert m["main"]["n_spectra"] == 3
    assert {s["source"]: s["duplicates"] for s in m["sources"]} == {"MassBank_2026.03": 0, "GNPS_ALL_NO_PROPAGATED": 1}
    assert "MoNA_LCMSMS_20241014" in m["missing_sources"]
    assert os.path.exists(str(db).replace(".sqlite", "_prev.sqlite")) and m["previous"]["n_spectra"] == 1
    sanity = {c["compound"]: c for c in m["sanity_searches"]}
    assert sanity["carbamazepine"]["top1_ok"]
    with open(db.parent / "spectral_library_manifest.json", encoding="utf-8") as f:
        assert json.load(f)["main"]["n_spectra"] == 3
