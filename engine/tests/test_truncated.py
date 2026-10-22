"""A truncated mzML (interrupted conversion/copy) is read up to the break and flagged, not silently accepted or fatal."""
import os

from mzml_writer import write_mzml, synthetic_scans
from hrms_nts.mzml_io import read_mzml


def test_truncated_mzml_is_flagged(tmp_path):
    full = tmp_path / "demo_POS_full.mzML"
    write_mzml(str(full), synthetic_scans())
    whole = read_mzml(str(full), use_cache=False)
    assert not whole.truncated

    cut = tmp_path / "demo_POS_cut.mzML"
    data = full.read_bytes()
    cut.write_bytes(data[: int(len(data) * 0.6)])
    part = read_mzml(str(cut), use_cache=False)
    assert part.truncated
    assert 0 < part.n_ms1 < whole.n_ms1
    assert part.ms1_rt.max() < whole.ms1_rt.max()
    assert os.path.exists(cut)
