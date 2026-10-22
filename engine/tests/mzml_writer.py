"""Minimal mzML writer + synthetic LC-MS dataset used by the tests and by `examples/make_demo_mzml.py`
(a tiny, fully-known dataset so anyone can test the whole pipeline without instrument data)."""
from __future__ import annotations
import base64, os
import numpy as np

CBZ = 236.094963 + 1.007276          # carbamazepine [M+H]+


def _b64(a):
    return base64.b64encode(np.asarray(a, dtype="<f8").tobytes()).decode()


def _arrays(mz, it):
    return ('<binaryDataArrayList count="2">'
            f'<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>'
            f'<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>'
            f'<binary>{_b64(mz)}</binary></binaryDataArray>'
            f'<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>'
            f'<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>'
            f'<binary>{_b64(it)}</binary></binaryDataArray></binaryDataArrayList>')


def write_mzml(path, scans, polarity="positive"):
    """scans: list of dict(level, rt_min, mz, int, [precursor_mz]).  Writes a plain (non-indexed) centroided mzML."""
    pol = ('<cvParam cvRef="MS" accession="MS:1000130" name="positive scan"/>' if polarity == "positive"
           else '<cvParam cvRef="MS" accession="MS:1000129" name="negative scan"/>')
    out = ['<?xml version="1.0" encoding="utf-8"?>',
           '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0"><cvList count="1"><cv id="MS" fullName="PSI-MS" URI="x" version="4"/></cvList>',
           f'<run id="synthetic"><spectrumList count="{len(scans)}">']
    for i, s in enumerate(scans):
        mz = np.asarray(s["mz"], float); it = np.asarray(s["int"], float)
        o = np.argsort(mz); mz, it = mz[o], it[o]
        pre = ""
        if s["level"] == 2:
            pre = (f'<precursorList count="1"><precursor><isolationWindow><cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{s["precursor_mz"]}"/>'
                   f'<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="0.5"/><cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.5"/></isolationWindow>'
                   f'<selectedIonList count="1"><selectedIon><cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{s["precursor_mz"]}"/></selectedIon></selectedIonList>'
                   f'<activation><cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation"/><cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="30"/></activation></precursor></precursorList>')
        out.append(f'<spectrum index="{i}" id="scan={i + 1}" defaultArrayLength="{len(mz)}">'
                   f'<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{s["level"]}"/>'
                   f'<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>{pol}'
                   f'<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="{float(it.sum())}"/>'
                   f'<scanList count="1"><scan><cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{s["rt_min"]}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/></scan></scanList>'
                   f'{pre}{_arrays(mz, it)}</spectrum>')
    out.append('</spectrumList></run></mzML>')
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))
    return path


def synthetic_scans(n=400, seed=0, with_ms2=True):
    """Carbamazepine peak (8.0 min, 5e7 counts, M+1 at 16.6 %), a continuous background ion with two bumps, 30 noise peaks/scan."""
    rng = np.random.default_rng(seed)
    scans = []
    for t in np.linspace(0, 20, n):
        peaks = []
        g = 5e7 * np.exp(-0.5 * ((t - 8.0) / 0.08) ** 2)
        if g > 2e4:
            peaks += [(CBZ + rng.normal(0, 2e-4), g), (CBZ + 1.00335, g * 0.166)]
        b = 3e5 * (1 + 3 * np.exp(-0.5 * ((t - 5) / 0.1) ** 2) + 3 * np.exp(-0.5 * ((t - 12) / 0.1) ** 2))
        peaks.append((150.0 + 0.0003 * rng.normal(), b))
        peaks += [(rng.uniform(100, 600), rng.uniform(2e3, 9e3)) for _ in range(30)]
        a = np.array(sorted(peaks))
        scans.append(dict(level=1, rt_min=float(t), mz=a[:, 0], int=a[:, 1]))
        if with_ms2 and abs(t - 8.0) < 0.12 and len(scans) % 2 == 0:
            scans.append(dict(level=2, rt_min=float(t) + 0.001, precursor_mz=CBZ,
                              mz=[167.0730, 180.0808, 192.0808, 194.0964, 237.1022], int=[1e4, 2e4, 3e4, 1e5, 5e3]))
    return scans
