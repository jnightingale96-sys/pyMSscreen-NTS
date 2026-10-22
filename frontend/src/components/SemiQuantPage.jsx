/**
 * @file SemiQuantPage.jsx
 * @summary Semi-quantification of screening results: re-calculate concentrations for an existing results
 * folder (single reference / logIE + elucidation distance / own calibration curve), matrix conversion,
 * calibrant table editing and calibrant-area measurement from a standard mzML.
 */
import { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';
import { ntsInfo, ntsSemiquant, ntsGetCalibrants, ntsSaveCalibrants, ntsStartCalibrate, apiError } from '../api/api';
import {
  useTheme, useSticky, Page, Panel, Btn, Field, TextInput, Select, Row, ErrorBox, LevelChips, DataTable, RunPicker,
  JobMonitor, downloadCsv,
} from './nts/ui';

const MODE_TEXT = {
  single: 'Every feature is quantified against one reference compound: conc = conc_ref × area / area_ref.',
  logie: 'Each feature is matched to the nearest calibrant by elucidation distance (logIE, RT, MW, logP) and corrected by the difference in predicted ionisation efficiency.',
  curve: 'Compounds with ≥ 3 calibration levels use their own regression; all others fall back to logIE / elucidation distance.',
};
const CAL_COLS = ['name', 'smiles', 'rt', 'conc_ug_ml', 'area', 'polarity'];

function CalibrantEditor({ theme, onUseArea }) {
  const [rows, setRows] = useState([]);
  const [path, setPath] = useState('');
  const [err, setErr] = useState('');
  const [msg, setMsg] = useState('');
  const [mzml, setMzml] = useSticky('calMzml', '');
  const [jobId, setJobId] = useState(null);
  const [measured, setMeasured] = useState([]);

  const load = () => ntsGetCalibrants().then((r) => { setRows(r.data.rows); setPath(r.data.path); }).catch((e) => setErr(apiError(e)));
  useEffect(() => { load(); }, []);

  const edit = (i, k, v) => setRows((rs) => rs.map((r, j) => (j === i ? { ...r, [k]: v } : r)));
  const save = async () => {
    setErr(''); setMsg('');
    try {
      const r = await ntsSaveCalibrants(rows);
      setMsg(`Saved ${r.data.n} calibrant row(s) to ${r.data.path}`);
    } catch (e) {
      setErr(apiError(e));
    }
  };
  const measure = async () => {
    setErr('');
    try {
      const r = await ntsStartCalibrate({ mzml });
      setJobId(r.data.job_id);
    } catch (e) {
      setErr(apiError(e));
    }
  };
  const applyMeasured = () => {
    setRows((rs) => rs.map((r) => {
      const m = measured.find((x) => String(x.name).toLowerCase() === String(r.name).toLowerCase());
      return m && m.area ? { ...r, area: m.area, rt: r.rt || m.rt_found } : r;
    }));
    setMsg('Measured areas copied into the table - click Save calibrants to keep them.');
  };

  return (
    <Panel theme={theme} title="Calibrants on this method" right={<span style={{ color: theme.muted, fontSize: '0.8rem' }}>{path}</span>}>
      <p style={{ color: theme.muted, fontSize: '0.88rem', marginTop: 0 }}>
        Needed for <b>logie</b> and <b>curve</b> modes (name, SMILES, RT, concentration in µg/mL extract, peak area). Areas must be measured
        with <i>this</i> software (not Compound Discoverer) — use “Measure areas” below on a calibrant standard mzML.
      </p>
      <div style={{ overflow: 'auto' }}>
        <table style={{ borderCollapse: 'collapse', fontSize: '0.85rem', width: '100%' }}>
          <thead><tr>{[...CAL_COLS, ''].map((c) => <th key={c} style={{ textAlign: 'left', padding: 4, background: theme.tableHeader }}>{c}</th>)}</tr></thead>
          <tbody>
            {rows.map((r, i) => (
              <tr key={i}>
                {CAL_COLS.map((c) => (
                  <td key={c} style={{ padding: 2, minWidth: c === 'smiles' ? 260 : 90 }}>
                    <TextInput theme={theme} value={r[c] ?? ''} onChange={(v) => edit(i, c, v)} style={{ padding: '4px 6px', fontSize: '0.85rem' }} />
                  </td>
                ))}
                <td>
                  <Btn theme={theme} danger onClick={() => setRows((rs) => rs.filter((_, j) => j !== i))} style={{ padding: '3px 8px' }}>✕</Btn>
                  {r.area && <Btn theme={theme} onClick={() => onUseArea(r)} style={{ padding: '3px 8px', marginLeft: 4 }} title="Use as single reference">ref</Btn>}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <Row style={{ marginTop: 10 }}>
        <Btn theme={theme} onClick={() => setRows((rs) => [...rs, { polarity: 'positive' }])}>+ Add row</Btn>
        <Btn theme={theme} primary onClick={save}>Save calibrants</Btn>
        <Btn theme={theme} onClick={load}>Reload</Btn>
      </Row>
      <Row style={{ marginTop: 16 }}>
        <Field theme={theme} label="Measure calibrant areas from a standard mzML (finds each calibrant by m/z and RT)">
          <TextInput theme={theme} value={mzml} onChange={setMzml} placeholder="C:\\data\\Calibrants_POS_1.mzML" />
        </Field>
        <Btn theme={theme} onClick={measure} disabled={!mzml}>Measure areas</Btn>
      </Row>
      <JobMonitor theme={theme} jobId={jobId} onDone={(j) => setMeasured(j.table || [])} />
      {measured.length > 0 && (
        <>
          <DataTable theme={theme} columns={Object.keys(measured[0])} rows={measured} maxHeight={260} />
          <Btn theme={theme} onClick={applyMeasured} style={{ marginTop: 8 }}>Copy measured areas into calibrant table</Btn>
        </>
      )}
      <ErrorBox theme={theme} error={err} />
      {msg && <p style={{ color: theme.ok }}>{msg}</p>}
    </Panel>
  );
}

export default function SemiQuantPage() {
  const theme = useTheme();
  const [info, setInfo] = useState(null);
  const [root, setRoot] = useSticky('resultsRoot', '');
  const [runDir, setRunDir] = useSticky('runDir', '');
  const [mode, setMode] = useSticky('sqMode', 'single');
  const [matrix, setMatrix] = useSticky('sqMatrix', 'water');
  const [levels, setLevels] = useSticky('sqLevels', ['1', '2a', '2b']);
  const [ref, setRef] = useState({ name: '', conc_ug_mL: '', area: '' });
  const [unc, setUnc] = useState('');
  const [res, setRes] = useState(null);
  const [err, setErr] = useState('');
  const [busy, setBusy] = useState(false);

  useEffect(() => {
    ntsInfo().then((r) => {
      setInfo(r.data);
      setRef({ name: r.data.single_reference.name, conc_ug_mL: r.data.single_reference.conc_ug_mL, area: r.data.single_reference.area });
      setUnc(r.data.uncertainty_factor);
    }).catch((e) => setErr(apiError(e)));
  }, []);

  const calc = async () => {
    setErr(''); setBusy(true);
    try {
      const r = await ntsSemiquant({ run_dir: runDir, mode, matrix, levels, uncertainty_factor: unc,
        single_reference: mode === 'single' ? ref : undefined });
      setRes(r.data);
    } catch (e) {
      setErr(apiError(e)); setRes(null);
    } finally {
      setBusy(false);
    }
  };

  const c = `conc_${matrix}`;
  const plotRows = (res?.rows || []).filter((r) => r[c] > 0).slice(0, 40);

  return (
    <Page theme={theme} title="Semi-quantification" subtitle="Concentrations for identified compounds without authentic standards, with × / ÷ uncertainty factor and matrix conversion.">
      <Panel theme={theme} title="1 · Choose screening results">
        <RunPicker theme={theme} root={root} setRoot={setRoot} runDir={runDir} setRunDir={setRunDir} refreshKey={0} />
      </Panel>
      <Panel theme={theme} title="2 · Method">
        <Row>
          <Field theme={theme} label="Mode" width={160}><Select theme={theme} value={mode} onChange={setMode} options={['single', 'logie', 'curve']} /></Field>
          <Field theme={theme} label="Matrix / unit" width={300}>
            <Select theme={theme} value={matrix} onChange={setMatrix}
              options={Object.entries(info?.matrices || { water: '' }).map(([k, u]) => ({ value: k, label: `${k} — ${u}` }))} />
          </Field>
          <Field theme={theme} label="Uncertainty factor (×/÷)" width={170}><TextInput theme={theme} type="number" value={unc} onChange={setUnc} /></Field>
          <Field theme={theme} label="Levels to quantify" width={300}><LevelChips theme={theme} levels={['1', '2a', '2b', '3', '4']} value={levels} onChange={setLevels} /></Field>
        </Row>
        <p style={{ color: theme.muted, fontSize: '0.88rem' }}>{MODE_TEXT[mode]}</p>
        {mode === 'single' && (
          <Row>
            <Field theme={theme} label="Reference compound"><TextInput theme={theme} value={ref.name} onChange={(v) => setRef({ ...ref, name: v })} /></Field>
            <Field theme={theme} label="Reference conc. (µg/mL extract)" width={200}><TextInput theme={theme} type="number" value={ref.conc_ug_mL} onChange={(v) => setRef({ ...ref, conc_ug_mL: v })} /></Field>
            <Field theme={theme} label="Reference peak area (this software)" width={240}
              hint="The default comes from Compound Discoverer — re-measure it with ‘Measure areas’ below.">
              <TextInput theme={theme} type="number" value={ref.area} onChange={(v) => setRef({ ...ref, area: v })} />
            </Field>
          </Row>
        )}
        <ErrorBox theme={theme} error={err} />
        <Btn theme={theme} primary disabled={!runDir || busy} onClick={calc} style={{ marginTop: 10 }}>{busy ? 'Calculating…' : 'Calculate concentrations'}</Btn>
      </Panel>
      {res && (
        <Panel theme={theme} title={`3 · Concentrations (${res.unit})`}
          right={<Row><span style={{ color: theme.muted, fontSize: '0.85rem' }}>Saved: {res.path}</span>
            <Btn theme={theme} onClick={() => downloadCsv(`semiquant_${mode}_${matrix}.csv`, res.columns, res.rows)}>CSV</Btn></Row>}>
          {plotRows.length > 0 && (
            <Plot
              data={[{ type: 'scatter', mode: 'markers', y: plotRows.map((r) => `${r.name} (${r.sample?.slice(-12) ?? ''})`), x: plotRows.map((r) => r[c]),
                error_x: { type: 'data', symmetric: false, array: plotRows.map((r) => r[`${c}_high`] - r[c]), arrayminus: plotRows.map((r) => r[c] - r[`${c}_low`]) },
                marker: { size: 9, color: theme.accent } }]}
              layout={{ height: Math.max(260, plotRows.length * 22 + 80), margin: { l: 280, r: 20, t: 30, b: 50 }, xaxis: { type: 'log', title: res.unit },
                yaxis: { autorange: 'reversed', automargin: true }, title: { text: `Top ${plotRows.length} (point = estimate, bar = ×/÷ ${unc})`, font: { size: 13 } },
                paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { color: theme.text } }}
              config={{ displaylogo: false, responsive: true }} style={{ width: '100%' }} />
          )}
          <DataTable theme={theme} columns={res.columns.filter((x) => x !== 'smiles')} rows={res.rows} levelColumn="level" />
          <p style={{ color: theme.muted, fontSize: '0.85rem' }}>
            Absolute semi-quantification is not yet validated for this software (see engine/VALIDATION.md). Ionisation efficiency differs by orders of magnitude
            between compound classes; keep the uncertainty factor and compare against targeted results where available.
          </p>
        </Panel>
      )}
      <CalibrantEditor theme={theme} onUseArea={(r) => { setMode('single'); setRef({ name: r.name, conc_ug_mL: r.conc_ug_ml, area: r.area }); }} />
    </Page>
  );
}
