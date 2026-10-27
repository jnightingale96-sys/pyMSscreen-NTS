/**
 * @file TargetedPage.jsx
 * @summary Targeted screen with the HRMS NTS engine (`hrms_nts targeted`): internal standards, targets and calibrants.
 * Q1 (precursor) and Q4 (product-ion) mass errors, confidence level, peak integration QC, blank check, calibration curves
 * and RT anchors. Works on DDA and DIA/SWATH data, centroided or profile (SCIEX profile data are centroided on reading).
 */
import { useState } from 'react';
import { ntsListMzml, ntsStartTargeted, ntsTargetedResults, ntsFileUrl, apiError } from '../api/api';
import { useTheme, useSticky, Page, Panel, Btn, Field, TextInput, Select, Row, ErrorBox, DataTable, JobMonitor, downloadCsv } from './nts/ui';

const ROLES = ['sample', 'blank', 'qc', 'cal'];
const TABS = [
  ['Detected', 'Detected'], ['Blank_check', 'Blank check'], ['Not_detected', 'Not detected'], ['Calibration', 'Calibration'],
  ['Quantification', 'Quantification'], ['Native_at_IS_RT', 'Natives at IS RT'], ['RT_anchors', 'RT anchors'], ['Files', 'Files'], ['Notes', 'Notes'],
];
const DETECTED_COLS = ['File', 'Compound', 'Adduct', 'Q1 theoretical', 'Q1 error (ppm)', 'Q1 <= 5 ppm', 'RT apex (min)', 'Peak area (±10 ppm XIC)',
  'S/N', 'Peak rating', 'Blank as % of this peak', 'Confidence level', 'Q4 ions <= 10 ppm', 'Q4 <= 10 ppm', 'Library match', 'Isobaric group'];

function RunPanel({ theme, onFinished }) {
  const [folder, setFolder] = useSticky('targetedFolder', '');
  const [files, setFiles] = useState([]);
  const [roles, setRoles] = useState({});
  const [opts, setOpts] = useSticky('targetedOpts', { list: '', q1_ppm: '5', q4_ppm: '10', polarity: 'auto', min_apex: '', update_standards: false });
  const [err, setErr] = useState('');
  const [jobId, setJobId] = useState(null);
  const [running, setRunning] = useState(false);
  const set = (k) => (v) => setOpts((o) => ({ ...o, [k]: v }));
  const setRole = (name, k, v) => setRoles((r) => ({ ...r, [name]: { ...(r[name] || { use: true, type: 'sample', conc: '' }), [k]: v } }));

  const list = async () => {
    setErr('');
    try {
      const r = await ntsListMzml(folder);
      setFiles(r.data.files);
      setRoles(Object.fromEntries(r.data.files.map((f) => [f.name, { use: true, type: /bl(an)?k|eq\d|solvent/i.test(f.name) ? 'blank' : 'sample', conc: '' }])));
      if (!r.data.files.length) setErr('No .mzML files in this folder.');
    } catch (e) {
      setFiles([]); setErr(apiError(e));
    }
  };
  const chosen = files.map((f) => f.name).filter((n) => roles[n]?.use);
  const start = async () => {
    setErr('');
    try {
      const r = await ntsStartTargeted({ folder, files: chosen, ...opts,
        roles: Object.fromEntries(chosen.map((n) => [n, { type: roles[n].type, conc: roles[n].conc }])) });
      setJobId(r.data.job_id); setRunning(true);
    } catch (e) {
      setErr(apiError(e));
    }
  };
  return (
    <Panel theme={theme} title="1 · Screen internal standards / targets / calibrants">
      <Row>
        <Field theme={theme} label="Folder with mzML files (DDA or DIA/SWATH; centroid or profile)">
          <TextInput theme={theme} value={folder} onChange={setFolder} placeholder="C:\\data\\run_mzML" />
        </Field>
        <Btn theme={theme} onClick={list}>List mzML files</Btn>
      </Row>
      {files.length > 0 && (
        <div style={{ margin: '12px 0', maxHeight: 260, overflow: 'auto', border: `1px solid ${theme.border}`, borderRadius: 6, padding: 8 }}>
          <table style={{ fontSize: '0.9rem', borderCollapse: 'collapse' }}>
            <thead><tr><th /><th style={{ textAlign: 'left' }}>File</th><th>Role</th><th>Cal. conc.</th></tr></thead>
            <tbody>
              {files.map((f) => (
                <tr key={f.name}>
                  <td><input type="checkbox" checked={!!roles[f.name]?.use} onChange={(e) => setRole(f.name, 'use', e.target.checked)} /></td>
                  <td style={{ paddingRight: 12 }}>{f.name} <span style={{ color: theme.muted }}>({f.size_mb} MB)</span></td>
                  <td style={{ width: 110 }}><Select theme={theme} value={roles[f.name]?.type || 'sample'} onChange={(v) => setRole(f.name, 'type', v)} options={ROLES} /></td>
                  <td style={{ width: 110 }}>{roles[f.name]?.type === 'cal' && (
                    <TextInput theme={theme} type="number" value={roles[f.name]?.conc} onChange={(v) => setRole(f.name, 'conc', v)} placeholder="e.g. 10" />)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
      <Row style={{ marginTop: 10 }}>
        <Field theme={theme} label="Compound list (CSV/XLSX: name + formula with [2H]/[13C], or SMILES + label; optional Q1, rt, role, is_name)">
          <TextInput theme={theme} value={opts.list} onChange={set('list')} placeholder="C:\\data\\DrugIS_pos.xlsx  or  engine\\examples\\drug_IS_pos_75.csv" />
        </Field>
      </Row>
      <Row style={{ marginTop: 10 }}>
        <Field theme={theme} label="Q1 tolerance (ppm)" width={150}><TextInput theme={theme} type="number" value={opts.q1_ppm} onChange={set('q1_ppm')} /></Field>
        <Field theme={theme} label="Q4 tolerance (ppm)" width={150}><TextInput theme={theme} type="number" value={opts.q4_ppm} onChange={set('q4_ppm')} /></Field>
        <Field theme={theme} label="Polarity" width={140}><Select theme={theme} value={opts.polarity} onChange={set('polarity')} options={['auto', 'positive', 'negative']} /></Field>
        <Field theme={theme} label="Min. peak height (blank = instrument default)" width={260}>
          <TextInput theme={theme} type="number" value={opts.min_apex} onChange={set('min_apex')} placeholder="SCIEX 600 · Orbitrap 3000" />
        </Field>
        <label style={{ alignSelf: 'center' }}><input type="checkbox" checked={opts.update_standards} onChange={(e) => set('update_standards')(e.target.checked)} />{' '}
          Add detected IS RTs to standards.csv (RT models / Level 1)</label>
      </Row>
      <ErrorBox theme={theme} error={err} />
      <Row style={{ marginTop: 12 }}>
        <Btn theme={theme} primary disabled={!chosen.length || !opts.list || running} onClick={start}>Start targeted screen ({chosen.length} file{chosen.length === 1 ? '' : 's'})</Btn>
        <span style={{ color: theme.muted, fontSize: '0.85rem' }}>~2–4 min per file for 75 compounds. Results go to a new targeted_out_… folder next to the mzML files.</span>
      </Row>
      <JobMonitor theme={theme} jobId={jobId} onDone={(j) => { setRunning(false); if (j.state === 'done') onFinished(j.out); }} />
    </Panel>
  );
}

function ResultsPanel({ theme, runDir, setRunDir }) {
  const [data, setData] = useState(null);
  const [tab, setTab] = useSticky('targetedTab', 'Detected');
  const [err, setErr] = useState('');
  const load = async (dir) => {
    setErr('');
    try {
      const r = await ntsTargetedResults(dir || runDir);
      setData(r.data);
    } catch (e) {
      setData(null); setErr(apiError(e));
    }
  };
  const rows = data?.[tab] || [];
  const cols = tab === 'Detected' && rows.length ? DETECTED_COLS.filter((c) => c in rows[0]) : (rows.length ? Object.keys(rows[0]) : []);
  return (
    <Panel theme={theme} title="2 · Results">
      <Row>
        <Field theme={theme} label="Results folder (targeted_out_…)"><TextInput theme={theme} value={runDir} onChange={setRunDir} /></Field>
        <Btn theme={theme} onClick={() => load()}>Load</Btn>
      </Row>
      <ErrorBox theme={theme} error={err} />
      {data && (
        <>
          <Row style={{ margin: '12px 0' }}>
            {data.files.map((f) => <a key={f} href={ntsFileUrl(runDir, f, true)} style={{ color: theme.accent, fontWeight: 600 }}>⬇ {f}</a>)}
          </Row>
          <Row style={{ marginBottom: 10 }}>
            {TABS.filter(([k]) => data[k]?.length).map(([k, label]) => (
              <Btn key={k} theme={theme} primary={tab === k} onClick={() => setTab(k)}>{label} ({data[k].length})</Btn>
            ))}
            {rows.length > 0 && <Btn theme={theme} onClick={() => downloadCsv(`${tab}.csv`, Object.keys(rows[0]), rows)}>Download table</Btn>}
          </Row>
          <DataTable theme={theme} columns={cols} rows={rows} maxHeight={600} />
        </>
      )}
    </Panel>
  );
}

export default function TargetedPage() {
  const theme = useTheme();
  const [runDir, setRunDir] = useSticky('targetedRunDir', '');
  return (
    <Page theme={theme} title="Targeted screen — internal standards, targets, calibrants"
      subtitle="Q1 / Q4 mass errors, confidence level (1/2a/2b/3/4), peak integration QC, blank check, calibration curves and RT anchors for the RT models.">
      <RunPanel theme={theme} onFinished={(out) => setRunDir(out)} />
      <ResultsPanel theme={theme} runDir={runDir} setRunDir={setRunDir} />
    </Page>
  );
}
