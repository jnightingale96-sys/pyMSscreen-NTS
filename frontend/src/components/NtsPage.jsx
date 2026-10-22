/**
 * @file NtsPage.jsx
 * @summary Non-target screening with the HRMS NTS engine: run Level 5 -> 1 identification on mzML files,
 * then browse level summary, compounds, candidates, mirror plots, and hand identified compounds to the
 * pyMSscreen prescreen/plotting workflow.
 */
import { useEffect, useMemo, useState } from 'react';
import Plot from 'react-plotly.js';
import {
  ntsInfo, ntsListMzml, ntsStartRun, ntsResults, ntsPlots, ntsFileUrl, ntsExportToPrescreen, ntsStructureUrl, apiError,
} from '../api/api';
import {
  useTheme, useSticky, Page, Panel, Btn, Field, TextInput, Select, Row, ErrorBox, LevelChips, DataTable, RunPicker,
  JobMonitor, downloadCsv, LEVEL_COLORS,
} from './nts/ui';

const LEVELS = ['1', '2a', '2b', '3', '4', '5'];
const LEVEL_TEXT = {
  '1': 'Authentic standard: RT within ±0.2 min and MS2 library match',
  '2a': 'Library spectrum: cosine ≥ 0.70, ≥ 2 matched fragments, same adduct, RT not contradicting',
  '2b': 'Diagnostic evidence: ≥ 3 decoy-adjusted explained product ions, clear top candidate (no library spectrum)',
  '3': 'Tentative: MS2 with ≥ 2 explained ions or weak library match, RT plausible',
  '4': 'Unequivocal molecular formula / exact-mass match (≤ 5 ppm)',
  '5': 'Exact mass of interest, no candidate',
};
const COMPOUND_COLS = ['best_level', 'name', 'formula', 'confidence_score', 'n_samples_detected', 'rt_apex', 'ppm_error', 'max_int', 'area'];

function EngineStatus({ theme, info, error }) {
  const [open, setOpen] = useState(false);
  if (error) return <ErrorBox theme={theme} error={`Engine not available: ${error}`} />;
  if (!info) return null;
  const missing = info.checks.filter((c) => !c.found && !['paths.user_compounds', 'paths.inhouse_msp'].includes(c.item));
  return (
    <Panel theme={theme} title={`Engine: hrms_nts ${info.version}`}
      right={<Btn theme={theme} onClick={() => setOpen(!open)}>{open ? 'Hide' : 'Show'} data sources</Btn>}>
      <div style={{ color: missing.length ? theme.dangerText : theme.ok }}>
        {missing.length ? `Missing: ${missing.map((m) => m.item).join(', ')} (edit engine/config.local.yaml)` : 'All data sources found.'}
      </div>
      {open && (
        <table style={{ marginTop: 10, fontSize: '0.85rem', borderCollapse: 'collapse' }}>
          <tbody>
            {info.checks.map((c) => (
              <tr key={c.item}>
                <td style={{ padding: '2px 10px 2px 0', color: c.found ? theme.ok : theme.muted }}>{c.found ? '✔' : '—'}</td>
                <td style={{ padding: '2px 10px 2px 0', fontWeight: 600 }}>{c.item}</td>
                <td style={{ padding: 2, color: theme.muted, wordBreak: 'break-all' }}>{c.path || 'not set'}</td>
              </tr>
            ))}
            <tr><td /><td style={{ fontWeight: 600 }}>config</td><td style={{ color: theme.muted }}>{info.config_layers.join('  ←  ')}</td></tr>
          </tbody>
        </table>
      )}
    </Panel>
  );
}

function RunPanel({ theme, info, onFinished }) {
  const [folder, setFolder] = useSticky('mzmlFolder', '');
  const [files, setFiles] = useState([]);
  const [chosen, setChosen] = useState([]);
  const [opts, setOpts] = useSticky('runOpts', { matrix: 'water', polarity: 'auto', semiquant: 'single', insilico: 'builtin', ppm: '', min_intensity: '', no_library: false, overrides: '' });
  const [err, setErr] = useState('');
  const [jobId, setJobId] = useState(null);
  const [running, setRunning] = useState(false);
  const set = (k) => (v) => setOpts((o) => ({ ...o, [k]: v }));

  const list = async () => {
    setErr('');
    try {
      const r = await ntsListMzml(folder);
      setFiles(r.data.files);
      setChosen(r.data.files.map((f) => f.name));
      if (!r.data.files.length) setErr('No .mzML files in this folder. Convert Thermo .raw files first (msconvert or ThermoRawFileParser, centroided).');
    } catch (e) {
      setFiles([]);
      setErr(apiError(e));
    }
  };

  const start = async () => {
    setErr('');
    try {
      const r = await ntsStartRun({
        folder, files: chosen, ...opts,
        overrides: (opts.overrides || '').split('\n').map((s) => s.trim()).filter(Boolean),
      });
      setJobId(r.data.job_id);
      setRunning(true);
    } catch (e) {
      setErr(apiError(e));
    }
  };

  return (
    <Panel theme={theme} title="1 · Run non-target screening">
      <Row>
        <Field theme={theme} label="Folder with centroided mzML files (DDA: MS1 + MS2)">
          <TextInput theme={theme} value={folder} onChange={setFolder} placeholder="C:\\Users\\me\\data\\maize_mzML" />
        </Field>
        <Btn theme={theme} onClick={list}>List mzML files</Btn>
      </Row>
      {files.length > 0 && (
        <div style={{ margin: '12px 0', maxHeight: 200, overflow: 'auto', border: `1px solid ${theme.border}`, borderRadius: 6, padding: 8 }}>
          <label style={{ display: 'block', fontWeight: 600, marginBottom: 4 }}>
            <input type="checkbox" checked={chosen.length === files.length}
              onChange={(e) => setChosen(e.target.checked ? files.map((f) => f.name) : [])} /> All ({files.length})
          </label>
          {files.map((f) => (
            <label key={f.name} style={{ display: 'block', fontSize: '0.9rem' }}>
              <input type="checkbox" checked={chosen.includes(f.name)}
                onChange={(e) => setChosen((c) => (e.target.checked ? [...c, f.name] : c.filter((x) => x !== f.name)))} />{' '}
              {f.name} <span style={{ color: theme.muted }}>({f.size_mb} MB)</span>
            </label>
          ))}
        </div>
      )}
      <Row style={{ marginTop: 10 }}>
        <Field theme={theme} label="Matrix (concentration units)" width={230}>
          <Select theme={theme} value={opts.matrix} onChange={set('matrix')}
            options={Object.entries(info?.matrices || { water: '' }).map(([k, u]) => ({ value: k, label: `${k} — ${u}` }))} />
        </Field>
        <Field theme={theme} label="Polarity" width={140}>
          <Select theme={theme} value={opts.polarity} onChange={set('polarity')} options={['auto', 'positive', 'negative']} />
        </Field>
        <Field theme={theme} label="Semi-quant" width={140}>
          <Select theme={theme} value={opts.semiquant} onChange={set('semiquant')} options={info?.semiquant_modes || ['single']} />
        </Field>
        <Field theme={theme} label="In-silico (product ions)" width={170}>
          <Select theme={theme} value={opts.insilico} onChange={set('insilico')} options={info?.insilico_engines || ['builtin']} />
        </Field>
        <Field theme={theme} label={`Mass tolerance ppm (default ${info?.ppm_default ?? 5})`} width={170}>
          <TextInput theme={theme} type="number" value={opts.ppm} onChange={set('ppm')} placeholder={String(info?.ppm_default ?? 5)} />
        </Field>
        <Field theme={theme} label="Min. feature intensity" width={170}>
          <TextInput theme={theme} type="number" value={opts.min_intensity} onChange={set('min_intensity')} placeholder={String(info?.min_intensity_default ?? 100000)} />
        </Field>
      </Row>
      <details style={{ marginTop: 10 }}>
        <summary style={{ cursor: 'pointer', color: theme.muted }}>Advanced</summary>
        <Row style={{ marginTop: 8 }}>
          <label><input type="checkbox" checked={opts.no_library} onChange={(e) => set('no_library')(e.target.checked)} /> Skip spectral library (faster, no Level 2a/1)</label>
        </Row>
        <Field theme={theme} label="Config overrides, one per line (section.key=value), e.g. ms1.bg_occupancy=0.5">
          <textarea value={opts.overrides} onChange={(e) => set('overrides')(e.target.value)} rows={3}
            style={{ fontFamily: 'monospace', padding: 8, borderRadius: 6, border: `1px solid ${theme.border}`, background: theme.inputBg, color: theme.inputText }} />
        </Field>
      </details>
      <ErrorBox theme={theme} error={err} />
      <Row style={{ marginTop: 12 }}>
        <Btn theme={theme} primary disabled={!chosen.length || running} onClick={start}>
          Start screening ({chosen.length} file{chosen.length === 1 ? '' : 's'})
        </Btn>
        <span style={{ color: theme.muted, fontSize: '0.85rem' }}>Roughly 1–6 min per sample. Results go to a new NTS_out_… folder next to the mzML files.</span>
      </Row>
      <JobMonitor theme={theme} jobId={jobId} onDone={(j) => { setRunning(false); if (j.state === 'done') onFinished(folder, j.out); }} />
    </Panel>
  );
}

function CompoundDetail({ theme, runDir, row, plots }) {
  if (!row) return <p style={{ color: theme.muted }}>Click a row to see its structure and MS2 mirror plot.</p>;
  const lvl = String(row.best_level ?? row.level);
  const fid = row.feature_id ? row.feature_id.split('_').pop() : null;
  const name = String(row.name || '').slice(0, 30).replace(/[/ :]/g, (c) => (c === '/' ? '-' : c === ' ' ? '_' : ''));
  const matches = plots.filter((p) => (fid ? p.includes(`_${fid}_`) : p.includes(name)));
  return (
    <div>
      <h4 style={{ margin: '0 0 4px' }}>{row.name}</h4>
      <div style={{ color: theme.muted, fontSize: '0.85rem', marginBottom: 6 }}>
        <span style={{ background: LEVEL_COLORS[lvl], color: '#fff', borderRadius: 10, padding: '1px 8px', fontWeight: 700 }}>Level {lvl}</span>{' '}
        {LEVEL_TEXT[lvl]}
      </div>
      {row.smiles && <img alt="structure" src={ntsStructureUrl(row.smiles, 240)} style={{ background: '#fff', borderRadius: 6, border: `1px solid ${theme.border}` }} />}
      <table style={{ fontSize: '0.82rem', marginTop: 6 }}>
        <tbody>
          {['formula', 'adduct', 'mz_apex', 'ppm_error', 'rt_apex', 'std_rt_diff', 'lib_name', 'lib_score', 'lib_matched', 'bi_score', 'bi_expl_adj_n', 'mf_score', 'level_basis', 'smiles']
            .filter((k) => row[k] !== undefined && row[k] !== null && row[k] !== '')
            .map((k) => <tr key={k}><td style={{ color: theme.muted, paddingRight: 10 }}>{k}</td><td style={{ wordBreak: 'break-all' }}>{typeof row[k] === 'number' ? Number(row[k].toPrecision(6)) : String(row[k])}</td></tr>)}
        </tbody>
      </table>
      {matches.slice(0, 3).map((p) => (
        <a key={p} href={ntsFileUrl(runDir, `plots/${p}`)} target="_blank" rel="noreferrer">
          <img alt={p} src={ntsFileUrl(runDir, `plots/${p}`)} style={{ width: '100%', marginTop: 8, borderRadius: 6, background: '#fff' }} />
        </a>
      ))}
    </div>
  );
}

function ResultsPanel({ theme, root, setRoot, runDir, setRunDir, refreshKey }) {
  const [summary, setSummary] = useState([]);
  const [tab, setTab] = useSticky('resultsTab', 'compounds');
  const [levels, setLevels] = useSticky('resultsLevels', ['1', '2a', '2b', '3']);
  const [q, setQ] = useState('');
  const [data, setData] = useState({ columns: [], rows: [], total: 0 });
  const [plots, setPlots] = useState([]);
  const [selected, setSelected] = useState(null);
  const [err, setErr] = useState('');
  const [msg, setMsg] = useState('');

  useEffect(() => {
    if (!runDir) return;
    setSelected(null);
    ntsResults({ run_dir: runDir, table: 'levels' }).then((r) => setSummary(r.data.rows)).catch((e) => setErr(apiError(e)));
    ntsPlots(runDir).then((r) => setPlots(r.data.plots)).catch(() => setPlots([]));
  }, [runDir]);

  useEffect(() => {
    if (!runDir) return undefined;
    const t = setTimeout(() => {
      ntsResults({ run_dir: runDir, table: tab, levels: levels.join(','), q, limit: 3000 })
        .then((r) => { setData(r.data); setErr(''); })
        .catch((e) => setErr(apiError(e)));
    }, 250);
    return () => clearTimeout(t);
  }, [runDir, tab, levels, q]);

  const columns = useMemo(() => {
    if (tab === 'compounds') {
      return [...COMPOUND_COLS.filter((c) => data.columns.includes(c)), ...data.columns.filter((c) => c.startsWith('conc_') && !COMPOUND_COLS.includes(c))];
    }
    return data.columns.filter((c) => !['smiles', 'inchikey', 'level_basis'].includes(c));
  }, [data.columns, tab]);

  const totals = LEVELS.map((l) => summary.reduce((s, r) => s + (Number(r[l]) || 0), 0));

  const exportPrescreen = async () => {
    setMsg('');
    try {
      const r = await ntsExportToPrescreen({ run_dir: runDir, levels: levels.filter((l) => l !== '5') });
      setMsg(`Wrote ${r.data.n} compounds to ${r.data.path}. On the Prescreen page, set that folder as working directory and choose this CSV to extract EICs / MS2 and review them in Plotting.`);
    } catch (e) {
      setErr(apiError(e));
    }
  };

  return (
    <Panel theme={theme} title="2 · Results">
      <RunPicker theme={theme} root={root} setRoot={setRoot} runDir={runDir} setRunDir={setRunDir} refreshKey={refreshKey} />
      {runDir && (
        <>
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(340px, 1fr))', gap: 16, marginTop: 8 }}>
            <Plot
              data={[{ type: 'bar', x: LEVELS.map((l) => `L${l}`), y: totals, marker: { color: LEVELS.map((l) => LEVEL_COLORS[l]) }, text: totals.map(String), textposition: 'outside' }]}
              layout={{ title: { text: 'Features per confidence level (all samples)', font: { size: 14 } }, height: 260, margin: { t: 40, b: 30, l: 50, r: 10 },
                yaxis: { type: 'log', title: 'features' }, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { color: theme.text } }}
              config={{ displaylogo: false, responsive: true }} style={{ width: '100%' }} />
            <div>
              <DataTable theme={theme} columns={['sample', ...LEVELS]} rows={summary} maxHeight={240} />
              <Row style={{ marginTop: 10 }}>
                <a href={ntsFileUrl(runDir, 'NTS_results.xlsx', true)}><Btn theme={theme}>Download Excel workbook</Btn></a>
                <a href={ntsFileUrl(runDir, 'run.log')} target="_blank" rel="noreferrer"><Btn theme={theme}>run.log</Btn></a>
                <Btn theme={theme} onClick={exportPrescreen} title="Write identified compounds as a pyMSscreen compounds CSV">Send to pyMSscreen prescreen</Btn>
              </Row>
              {msg && <p style={{ color: theme.ok, fontSize: '0.85rem' }}>{msg}</p>}
            </div>
          </div>

          <Row style={{ marginTop: 14, alignItems: 'center' }}>
            {[['compounds', 'Compounds (best level)'], ['identified', 'Identified L1–2b'], ['top', 'All features (top candidate)']].map(([k, label]) => (
              <Btn key={k} theme={theme} primary={tab === k} onClick={() => setTab(k)}>{label}</Btn>
            ))}
            {tab !== 'identified' && <LevelChips theme={theme} levels={LEVELS} value={levels} onChange={setLevels} />}
            <div style={{ width: 220 }}><TextInput theme={theme} value={q} onChange={setQ} placeholder="Search name…" /></div>
            <span style={{ color: theme.muted }}>{data.total} row(s){data.total > data.rows.length ? `, showing ${data.rows.length}` : ''}</span>
            <Btn theme={theme} onClick={() => downloadCsv(`${tab}.csv`, data.columns, data.rows)}>CSV</Btn>
          </Row>
          <ErrorBox theme={theme} error={err} />
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: 16, marginTop: 10 }}>
            <div style={{ flex: '2 1 560px', minWidth: 0 }}>
            <DataTable theme={theme} columns={columns} rows={data.rows} levelColumn={tab === 'compounds' ? 'best_level' : 'level'}
              onRowClick={setSelected} rowKey={(r) => `${r.sample || ''}|${r.feature_id || ''}|${r.name}|${r.ik1 || ''}`}
              selectedKey={selected ? `${selected.sample || ''}|${selected.feature_id || ''}|${selected.name}|${selected.ik1 || ''}` : undefined} />
            </div>
            <div style={{ flex: '1 1 300px', border: `1px solid ${theme.border}`, borderRadius: 6, padding: 10, maxHeight: 620, overflow: 'auto' }}>
              <CompoundDetail theme={theme} runDir={runDir} row={selected} plots={plots} />
            </div>
          </div>
          <details style={{ marginTop: 12 }}>
            <summary style={{ cursor: 'pointer', color: theme.muted }}>What the levels mean</summary>
            <ul style={{ fontSize: '0.88rem' }}>{LEVELS.map((l) => <li key={l}><strong>Level {l}</strong>: {LEVEL_TEXT[l]}</li>)}</ul>
            <p style={{ fontSize: '0.85rem', color: theme.muted }}>Features carrying a background flag cannot reach 2b/3; neutral mass &lt; 120 Da cannot reach 2b/3 from in-silico evidence alone. confidence_score is a ranking aid, not a probability.</p>
          </details>
        </>
      )}
    </Panel>
  );
}

export default function NtsPage() {
  const theme = useTheme();
  const [info, setInfo] = useState(null);
  const [infoErr, setInfoErr] = useState('');
  const [root, setRoot] = useSticky('resultsRoot', '');
  const [runDir, setRunDir] = useSticky('runDir', '');
  const [refreshKey, setRefreshKey] = useState(0);

  useEffect(() => { ntsInfo().then((r) => setInfo(r.data)).catch((e) => setInfoErr(apiError(e))); }, []);

  return (
    <Page theme={theme} title="Non-target screening" subtitle="Level 5 → 4 → 3 → 2b → 2a → 1 identification with NORMAN SusDat, PubChemLite, MS/MS libraries, RT prediction and product-ion evidence.">
      <EngineStatus theme={theme} info={info} error={infoErr} />
      <RunPanel theme={theme} info={info} onFinished={(folder, out) => { setRoot(folder); setRunDir(out); setRefreshKey((k) => k + 1); }} />
      <ResultsPanel theme={theme} root={root} setRoot={setRoot} runDir={runDir} setRunDir={setRunDir} refreshKey={refreshKey} />
    </Page>
  );
}
