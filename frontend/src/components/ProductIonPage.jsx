/**
 * @file ProductIonPage.jsx
 * @summary Product-ion prediction for a structure (built-in fragmenter: bond cleavages + neutral losses),
 * suggested MRM/PRM transitions, and optional scoring against an experimental MS2 spectrum.
 */
import { useState } from 'react';
import Plot from 'react-plotly.js';
import { ntsPredictIons, ntsStructureUrl, apiError } from '../api/api';
import { useTheme, useSticky, Page, Panel, Btn, Field, TextInput, Row, ErrorBox, DataTable, downloadCsv } from './nts/ui';

const ADDUCTS = { positive: ['[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+'], negative: ['[M-H]-', '[M+Cl]-', '[M+HCOO]-', '[M+CH3COO]-'] };

function SpectrumPlot({ theme, result }) {
  const ions = result.ions || [];
  const maxScore = Math.max(...ions.map((i) => i.rank_score), 1e-9);
  const traces = [{
    type: 'bar', x: ions.map((i) => i.mz), y: ions.map((i) => (-100 * i.rank_score) / maxScore), width: 0.6, name: 'predicted (rank)',
    marker: { color: ions.map((i) => (i.kind === 'neutral_loss' ? '#f97316' : '#5daeff')) },
    text: ions.map((i) => `${i.mz.toFixed(4)}  ${i.label}`), hoverinfo: 'text',
  }];
  const exp = result.experimental?.peaks;
  if (exp?.length) {
    const mx = Math.max(...exp.map((p) => p[1]));
    traces.unshift({ type: 'bar', x: exp.map((p) => p[0]), y: exp.map((p) => (100 * p[1]) / mx), width: 0.6, name: 'experimental',
      marker: { color: theme.text }, hovertemplate: '%{x:.4f}  %{y:.0f}%<extra></extra>' });
  }
  return (
    <Plot data={traces}
      layout={{ height: 330, barmode: 'overlay', margin: { t: 30, l: 60, r: 10, b: 45 }, showlegend: true, legend: { orientation: 'h', y: 1.12 },
        xaxis: { title: 'm/z', range: [0, result.precursor_mz + 10] },
        yaxis: { title: exp?.length ? 'experimental ↑ / predicted ↓' : 'predicted rank', zeroline: true },
        shapes: [{ type: 'line', x0: result.precursor_mz, x1: result.precursor_mz, y0: -100, y1: 100, line: { dash: 'dot', color: theme.muted } }],
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { color: theme.text } }}
      config={{ displaylogo: false, responsive: true }} style={{ width: '100%' }} />
  );
}

export default function ProductIonPage() {
  const theme = useTheme();
  const [name, setName] = useSticky('piName', 'Carbamazepine');
  const [smiles, setSmiles] = useSticky('piSmiles', 'NC(=O)N1c2ccccc2C=Cc2ccccc12');
  const [adducts, setAdducts] = useSticky('piAdducts', ['[M+H]+']);
  const [top, setTop] = useState(20);
  const [peaks, setPeaks] = useState('');
  const [res, setRes] = useState(null);
  const [err, setErr] = useState('');
  const [busy, setBusy] = useState(false);

  const toggle = (a) => setAdducts((v) => (v.includes(a) ? v.filter((x) => x !== a) : [...v, a]));
  const run = async () => {
    setErr(''); setBusy(true);
    try {
      const r = await ntsPredictIons({ smiles, adducts, top, experimental_peaks: peaks });
      setRes(r.data);
    } catch (e) {
      setErr(apiError(e)); setRes(null);
    } finally {
      setBusy(false);
    }
  };

  return (
    <Page theme={theme} title="Product-ion prediction" subtitle="Predicted MS2 fragments and neutral losses for a structure — for compounds with no library spectrum (e.g. NORMAN suspects), method building and checking an experimental spectrum.">
      <Panel theme={theme} title="Structure">
        <Row>
          <Field theme={theme} label="Name (for file names)" width={220}><TextInput theme={theme} value={name} onChange={setName} /></Field>
          <Field theme={theme} label="SMILES"><TextInput theme={theme} value={smiles} onChange={setSmiles} /></Field>
          <Field theme={theme} label="Ions to show" width={110}><TextInput theme={theme} type="number" value={top} onChange={(v) => setTop(Number(v) || 20)} /></Field>
        </Row>
        <div style={{ display: 'flex', gap: 30, marginTop: 12, flexWrap: 'wrap' }}>
          {Object.entries(ADDUCTS).map(([pol, list]) => (
            <div key={pol}>
              <div style={{ fontWeight: 600, color: theme.muted, fontSize: '0.85rem' }}>{pol}</div>
              {list.map((a) => (
                <label key={a} style={{ marginRight: 12 }}><input type="checkbox" checked={adducts.includes(a)} onChange={() => toggle(a)} /> {a}</label>
              ))}
            </div>
          ))}
        </div>
        <Field theme={theme} label="Optional: experimental MS2 peaks to compare (one 'm/z intensity' pair per line)">
          <textarea value={peaks} onChange={(e) => setPeaks(e.target.value)} rows={4} placeholder={'194.0964 100\n192.0808 35\n179.0730 12'}
            style={{ fontFamily: 'monospace', padding: 8, borderRadius: 6, border: `1px solid ${theme.border}`, background: theme.inputBg, color: theme.inputText }} />
        </Field>
        <ErrorBox theme={theme} error={err} />
        <Btn theme={theme} primary onClick={run} disabled={!smiles || !adducts.length || busy} style={{ marginTop: 10 }}>{busy ? 'Predicting…' : 'Predict product ions'}</Btn>
      </Panel>

      {res && (
        <Panel theme={theme} title={`${name || 'Compound'} — ${res.formula}, monoisotopic mass ${res.exact_mass.toFixed(5)}`}>
          <div style={{ display: 'flex', gap: 16, flexWrap: 'wrap', alignItems: 'flex-start' }}>
            <img alt="structure" src={ntsStructureUrl(res.smiles, 260)} style={{ background: '#fff', borderRadius: 8, border: `1px solid ${theme.border}` }} />
            <p style={{ color: theme.muted, fontSize: '0.85rem', maxWidth: 520 }}>
              Built-in predictor: even-electron product ions from one or two bond cleavages, plus common neutral losses. It does not model
              rearrangements. Benchmark on MassBank: top-1 ≈ 60 %, AUC 0.80–0.84. Treat predicted transitions as suggestions and confirm with a
              standard; experimentally derived ions take precedence.
            </p>
          </div>
          {res.results.map((r) => (
            <div key={r.adduct} style={{ borderTop: `1px solid ${theme.border}`, marginTop: 14, paddingTop: 10 }}>
              <h4 style={{ margin: '0 0 6px' }}>{r.adduct}{r.precursor_mz ? `  ·  precursor m/z ${r.precursor_mz.toFixed(4)}` : ''}</h4>
              {r.error ? <ErrorBox theme={theme} error={r.error} /> : (
                <>
                  {r.experimental && (
                    <p style={{ margin: '4px 0' }}>
                      Experimental spectrum: <b>{r.experimental.bi_expl_n}</b> of {r.experimental.bi_n_considered} peaks explained
                      (decoy-adjusted {r.experimental.bi_expl_adj_n.toFixed(1)}); explained intensity {(100 * r.experimental.bi_expl_frac).toFixed(0)} %
                      vs decoy {(100 * r.experimental.bi_decoy_frac).toFixed(0)} % → score <b>{r.experimental.bi_score.toFixed(2)}</b>
                      <span style={{ color: theme.muted }}> (Level 2b needs ≥ 3 adjusted ions and score ≥ 0.20)</span>
                    </p>
                  )}
                  <SpectrumPlot theme={theme} result={r} />
                  <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(320px, 1fr))', gap: 16 }}>
                    <div>
                      <Row style={{ marginBottom: 6 }}>
                        <strong>Predicted ions</strong>
                        <Btn theme={theme} onClick={() => downloadCsv(`${name || 'compound'}_${r.adduct}_product_ions.csv`, ['mz', 'label', 'kind', 'rank_score'], r.ions)}>CSV</Btn>
                      </Row>
                      <DataTable theme={theme} columns={['mz', 'label', 'kind', 'rank_score']} rows={r.ions} maxHeight={320} />
                    </div>
                    <div>
                      <strong>Suggested transitions (MRM/PRM)</strong>
                      <DataTable theme={theme} columns={['precursor_mz', 'mz', 'label', 'kind']} rows={r.transitions} maxHeight={320} />
                    </div>
                  </div>
                </>
              )}
            </div>
          ))}
        </Panel>
      )}
    </Page>
  );
}
