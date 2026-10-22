/**
 * @file ui.jsx
 * @summary Shared building blocks for the HRMS NTS pages (theme, panels, inputs, sortable table,
 * results-folder picker, job monitor). Same palette as the original pyMSscreen pages.
 */
import { useEffect, useMemo, useState } from 'react';
import { ntsRuns, ntsJob, ntsCancel, apiError } from '../../api/api';

export const usePrefersDark = () => {
  const mq = typeof window !== 'undefined' && window.matchMedia ? window.matchMedia('(prefers-color-scheme: dark)') : null;
  const [isDark, setIsDark] = useState(mq ? mq.matches : false);
  useEffect(() => {
    if (!mq) return;
    const onChange = (e) => setIsDark(e.matches);
    mq.addEventListener('change', onChange);
    return () => mq.removeEventListener('change', onChange);
  }, [mq]);
  return isDark;
};

export const makeTheme = (isDark) =>
  isDark
    ? { bg: '#0f172a', panel: '#0b1220', panelLite: '#111827', text: '#e5e7eb', muted: '#9ca3af', border: '#334155',
        tableHeader: '#111827', inputBg: '#0b1220', inputText: '#e5e7eb', btnBg: '#1f2937', btnText: '#e5e7eb',
        accent: '#5daeff', dangerText: '#ff7b7b', ok: '#4ade80', rowHover: '#1e293b' }
    : { bg: '#e0f4f7', panel: '#ffffff', panelLite: '#f0faff', text: '#111827', muted: '#374151', border: '#D1D5DB',
        tableHeader: '#f2f2f2', inputBg: '#ffffff', inputText: '#111827', btnBg: '#a7cce5', btnText: '#111827',
        accent: '#0077ff', dangerText: '#c00', ok: '#15803d', rowHover: '#eef6fb' };

export const useTheme = () => makeTheme(usePrefersDark());

/** useState remembered in localStorage (folders, last results folder) so pages share them across visits. */
export function useSticky(key, initial) {
  const [v, setV] = useState(() => {
    try {
      const s = localStorage.getItem(`nts.${key}`);
      return s === null ? initial : JSON.parse(s);
    } catch {
      return initial;
    }
  });
  useEffect(() => {
    try { localStorage.setItem(`nts.${key}`, JSON.stringify(v)); } catch { /* storage unavailable */ }
  }, [key, v]);
  return [v, setV];
}

export const LEVEL_COLORS = { '1': '#15803d', '2a': '#22c55e', '2b': '#84cc16', '3': '#eab308', '4': '#f97316', '5': '#9ca3af' };

export function Page({ theme, title, subtitle, children }) {
  return (
    <div style={{ background: theme.bg, minHeight: '100vh', color: theme.text, padding: '20px 16px 60px' }}>
      <div style={{ maxWidth: 1400, margin: '0 auto' }}>
        <h1 style={{ fontSize: '2rem', fontWeight: 800, margin: '6px 0 2px' }}>{title}</h1>
        {subtitle && <p style={{ color: theme.muted, margin: '0 0 18px' }}>{subtitle}</p>}
        {children}
      </div>
    </div>
  );
}

export function Panel({ theme, title, right, children, style }) {
  return (
    <section style={{ background: theme.panel, border: `1px solid ${theme.border}`, borderRadius: 10, padding: 16, marginBottom: 16, ...style }}>
      {(title || right) && (
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', gap: 12, marginBottom: 12, flexWrap: 'wrap' }}>
          {title && <h3 style={{ margin: 0, fontSize: '1.15rem' }}>{title}</h3>}
          {right}
        </div>
      )}
      {children}
    </section>
  );
}

export function Btn({ theme, primary, danger, disabled, onClick, children, title, style }) {
  return (
    <button type="button" disabled={disabled} onClick={onClick} title={title}
      style={{ padding: '8px 14px', borderRadius: 6, border: `1px solid ${theme.border}`, fontWeight: 600,
        background: primary ? theme.accent : theme.btnBg, color: primary ? '#fff' : (danger ? theme.dangerText : theme.btnText),
        opacity: disabled ? 0.5 : 1, cursor: disabled ? 'not-allowed' : 'pointer', ...style }}>
      {children}
    </button>
  );
}

const inputStyle = (theme) => ({ padding: '7px 9px', borderRadius: 6, border: `1px solid ${theme.border}`, background: theme.inputBg,
  color: theme.inputText, fontSize: '0.95rem', boxSizing: 'border-box', width: '100%' });

export function Field({ theme, label, hint, children, width }) {
  return (
    <label style={{ display: 'flex', flexDirection: 'column', gap: 4, minWidth: width || 160, flex: width ? `0 0 ${width}px` : 1 }}>
      <span style={{ fontSize: '0.85rem', fontWeight: 600, color: theme.muted }}>{label}</span>
      {children}
      {hint && <span style={{ fontSize: '0.78rem', color: theme.muted }}>{hint}</span>}
    </label>
  );
}

export function TextInput({ theme, value, onChange, placeholder, type = 'text', style }) {
  return <input type={type} value={value ?? ''} placeholder={placeholder} onChange={(e) => onChange(e.target.value)}
    style={{ ...inputStyle(theme), ...style }} />;
}

export function Select({ theme, value, onChange, options }) {
  return (
    <select value={value} onChange={(e) => onChange(e.target.value)} style={inputStyle(theme)}>
      {options.map((o) => (typeof o === 'string'
        ? <option key={o} value={o}>{o}</option>
        : <option key={o.value} value={o.value}>{o.label}</option>))}
    </select>
  );
}

export function Row({ children, gap = 12, style }) {
  return <div style={{ display: 'flex', gap, flexWrap: 'wrap', alignItems: 'flex-end', ...style }}>{children}</div>;
}

export function ErrorBox({ theme, error }) {
  if (!error) return null;
  return <div style={{ color: theme.dangerText, border: `1px solid ${theme.dangerText}`, borderRadius: 6, padding: '8px 10px', margin: '10px 0', whiteSpace: 'pre-wrap' }}>{error}</div>;
}

export function LevelChips({ theme, levels, value, onChange }) {
  const toggle = (l) => onChange(value.includes(l) ? value.filter((x) => x !== l) : [...value, l]);
  return (
    <div style={{ display: 'flex', gap: 6, flexWrap: 'wrap' }}>
      {levels.map((l) => {
        const on = value.includes(l);
        return (
          <button key={l} type="button" onClick={() => toggle(l)}
            style={{ padding: '4px 10px', borderRadius: 14, fontWeight: 700, cursor: 'pointer', border: `2px solid ${LEVEL_COLORS[l]}`,
              background: on ? LEVEL_COLORS[l] : 'transparent', color: on ? '#fff' : theme.text }}>
            L{l}
          </button>
        );
      })}
    </div>
  );
}

const fmt = (v, col = '') => {
  if (v === null || v === undefined) return '';
  if (typeof v === 'number') {
    if (!Number.isFinite(v)) return '';
    if (/(^|_)mz($|_)|exact_mass/.test(col)) return v.toFixed(4);
    if (/rt/.test(col) && !/rt_(agree|conflict|available)/.test(col)) return v.toFixed(2);
    if (Number.isInteger(v)) return Math.abs(v) >= 1e6 ? v.toExponential(2) : String(v);
    const a = Math.abs(v);
    if (a !== 0 && (a >= 1e6 || a < 1e-3)) return v.toExponential(2);
    return v.toFixed(a >= 100 ? 1 : a >= 1 ? 3 : 4);
  }
  if (typeof v === 'boolean') return v ? 'yes' : '';
  return String(v);
};

/** Sortable, scrollable table. `onRowClick(row)` optional; `selected` highlights a row by key. */
export function DataTable({ theme, columns, rows, onRowClick, selectedKey, rowKey, maxHeight = 520, levelColumn }) {
  const [sort, setSort] = useState({ col: null, asc: true });
  const sorted = useMemo(() => {
    if (!sort.col) return rows;
    const r = [...rows];
    r.sort((a, b) => {
      const x = a[sort.col], y = b[sort.col];
      if (x === y) return 0;
      if (x === null || x === undefined) return 1;
      if (y === null || y === undefined) return -1;
      const c = typeof x === 'number' && typeof y === 'number' ? x - y : String(x).localeCompare(String(y), undefined, { numeric: true });
      return sort.asc ? c : -c;
    });
    return r;
  }, [rows, sort]);
  if (!rows?.length) return <p style={{ color: theme.muted }}>No rows.</p>;
  return (
    <div style={{ overflow: 'auto', maxHeight, border: `1px solid ${theme.border}`, borderRadius: 6 }}>
      <table style={{ borderCollapse: 'collapse', width: '100%', fontSize: '0.85rem' }}>
        <thead>
          <tr>
            {columns.map((c) => (
              <th key={c} onClick={() => setSort((s) => ({ col: c, asc: s.col === c ? !s.asc : true }))}
                style={{ position: 'sticky', top: 0, background: theme.tableHeader, color: theme.text, textAlign: 'left',
                  padding: '6px 8px', borderBottom: `1px solid ${theme.border}`, cursor: 'pointer', whiteSpace: 'nowrap', zIndex: 1 }}>
                {c}{sort.col === c ? (sort.asc ? ' ▲' : ' ▼') : ''}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {sorted.map((r, i) => {
            const key = rowKey ? rowKey(r) : i;
            const sel = selectedKey !== undefined && key === selectedKey;
            return (
              <tr key={key} onClick={onRowClick ? () => onRowClick(r) : undefined}
                style={{ cursor: onRowClick ? 'pointer' : 'default', background: sel ? theme.rowHover : 'transparent' }}>
                {columns.map((c) => {
                  const v = r[c];
                  const isLvl = c === levelColumn;
                  return (
                    <td key={c} title={typeof v === 'string' && v.length > 40 ? v : undefined}
                      style={{ padding: '5px 8px', borderBottom: `1px solid ${theme.border}`, whiteSpace: 'nowrap',
                        maxWidth: 320, overflow: 'hidden', textOverflow: 'ellipsis' }}>
                      {isLvl
                        ? <span style={{ background: LEVEL_COLORS[v] || theme.border, color: '#fff', borderRadius: 10, padding: '1px 8px', fontWeight: 700 }}>{v}</span>
                        : fmt(v, c)}
                    </td>
                  );
                })}
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

/** Pick a results folder: type a parent folder, list result folders inside it, choose one. */
export function RunPicker({ theme, root, setRoot, runDir, setRunDir, refreshKey }) {
  const [runs, setRuns] = useState([]);
  const [err, setErr] = useState('');
  const load = async (folder) => {
    setErr('');
    if (!folder) return;
    try {
      const r = await ntsRuns(folder);
      setRuns(r.data.runs || []);
      if (!r.data.runs?.length) setErr('No results folders (containing level_summary.csv) found here.');
      else if (!r.data.runs.find((x) => x.run_dir === runDir)) setRunDir(r.data.runs[0].run_dir);
    } catch (e) {
      setErr(apiError(e));
    }
  };
  useEffect(() => { load(root); }, [refreshKey]); // eslint-disable-line react-hooks/exhaustive-deps
  return (
    <>
      <Row>
        <Field theme={theme} label="Folder containing results (NTS_out_* folders)">
          <TextInput theme={theme} value={root} onChange={setRoot} placeholder="C:\\data\\my_mzML_folder" />
        </Field>
        <Btn theme={theme} onClick={() => load(root)}>Find results</Btn>
        <Field theme={theme} label="Results folder" width={420}>
          <Select theme={theme} value={runDir || ''} onChange={setRunDir}
            options={[{ value: '', label: runs.length ? '— choose —' : '(none found)' },
              ...runs.map((r) => ({ value: r.run_dir, label: `${r.name}   (${r.modified.replace('T', ' ')})` }))]} />
        </Field>
      </Row>
      <ErrorBox theme={theme} error={err} />
    </>
  );
}

/** Polls a background engine job and shows progress + log tail. */
export function JobMonitor({ theme, jobId, onDone }) {
  const [job, setJob] = useState(null);
  useEffect(() => {
    if (!jobId) return undefined;
    let stop = false;
    let timer;
    const tick = async () => {
      try {
        const r = await ntsJob(jobId);
        if (stop) return;
        setJob(r.data);
        if (r.data.state === 'running') timer = setTimeout(tick, 2000);
        else onDone?.(r.data);
      } catch {
        if (!stop) timer = setTimeout(tick, 4000);
      }
    };
    tick();
    return () => { stop = true; clearTimeout(timer); };
  }, [jobId]); // eslint-disable-line react-hooks/exhaustive-deps
  if (!jobId || !job) return null;
  const colour = job.state === 'done' ? theme.ok : job.state === 'running' ? theme.accent : theme.dangerText;
  return (
    <div style={{ marginTop: 12 }}>
      <Row style={{ alignItems: 'center' }}>
        <strong style={{ color: colour }}>
          {job.state === 'running' ? 'Running…' : job.state === 'done' ? 'Finished' : job.state === 'cancelled' ? 'Cancelled' : 'Failed'}
        </strong>
        {job.n_files > 0 && <span style={{ color: theme.muted }}>{job.samples_done} of {job.n_files} sample(s) processed{job.samples_failed ? `, ${job.samples_failed} failed` : ''}</span>}
        {job.state === 'running' && <Btn theme={theme} danger onClick={() => ntsCancel(jobId)}>Cancel</Btn>}
      </Row>
      {job.progress !== null && job.progress !== undefined && (
        <div style={{ height: 10, background: theme.panelLite, border: `1px solid ${theme.border}`, borderRadius: 6, marginTop: 8, overflow: 'hidden' }}>
          <div style={{ width: `${job.progress}%`, height: '100%', background: colour, transition: 'width .5s' }} />
        </div>
      )}
      {job.warnings?.length > 0 && (
        <div style={{ color: theme.dangerText, fontSize: '0.85rem', marginTop: 8 }}>
          <strong>Warnings — check before using these results:</strong>
          <ul style={{ margin: '4px 0', paddingLeft: 18 }}>
            {job.warnings.map((w, i) => <li key={i}>{w.replace(/^\S+ +(WARNING|ERROR) +/, '')}</li>)}
          </ul>
        </div>
      )}
      <pre style={{ background: theme.panelLite, border: `1px solid ${theme.border}`, borderRadius: 6, padding: 8, marginTop: 8,
        maxHeight: 220, overflow: 'auto', fontSize: '0.78rem', whiteSpace: 'pre-wrap' }}>
        {(job.log || []).join('\n') || 'Starting…'}
      </pre>
    </div>
  );
}

/** Client-side CSV download of table rows. */
export function downloadCsv(filename, columns, rows) {
  const esc = (v) => {
    if (v === null || v === undefined) return '';
    const s = String(v);
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  const text = [columns.join(','), ...rows.map((r) => columns.map((c) => esc(r[c])).join(','))].join('\n');
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([text], { type: 'text/csv' }));
  a.download = filename;
  a.click();
  URL.revokeObjectURL(a.href);
}
