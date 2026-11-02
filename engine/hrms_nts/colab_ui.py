"""ipywidgets front-ends for the Colab notebook: sample sheet, calibration list, chromatogram explorer."""
from __future__ import annotations
import numpy as np
import pandas as pd

TYPES = ["sample", "blank", "qc", "cal", "exclude"]


def _w():
    import ipywidgets as w
    from IPython.display import display
    return w, display


def sample_sheet_editor(P, sheet: pd.DataFrame | None = None):
    """One row per mzML: type, group (replicates of one sample share a group), replicate, conc (calibration level)."""
    w, display = _w()
    s = sheet if sheet is not None else (P.sheet if P.sheet is not None else P.guess_sample_sheet())
    rows = []
    head = w.HBox([w.HTML(f"<b>{t}</b>", layout=w.Layout(width=wd)) for t, wd in
                   (("file", "300px"), ("type", "110px"), ("group (replicates share it)", "240px"), ("replicate", "90px"), ("cal. conc", "110px"))])
    for r in s.itertuples():
        rows.append(dict(file=r.file,
                         type=w.Dropdown(options=TYPES, value=r.type if r.type in TYPES else "sample", layout=w.Layout(width="110px")),
                         group=w.Text(value=str(r.group), layout=w.Layout(width="240px")),
                         replicate=w.Text(value=str(getattr(r, "replicate", "")), layout=w.Layout(width="90px")),
                         conc=w.FloatText(value=float(r.conc) if pd.notna(r.conc) else 0.0, layout=w.Layout(width="110px"))))
    box = w.VBox([head] + [w.HBox([w.Label(x["file"], layout=w.Layout(width="300px")), x["type"], x["group"], x["replicate"], x["conc"]]) for x in rows])
    btn = w.Button(description="Save sample sheet", button_style="success"); out = w.Output()

    def save(_):
        df = pd.DataFrame([dict(file=x["file"], type=x["type"].value, group=x["group"].value, replicate=x["replicate"].value,
                                conc=x["conc"].value if x["type"].value == "cal" else np.nan) for x in rows])
        df = df[df.type != "exclude"]
        with out:
            out.clear_output()
            try:
                P.set_sample_sheet(df)
                print("saved:", df.type.value_counts().to_dict())
                g = df[df.type == "sample"].groupby("group").size()
                print("sample groups (replicate injections):", g.to_dict())
            except Exception as e:
                print("ERROR:", e)
    btn.on_click(save)
    display(box, btn, out)


def calibration_editor(P, n_rows=6):
    """Parents (targets) and their internal standards; concentrations come from the sample sheet ('cal' rows)."""
    w, display = _w()
    rows = []
    box = w.VBox()

    def add_row(_=None, name="", smiles="", label="", role="target", isn=""):
        r = dict(name=w.Text(value=name, placeholder="Carbamazepine", layout=w.Layout(width="220px")),
                 smiles=w.Text(value=smiles, placeholder="SMILES of the unlabelled compound", layout=w.Layout(width="300px")),
                 formula=w.Text(value="", placeholder="or formula, e.g. C15H2[2H]10N2O", layout=w.Layout(width="220px")),
                 label=w.Text(value=label, placeholder="D10 / 13C3 (IS)", layout=w.Layout(width="110px")),
                 role=w.Dropdown(options=["target", "IS"], value=role, layout=w.Layout(width="90px")),
                 is_name=w.Text(value=isn, placeholder="IS used for this target", layout=w.Layout(width="220px")))
        rows.append(r)
        box.children = list(box.children) + [w.HBox([r["name"], r["smiles"], r["formula"], r["label"], r["role"], r["is_name"]])]

    head = w.HBox([w.HTML(f"<b>{t}</b>", layout=w.Layout(width=wd)) for t, wd in
                   (("compound", "220px"), ("SMILES", "300px"), ("formula (optional)", "220px"), ("label", "110px"), ("role", "90px"), ("internal standard", "220px"))])
    for _ in range(n_rows):
        add_row()
    add = w.Button(description="+ row"); add.on_click(add_row)
    out = w.Output()

    def table():
        recs = []
        for r in rows:
            if not r["name"].value.strip():
                continue
            recs.append(dict(name=r["name"].value.strip(), smiles=r["smiles"].value.strip(), formula=r["formula"].value.strip(),
                             label=r["label"].value.strip(), role=r["role"].value, is_name=r["is_name"].value.strip()))
        return pd.DataFrame(recs)

    display(head, box, add, out)
    return table


def chromatogram_explorer(P, min_level="3"):
    """Select a compound and an injection -> Q1/Q3 chromatograms, integration, peak rating, presence confidence."""
    w, display = _w()
    C = P.compounds(min_level)
    if C.empty:
        print("no identifications yet - run the screening and QA steps first"); return
    names = sorted(C.name.dropna().unique(), key=str.lower)
    samples = list(P.mzml_files())
    cn = w.Combobox(options=names, placeholder="type a compound", ensure_option=True, layout=w.Layout(width="420px"))
    sm = w.Dropdown(options=samples, layout=w.Layout(width="300px"))
    nq = w.IntSlider(value=6, min=1, max=10, description="Q3 ions")
    win = w.FloatSlider(value=1.0, min=0.3, max=3.0, step=0.1, description="± min")
    btn = w.Button(description="Show", button_style="primary"); out = w.Output()

    def upd_samples(ch):
        found = sorted(C.loc[C.name == cn.value, "sample"].unique())
        sm.options = found + [s for s in samples if s not in found]

    def show(_):
        with out:
            out.clear_output()
            if not cn.value:
                print("choose a compound"); return
            try:
                r = P.chromatogram(cn.value, sm.value, n_q3=nq.value, window=win.value, show=True)
                s = r["summary"]
                print(f"area {s['area']:,.0f} | S/N {s['sn']} | peak rating {s['peak_rating']}/10 ({s['integration']}) | "
                      f"Q1 {s['q1_ppm']} ppm | presence confidence {s['presence_confidence_pct']} % | Q3 ions with signal {s['n_q3_with_signal']}")
                if s["ion_ratios"]:
                    print("ion ratios vs the first Q3 ion:", s["ion_ratios"])
                print("product ions:", "; ".join(f"{x['label']} [{x['source']}]" for x in r["trace"]["q3"]))
            except Exception as e:
                print("ERROR:", e)
    cn.observe(upd_samples, names="value"); btn.on_click(show)
    display(w.HBox([cn, sm]), w.HBox([nq, win, btn]), out)


def show_calibration(res, max_plots=12):
    """Inline plotly panels: calibration levels, linear (dashed) and quadratic (dotted) fits, chosen model solid."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from .calibration import predict
    fits = res["fits"]; names = sorted(fits)[:max_plots]
    if not names:
        print("no calibration curves"); return
    cols = 3; rows = (len(names) + cols - 1) // cols
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=[n[:35] for n in names])
    for k, n in enumerate(names):
        F = fits[n]; ch = F["chosen"]; r, c = k // cols + 1, k % cols + 1
        xs = np.linspace(0, ch["x"].max() * 1.05, 100)
        fig.add_trace(go.Scatter(x=ch["x"], y=ch["y"], mode="markers", marker=dict(color="#2a78d6", size=8), name="levels", showlegend=k == 0), r, c)
        fig.add_trace(go.Scatter(x=xs, y=predict(F["linear"], xs), mode="lines", line=dict(color="#555", dash="dash"),
                                 name="linear", showlegend=k == 0), r, c)
        if F["quadratic"] is not None:
            fig.add_trace(go.Scatter(x=xs, y=predict(F["quadratic"], xs), mode="lines", line=dict(color="#eb6834", dash="dot"),
                                     name="quadratic", showlegend=k == 0), r, c)
        fig.add_trace(go.Scatter(x=xs, y=predict(ch, xs), mode="lines", line=dict(color="#1baf7a", width=3), opacity=0.5,
                                 name="chosen", showlegend=k == 0), r, c)
    fig.update_layout(height=320 * rows, template="plotly_white", title=f"Calibration curves ({res.get('units', '')})")
    fig.show()
