"""Two-panel figure: (a) PEC (bar = property range) vs MEC per chemical, MS1-only vs MS2-supported grain detections; (b) grain EIC apex intensity by identification status.
usage: python tools/ms1only_figure.py <pec_screen_out_dir> <png>"""
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

out, png = sys.argv[1], sys.argv[2]; MAT = sys.argv[3] if len(sys.argv) > 3 else "grain"
S = pd.read_csv(os.path.join(out, "pec_screen_summary.csv")); P = pd.read_csv(os.path.join(out, "pec_vs_mec_by_field.csv")); D = pd.read_csv(os.path.join(out, "detections_all_injections.csv"))
S = S[~S["class"].str.startswith(("endogenous", "polymer", "identity conflict", "fate model"))]
S["ev"] = np.where(S.grain_identity == "MS1 only", "MS1 only", np.where(S.grain_identity.isin(["MS2 confirmed", "MS2 weak"]), "MS2 confirmed / weak", ""))
gm = lambda v: float(np.exp(np.nanmean(np.log(v))))
P = P[P.MEC_ng_g.notna() & P.PEC_central.notna()]
rows = []
for cid, g in P.groupby("cand_id"):
    s = S[S.cand_id == cid]
    if len(s) and s.ev.iloc[0]:
        rows.append(dict(cand_id=cid, ev=s.ev.iloc[0], mec=gm(g.MEC_ng_g), pc=gm(g.PEC_central), pl=gm(g.PEC_low), ph=gm(g.PEC_high)))
R = pd.DataFrame(rows)
fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5.6), dpi=150, gridspec_kw=dict(width_ratios=[1.35, 1]))
lim = [10 ** np.floor(np.log10(min(R.mec.min(), R.pl.min()))), 10 ** np.ceil(np.log10(max(R.mec.max(), R.ph.max())))]
xs = np.array(lim); a.plot(xs, xs, color="#333", lw=1); a.fill_between(xs, xs / 10, xs * 10, color="#999", alpha=0.18, lw=0)
for ev, col, mk, fc in (("MS2 confirmed / weak", "#1b9aaa", "o", "#1b9aaa"), ("MS1 only", "#d9822b", "o", "white")):
    r = R[R.ev == ev]
    a.errorbar(r.mec, r.pc, yerr=[r.pc - r.pl, r.ph - r.pc], fmt="none", ecolor=col, alpha=0.25, lw=0.8)
    a.scatter(r.mec, r.pc, s=26, marker=mk, facecolor=fc, edgecolor=col, linewidth=1.2, label=f"{ev} (n={len(r)})", zorder=3)
a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel("MEC in " + MAT + ", ng/g dw (semi-quantified, geo-mean over fields)"); a.set_ylabel("PEC in " + MAT + ", ng/g dw (dot = central, line = property range)")
a.set_title("PEC vs MEC per chemical (non-endogenous, non-polymer)", loc="left", fontsize=10.5); a.legend(frameon=False, loc="upper left", fontsize=9)
a.text(lim[1] / 30, lim[1] / 30 * 14, "grey band = within x10", fontsize=8, color="#555", rotation=33, ha="right")
g = D[(D.matrix == "grain") & (D.role == "sample") & (D.area > 0)]
order = ["MS2 confirmed", "MS2 weak", "MS2 conflict", "MS1 only"]
data = [np.log10(g[g.id_status == k].apex.values) for k in order]
bp = b.boxplot(data, vert=False, tick_labels=[f"{k}\n(n={len(d)})" for k, d in zip(order, data)], showfliers=False, patch_artist=True, medianprops=dict(color="#222"))
for p_, k in zip(bp["boxes"], order): p_.set_facecolor("#f3c89a" if k == "MS1 only" else "#a9dbe1")
b.set_xlabel(MAT + " EIC apex intensity at the matched peak (log10 counts)"); b.set_title("MS1-only peaks are weaker (MS2 rarely triggered)", loc="left", fontsize=10.5)
b.invert_yaxis()
for ax in (a, b): ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(png); print(png, len(R), R.ev.value_counts().to_dict())
