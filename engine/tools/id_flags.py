"""Attach identification-confidence flags to a pec-screen summary:
  * rt_err_min / rt_plausible : retention time vs the RTI->RT map fitted on the robust Level 1/2a anchors (tools/rt_check_levels.py, work/rt_check_2b.csv);
                                anchors themselves are 'anchor'.  plausible = |error| <= 2 x anchor leave-one-out RMSE (default 2.1 min)
  * metfrag_rank_all / metfrag_drug_n : rank of the assigned structure among all SusDat+PubChemLite isomers (FragmenterScore + MassBank score) and the number of
                                structures of the drug-rich NORMAN collection at that mass (work/metfrag_verify.csv)
Writes <run>/pec_screen_summary_flagged.csv and prints how many of the chart chemicals survive.

    python tools/id_flags.py work/pec_screen_israel_pos [--tol 2.1]"""
import argparse, os, sys
import numpy as np, pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--tol", type=float, default=2.1); a = ap.parse_args()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); W = os.path.join(ROOT, "work")
S = pd.read_csv(os.path.join(a.run, "pec_screen_summary.csv"))
rt = pd.read_csv(os.path.join(W, "rt_check_2b.csv"))[["name", "err"]].rename(columns={"err": "rt_err_min"}).drop_duplicates("name")
S = S.merge(rt, on="name", how="left")
S["rt_status"] = np.where(S.rt_err_min.isna(), np.where(S.level.isin(["1", "2a"]), "anchor level (not tested)", "no RTI"), np.where(S.rt_err_min.abs() <= a.tol, "plausible", "IMPLAUSIBLE"))
mf_path = os.path.join(W, "metfrag_verify.csv")
if os.path.exists(mf_path):
    M = pd.read_csv(mf_path)[["cand_id", "all_n", "all_rank", "all_score", "all_top1", "drug_n", "drug_rank"]].rename(columns={"all_n": "metfrag_n_isomers", "all_rank": "metfrag_rank_all", "all_score": "metfrag_score",
                                                                                                                       "all_top1": "metfrag_top1", "drug_n": "metfrag_drug_collection_n", "drug_rank": "metfrag_drug_rank"})
    S = S.merge(M, on="cand_id", how="left")
S["id_flag"] = np.where(S.rt_status == "IMPLAUSIBLE", "RT implausible", "")
out = os.path.join(a.run, "pec_screen_summary_flagged.csv"); S.to_csv(out, index=False)
chart = S[(S.n_pairs > 0) & ~S["class"].str.startswith(("unconstrained", "identity", "endogenous"))]
print(f"chart chemicals: {len(chart)}; RT implausible among them: {(chart.rt_status == 'IMPLAUSIBLE').sum()}")
cols = [c for c in ["name", "level", "class", "rt_err_min", "rt_status", "metfrag_rank_all", "metfrag_n_isomers", "metfrag_top1", "metfrag_drug_collection_n"] if c in chart.columns]
pd.set_option("display.width", 250); pd.set_option("display.float_format", lambda v: f"{v:,.2f}")
print(chart[cols].assign(name=chart.name.str[:34], **{"class": chart["class"].str[:22]}, **({"metfrag_top1": chart.metfrag_top1.astype(str).str[:30]} if "metfrag_top1" in chart.columns else {})).to_string(index=False))
print("->", out)
