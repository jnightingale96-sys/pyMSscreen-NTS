"""Calibration curves: weighted linear and quadratic (second-order polynomial) fits, model choice, LOD/LOQ, back-calculation,
quantification and plots.

Model choice (per compound): the quadratic term is kept only when the curve is genuinely curved - at least `min_levels_quadratic`
levels, a significant extra-sum-of-squares F-test (p < 0.05), lower AICc, and a monotonic curve over the calibrated range.
Otherwise the linear model is used. Both fits are always reported. Response = area, or area / IS area when an IS is assigned."""
from __future__ import annotations
import math
import numpy as np
import pandas as pd


def _weights(x, weight):
    x = np.asarray(x, float)
    return {"1/x": 1 / x, "1/x2": 1 / x ** 2, "1/x^2": 1 / x ** 2}.get(weight, np.ones_like(x))


def _wls(X, y, w):
    W = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * W[:, None], y * W, rcond=None)
    return beta


def _aicc(rss_w, n, k):
    if n - k - 1 <= 0 or rss_w <= 0:
        return np.nan
    return n * math.log(rss_w / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)


def fit(x, y, weight="1/x", model="linear"):
    """-> dict(model, coef (c0 + c1 x [+ c2 x^2]), r2, n, rss_w, aicc, sres) or None."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x, y = x[ok], y[ok]
    k = 2 if model == "linear" else 3
    if len(np.unique(x)) < k:
        return None
    w = _weights(x, weight)
    X = np.c_[np.ones_like(x), x] if model == "linear" else np.c_[np.ones_like(x), x, x ** 2]
    b = _wls(X, y, w)
    pred = X @ b
    rss_w = float((w * (y - pred) ** 2).sum())
    ss = ((y - pred) ** 2).sum(); tot = ((y - y.mean()) ** 2).sum()
    return dict(model=model, coef=[float(c) for c in b], r2=float(1 - ss / tot) if tot > 0 else np.nan, n=int(len(x)),
                rss_w=rss_w, aicc=_aicc(rss_w, len(x), k), sres=float(np.sqrt(ss / max(len(x) - k, 1))), x=x, y=y, weight=weight)


def predict(f, x):
    c = f["coef"]; x = np.asarray(x, float)
    return c[0] + c[1] * x + (c[2] * x ** 2 if len(c) == 3 else 0.0)


def invert(f, y, lo=None, hi=None):
    """Concentration for a response; for the quadratic the root inside (or nearest to) the calibrated range."""
    c = f["coef"]; y = np.atleast_1d(np.asarray(y, float)); out = np.full(len(y), np.nan)
    lo = f["x"].min() if lo is None else lo; hi = f["x"].max() if hi is None else hi
    for i, v in enumerate(y):
        if not np.isfinite(v):
            continue
        if len(c) == 2:
            out[i] = (v - c[0]) / c[1] if c[1] else np.nan
            continue
        a, b, cc = c[2], c[1], c[0] - v
        if abs(a) < 1e-15:
            out[i] = -cc / b if b else np.nan; continue
        disc = b * b - 4 * a * cc
        if disc < 0:
            continue
        r = np.array([(-b + math.sqrt(disc)) / (2 * a), (-b - math.sqrt(disc)) / (2 * a)])
        r = r[r >= 0] if (r >= 0).any() else r
        mid = (lo + hi) / 2
        out[i] = r[np.argmin(np.abs(r - mid))]
    return out


def _monotonic(f, lo, hi):
    xs = np.linspace(lo, hi, 50)
    d = np.diff(predict(f, xs))
    return bool((d > 0).all() or (d < 0).all())


def choose(x, y, weight="1/x", min_levels_quadratic=5, alpha=0.05):
    """Fit linear and quadratic; return (chosen, linear, quadratic, reason)."""
    lin = fit(x, y, weight, "linear")
    quad = fit(x, y, weight, "quadratic") if lin is not None else None
    if lin is None:
        return None, None, None, "fewer than 2 calibration levels"
    if quad is None or len(np.unique(lin["x"])) < min_levels_quadratic:
        return lin, lin, quad, f"linear (quadratic needs >= {min_levels_quadratic} levels)"
    n = lin["n"]
    F = (lin["rss_w"] - quad["rss_w"]) / (quad["rss_w"] / max(n - 3, 1)) if quad["rss_w"] > 0 else np.inf
    try:
        from scipy.stats import f as fdist
        p = float(1 - fdist.cdf(F, 1, max(n - 3, 1)))
    except Exception:
        p = np.nan
    mono = _monotonic(quad, lin["x"].min(), lin["x"].max())
    if np.isfinite(p) and p < alpha and (not np.isfinite(lin["aicc"]) or quad["aicc"] < lin["aicc"]) and mono:
        return quad, lin, quad, f"quadratic (curvature significant, F-test p={p:.3g}, AICc {quad['aicc']:.1f} < {lin['aicc']:.1f})"
    why = "not monotonic over the range" if not mono else f"curvature not significant (p={p:.3g})"
    return lin, lin, quad, f"linear ({why})"


def figures_of_merit(f):
    """LOD = 3.3 s_res / slope, LOQ = 10 s_res / slope (slope at the low end for the quadratic); back-calculated accuracy (%)."""
    c = f["coef"]
    slope = c[1] + (2 * c[2] * f["x"].min() if len(c) == 3 else 0.0)
    back = invert(f, f["y"])
    acc = back / f["x"] * 100
    lod = 3.3 * f["sres"] / slope if slope > 0 else np.nan
    return dict(slope_low=slope, lod=lod, loq=lod / 3.3 * 10 if np.isfinite(lod) else np.nan,
                accuracy_pct=";".join(f"{v:.0f}" for v in acc), accuracy_ok=bool(np.all(np.abs(acc - 100) <= 20)) if len(acc) else False)


def calibrate_table(cal: pd.DataFrame, weight="1/x", min_levels_quadratic=5):
    """cal: rows name, conc, response [, file]. -> (summary DataFrame, dict name -> chosen fit)."""
    rows, fits = [], {}
    for n, g in cal.groupby("name"):
        ch, lin, quad, why = choose(g["conc"].values, g["response"].values, weight, min_levels_quadratic)
        if ch is None:
            rows.append(dict(name=n, model="none", reason=why, n_levels=g["conc"].nunique())); continue
        fom = figures_of_merit(ch)
        warn = []
        if len(np.unique(ch["x"])) < 3: warn.append(f"only {len(np.unique(ch['x']))} levels - not a valid calibration (>= 3, better >= 5)")
        if fom["slope_low"] <= 0: warn.append("slope <= 0 - response does not rise with concentration")
        if np.isfinite(ch["r2"]) and ch["r2"] < 0.99: warn.append(f"R2 {ch['r2']:.3f} < 0.99")
        if not fom["accuracy_ok"]: warn.append("back-calculated accuracy outside 80-120 %")
        rows.append(dict(name=n, model=ch["model"], reason=why, weighting=weight, n_levels=int(len(np.unique(ch["x"]))),
                         range_low=float(ch["x"].min()), range_high=float(ch["x"].max()),
                         intercept=ch["coef"][0], slope=ch["coef"][1], quad_term=ch["coef"][2] if len(ch["coef"]) == 3 else 0.0, r2=ch["r2"],
                         r2_linear=lin["r2"], r2_quadratic=quad["r2"] if quad else np.nan, aicc_linear=lin["aicc"],
                         aicc_quadratic=quad["aicc"] if quad else np.nan, **fom, warning="; ".join(warn)))
        fits[n] = dict(chosen=ch, linear=lin, quadratic=quad)
    return pd.DataFrame(rows), fits


def quantify(fits, samples: pd.DataFrame):
    """samples: rows name, file, response -> conc + flag (below LOQ / outside range / no curve)."""
    out = []
    for r in samples.itertuples():
        F = fits.get(r.name)
        if F is None:
            out.append(dict(name=r.name, file=r.file, response=r.response, conc=np.nan, flag="no calibration curve")); continue
        f = F["chosen"]; c = float(invert(f, [r.response])[0])
        fom = figures_of_merit(f)
        flag = ("below LOQ" if np.isfinite(fom["loq"]) and c < fom["loq"] else "above range" if c > f["x"].max() * 1.1
                else "below range" if c < f["x"].min() * 0.9 else "")
        out.append(dict(name=r.name, file=r.file, response=r.response, conc=c, model=f["model"], flag=flag))
    return pd.DataFrame(out)


def plot_curves(fits, path, units="", max_per_page=9):
    """One panel per compound: points, linear (dashed) and quadratic (dotted) fits, chosen model solid; residual inset text."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    names = sorted(fits)
    with PdfPages(path) as pdf:
        for k in range(0, max(len(names), 1), max_per_page):
            fig, axs = plt.subplots(3, 3, figsize=(11.7, 8.3)); axs = axs.ravel()
            for ax, n in zip(axs, names[k:k + max_per_page]):
                F = fits[n]; ch = F["chosen"]
                xs = np.linspace(0, ch["x"].max() * 1.05, 200)
                ax.plot(ch["x"], ch["y"], "o", color="#2a78d6", ms=5, label="levels")
                ax.plot(xs, predict(F["linear"], xs), "--", color="#555", lw=1, label=f"linear R2 {F['linear']['r2']:.4f}")
                if F["quadratic"] is not None:
                    ax.plot(xs, predict(F["quadratic"], xs), ":", color="#eb6834", lw=1.4, label=f"quadratic R2 {F['quadratic']['r2']:.4f}")
                ax.plot(xs, predict(ch, xs), "-", color="#1baf7a", lw=1.6, alpha=0.6, label=f"chosen: {ch['model']}")
                ax.set_title(n[:40], fontsize=8, loc="left"); ax.set_xlabel(f"concentration {units}", fontsize=7); ax.set_ylabel("response", fontsize=7)
                ax.tick_params(labelsize=6); ax.legend(fontsize=6, frameon=False)
                for s in ("top", "right"): ax.spines[s].set_visible(False)
            for ax in axs[len(names[k:k + max_per_page]):]:
                ax.axis("off")
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    return path
