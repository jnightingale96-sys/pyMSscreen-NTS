"""Run provenance: everything needed to reproduce / audit a result folder (versions, git commit, inputs, merged config)."""
from __future__ import annotations
import hashlib, json, os, platform, subprocess, sys, datetime
import yaml
from . import __version__


def _sha_partial(path, n=4 << 20):
    try:
        h = hashlib.sha256(); size = os.path.getsize(path)
        with open(path, "rb") as f:
            h.update(f.read(n))
            if size > 2 * n:
                f.seek(-n, 2); h.update(f.read(n))
        return h.hexdigest()[:16]
    except OSError:
        return None


def _git(root):
    def run(*a):
        try:
            return subprocess.run(["git", "-C", root, *a], capture_output=True, text=True, timeout=10).stdout.strip()
        except Exception:
            return ""
    c = run("rev-parse", "--short", "HEAD")
    return dict(commit=c or None, dirty=bool(run("status", "--porcelain")) if c else None, branch=run("rev-parse", "--abbrev-ref", "HEAD") or None)


def write_provenance(out_dir: str, cfg: dict, files: list[str], argv: list[str] | None = None, extra: dict | None = None):
    pk = {}
    for m in ("numpy", "pandas", "scipy", "pyteomics", "rdkit", "sklearn", "matplotlib", "openpyxl"):
        try:
            mod = __import__(m); pk[m] = getattr(mod, "__version__", "?")
        except Exception:
            pk[m] = None
    def finfo(p):
        if not p or not os.path.exists(p): return None
        st = os.stat(p); return dict(path=os.path.abspath(p), size=st.st_size, mtime=datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"), sha256_partial=_sha_partial(p))
    from .common import resolve
    P = cfg["paths"]
    data = {k: finfo(resolve(cfg, P.get(k))) for k in ("susdat", "pubchemlite", "user_compounds", "library_db", "standards", "calibrants")}
    prov = dict(toolkit_version=__version__, timestamp=datetime.datetime.now().isoformat(timespec="seconds"), python=sys.version.split()[0],
                platform=platform.platform(), packages=pk, git=_git(cfg["_root"]), command=argv or sys.argv,
                config_layers=cfg.get("_config_layers"), data_sources=data, inputs=[finfo(f) for f in files], **(extra or {}))
    with open(os.path.join(out_dir, "provenance.json"), "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=2, default=str)
    clean = {k: v for k, v in cfg.items() if not k.startswith("_")}
    with open(os.path.join(out_dir, "config_used.yaml"), "w", encoding="utf-8") as f:
        yaml.safe_dump(clean, f, sort_keys=False)
    return prov
