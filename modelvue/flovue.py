"""
modelvue.flovue - Read a MODFLOW model with FloPy and show its configuration.

Supports MODFLOW 6 simulations (every model in mfsim.nam: GWF, GWT, GWE, PRT ...)
and classic name-file models: MODFLOW-2005, MODFLOW-NWT, MODFLOW-USG,
MODFLOW-2000, MT3DMS / MT3D-USGS and SEAWAT. Everything is read through FloPy
(https://flopy.readthedocs.io), so any package FloPy can load is shown, with
a generic summary of its options, arrays and stress-period data.

Shows: grid and georeference, time discretisation, solver settings, every
package, FloPy's model check, and - if the model has been run - the listing
file (termination, convergence, water budget and discrepancy), head and budget
files. The web page adds figures drawn with FloPy's own plotting (needs
matplotlib).

Requires: flopy, numpy. Optional: matplotlib (figures), pandas (FloPy uses it
for list budgets).

Usage
-----
    flovue  path/to/mfsim.nam                  # MODFLOW 6 simulation
    flovue  path/to/model.nam                  # MF2005 / NWT / USG / MT3D / SEAWAT
    flovue  path/to/folder                     # finds every simulation / name file
    flovue  folder --model sagehen             # pick one when several are found
    flovue  model.nam --version mfnwt          # force the MODFLOW version
    flovue  mfsim.nam --all                    # terminal: every dataset of every package
    flovue  mfsim.nam --no-check               # skip FloPy's model check
    flovue  mfsim.nam --json config.json
    flovue  mfsim.nam --html [report.html] [--no-browser]
"""
from __future__ import annotations

import argparse
import base64
import contextlib
import datetime as _dt
import io
import json
import math
import re
import sys
import warnings
from pathlib import Path

try:
    import numpy as np
    import flopy
    from flopy.utils import Util2d, Util3d, Transient2d, Transient3d, MfList
except ImportError as _err:  # pragma: no cover
    raise ImportError(f"modelvue.flovue needs flopy and numpy ({_err}). "
                      "Install with: pip install 'modelvue[modflow]'") from _err

from ._common import _blk, _e, h1, h2, _kv, row, _table, _two

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Package names
# --------------------------------------------------------------------------- #
PKG_NAMES = {
    # MODFLOW 6
    "tdis": "Temporal discretisation", "ims": "Iterative model solution", "ems": "Explicit model solution",
    "pts": "Particle tracking solution", "gwfgwf": "GWF-GWF exchange", "gwfgwt": "GWF-GWT exchange",
    "gwfgwe": "GWF-GWE exchange", "gwfprt": "GWF-PRT exchange", "mvr": "Water mover", "gnc": "Ghost-node correction",
    "dis": "Structured discretisation", "disv": "Vertex discretisation", "disu": "Unstructured discretisation",
    "dis2d": "2D discretisation", "disv1d": "1D vertex discretisation",
    "ic": "Initial conditions", "npf": "Node property flow", "sto": "Storage", "csub": "Skeletal storage, compaction & subsidence",
    "buy": "Buoyancy", "vsc": "Viscosity", "hfb": "Horizontal flow barrier", "oc": "Output control", "obs": "Observations",
    "chd": "Constant head", "wel": "Well", "drn": "Drain", "riv": "River", "ghb": "General head boundary",
    "rch": "Recharge", "rcha": "Recharge (array)", "evt": "Evapotranspiration", "evta": "Evapotranspiration (array)",
    "maw": "Multi-aquifer well", "sfr": "Streamflow routing", "lak": "Lake", "uzf": "Unsaturated zone flow",
    "api": "API package", "adv": "Advection", "dsp": "Dispersion", "ssm": "Source-sink mixing", "mst": "Mobile storage & transfer",
    "ist": "Immobile storage & transfer", "src": "Mass source", "cnc": "Constant concentration", "fmi": "Flow model interface",
    "sft": "Streamflow transport", "lkt": "Lake transport", "mwt": "Multi-aquifer well transport", "uzt": "UZF transport",
    "mvt": "Mover transport", "cnd": "Conduction & dispersion", "est": "Energy storage & transfer", "ctp": "Constant temperature",
    "esl": "Energy source loading", "prp": "Particle release point", "mip": "Model input (PRT)",
    # classic
    "bas6": "Basic", "bcf6": "Block-centred flow", "lpf": "Layer property flow", "upw": "Upstream weighting",
    "huf2": "Hydrogeologic-unit flow", "nwt": "Newton solver", "pcg": "Preconditioned conjugate gradient",
    "pcgn": "PCG with improved nonlinear control", "gmg": "Geometric multigrid", "sip": "Strongly implicit procedure",
    "de4": "Direct solver", "sms": "Sparse matrix solver", "uzf": "Unsaturated zone flow", "sfr": "Streamflow routing",
    "str": "Stream", "lak": "Lake", "mnw2": "Multi-node well", "mnwi": "MNW information", "fhb": "Flow & head boundary",
    "evt": "Evapotranspiration", "ets": "ET segments", "sub": "Subsidence", "swt": "Subsidence & aquifer-system compaction",
    "swi2": "Seawater intrusion", "gage": "Gage", "hob": "Head observations", "lmt6": "Link to MT3D",
    "zone": "Zone arrays", "mult": "Multiplier arrays", "pval": "Parameter values", "list": "Listing file",
    "btn": "Basic transport", "rct": "Reactions",
    "gcg": "Generalised conjugate gradient", "vdf": "Variable-density flow", "vsc": "Viscosity", "disu": "Unstructured discretisation",
    "cln": "Connected linear network", "gnc": "Ghost-node correction",
}
STRESS_TYPES = {"chd", "wel", "drn", "riv", "ghb", "rch", "rcha", "evt", "evta", "maw", "sfr", "lak", "uzf", "cnc", "src",
                "ctp", "esl", "fhb", "str", "mnw2", "drt", "ets", "ssm", "prp", "api"}


def pkg_tag(p) -> str:
    t, n = p["type"].upper(), (p.get("name") or "").upper()
    return n if n and (n.startswith(t) or n == t) else (f"{t} ({n})" if n else t)


def pkg_label(ptype: str) -> str:
    t = ptype.lower().rstrip("0123456789") if ptype.lower() not in PKG_NAMES else ptype.lower()
    return PKG_NAMES.get(ptype.lower()) or PKG_NAMES.get(t) or ""


# --------------------------------------------------------------------------- #
# Generic helpers
# --------------------------------------------------------------------------- #
def g(v, nd=4):
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if math.isnan(f):
        return "nan"
    if f.is_integer() and abs(f) < 1e12:
        return str(int(f))
    return f"{f:.{nd}g}"


def fmt_time(t, units=None):
    if t is None:
        return None
    u = f" {units}" if units and units not in ("unknown", "undefined", "None") else ""
    return f"{g(t)}{u}"


@contextlib.contextmanager
def quiet():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        yield buf


def arr_stats(a, mask_big=True) -> dict | None:
    try:
        a = np.asarray(a, dtype=float)
    except (TypeError, ValueError):
        return None
    if a.size == 0:
        return None
    v = a[np.isfinite(a)]
    if mask_big:
        v = v[np.abs(v) < 1e20]
    if v.size == 0:
        return {"shape": a.shape, "n": 0}
    return {"shape": a.shape, "n": int(v.size), "min": float(v.min()), "max": float(v.max()), "mean": float(v.mean()),
            "constant": bool(v.min() == v.max())}


def stats_text(s, unit="") -> str:
    if not s:
        return ""
    shp = "x".join(str(x) for x in s["shape"])
    if not s.get("n"):
        return f"[{shp}] no valid values"
    if s["constant"]:
        return f"[{shp}] constant {g(s['min'])}{unit}"
    return f"[{shp}] {g(s['min'])} .. {g(s['max'])}{unit} (mean {g(s['mean'])})"


def rec_summary(ra, maxcols=8) -> str:
    """Summarise a numpy recarray: rows, numeric ranges per column."""
    try:
        n = len(ra)
    except TypeError:
        return repr(ra)[:120]
    if n == 0:
        return "0 rows"
    names = list(ra.dtype.names or [])
    parts = []
    for c in names[:maxcols]:
        col = ra[c]
        if col.dtype.kind in "fiu":
            s = arr_stats(col)
            if s and s.get("n"):
                parts.append(f"{c} {g(s['min'])}" + ("" if s["constant"] else f"..{g(s['max'])}"))
        elif c.lower() in ("cellid", "cellid1", "cellid2"):
            parts.append(f"{c} {len(set(map(str, col)))} unique")
        else:
            vals = sorted({str(x) for x in col})
            parts.append(f"{c} {', '.join(vals[:3])}{' ...' if len(vals) > 3 else ''}")
    more = f" (+{len(names) - maxcols} cols)" if len(names) > maxcols else ""
    return f"{n} rows; " + "; ".join(parts) + more


def period_summary(d: dict, nper=None) -> str:
    """dict {kper: data} -> text."""
    keys = sorted(k for k in d.keys() if isinstance(k, (int, np.integer)))
    if not keys:
        return ""
    first = d[keys[0]]
    if isinstance(first, np.ndarray) and first.dtype.names is None:
        body = stats_text(arr_stats(first))
    elif hasattr(first, "dtype") and getattr(first.dtype, "names", None):
        body = rec_summary(first)
    elif first is None:
        body = "none"
    else:
        body = repr(first)[:100]
    span = f"period {keys[0] + 1}" if len(keys) == 1 else f"{len(keys)} periods ({keys[0] + 1} .. {keys[-1] + 1})"
    sizes = [len(d[k]) for k in keys if hasattr(d[k], "dtype") and getattr(d[k].dtype, "names", None)]
    size = f", rows/period {min(sizes)}..{max(sizes)}" if sizes and min(sizes) != max(sizes) else ""
    return f"{span}{size}; first: {body}"


# --------------------------------------------------------------------------- #
# Discovery
# --------------------------------------------------------------------------- #
def _nam_entries(path: Path) -> list[tuple[str, str, str]]:
    out = []
    try:
        txt = path.read_text(errors="replace")
    except OSError:
        return out
    for ln in txt.splitlines():
        s = ln.split("#", 1)[0].strip()
        if not s:
            continue
        t = s.split()
        if len(t) >= 3 and re.fullmatch(r"-?\d+", t[1]):
            out.append((t[0].upper(), t[1], t[2]))
        elif t[0].upper() == "BEGIN":
            return [("BEGIN", "", "")]
    return out


def classic_version(path: Path) -> str:
    ft = {e[0] for e in _nam_entries(path)}
    if ft & {"VDF", "VSC"}:
        return "seawat"
    if "BTN" in ft:
        return "mt3d-usgs" if ft & {"SFT", "UZT", "LKT", "CTS", "TSO"} else "mt3dms"
    if ft & {"DISU", "SMS", "CLN"}:
        return "mfusg"
    if ft & {"NWT", "UPW"}:
        return "mfnwt"
    try:
        head = "\n".join(ln for ln in path.read_text(errors="replace").splitlines() if ln.lstrip().startswith("#"))
    except OSError:
        head = ""
    if re.search(r"\bUSG\b", head, re.I):
        return "mfusg"
    if re.search(r"NWT", head, re.I):
        return "mfnwt"
    for e in _nam_entries(path):              # NWT-style OPTIONS blocks in UZF / SFR / WEL / LAK ...
        if e[0] in ("UZF", "SFR", "WEL", "LAK", "MNW2", "AG"):
            f = path.parent / e[2]
            try:
                first = [ln.strip() for ln in f.read_text(errors="replace").splitlines()[:6]
                         if ln.strip() and not ln.lstrip().startswith("#")]
            except OSError:
                continue
            if first and first[0].upper().startswith("OPTIONS"):
                return "mfnwt"
    if "GLOBAL" in ft or "GLO" in ft:
        return "mf2k"
    return "mf2005"


def discover(arg: str) -> list[dict]:
    p = Path(arg).expanduser().resolve()
    if not p.exists():
        sys.exit(f"Path not found: {p}")
    if p.is_file():
        if p.name.lower() == "mfsim.nam":
            return [{"kind": "mf6", "path": p}]
        e = _nam_entries(p)
        if e and e[0][0] == "BEGIN":          # an MF6 model name file -> its simulation
            sim = p.parent / "mfsim.nam"
            if sim.exists():
                return [{"kind": "mf6", "path": sim, "only_model_nam": p.name}]
            sys.exit(f"{p.name} is a MODFLOW 6 model name file but no mfsim.nam was found next to it")
        return [{"kind": "classic", "path": p, "version": classic_version(p)}]
    items, mf6_dirs = [], set()
    for f in sorted(p.rglob("mfsim.nam")):
        items.append({"kind": "mf6", "path": f})
        mf6_dirs.add(f.parent.resolve())
    for f in sorted(p.rglob("*.nam")):
        if f.name.lower() == "mfsim.nam" or f.parent.resolve() in mf6_dirs:
            continue
        e = _nam_entries(f)
        if not e or e[0][0] == "BEGIN":
            continue
        items.append({"kind": "classic", "path": f, "version": classic_version(f)})
    if not items:
        sys.exit(f"No mfsim.nam or MODFLOW name file found under {p}")
    return items


def _label(item):
    return item["path"].parent.name if item["kind"] == "mf6" else item["path"].stem


# --------------------------------------------------------------------------- #
# Grid & time (both families, via modelgrid / modeltime)
# --------------------------------------------------------------------------- #
def grid_info(model, header=None) -> dict:
    out = {}
    try:
        mg = model.modelgrid
    except Exception as e:
        return {"error": f"modelgrid not available: {e}"}
    out["type"] = type(mg).__name__.replace("Grid", "") + " grid"
    out["nlay"] = mg.nlay
    if getattr(mg, "grid_type", "") == "structured":
        out["rows x columns"] = f"{mg.nrow} x {mg.ncol}"
        try:
            out["delr (column width)"] = stats_text(arr_stats(mg.delr))
            out["delc (row height)"] = stats_text(arr_stats(mg.delc))
        except Exception:
            pass
        if mg.nrow == 1 or mg.ncol == 1:
            out["layout"] = "single row/column - cross-section model"
    else:
        out["cells per layer"] = getattr(mg, "ncpl", None)
    out["nodes"] = mg.nnodes
    try:
        x0, x1, y0, y1 = mg.extent
        out["extent x"] = f"{g(x0, 10)} .. {g(x1, 10)} ({g(x1 - x0)})"
        out["extent y"] = f"{g(y0, 10)} .. {g(y1, 10)} ({g(y1 - y0)})"
    except Exception:
        pass
    out["offset / rotation"] = f"xoff {g(mg.xoffset)}, yoff {g(mg.yoffset)}, angrot {g(mg.angrot)}"
    crs = None
    try:
        crs = mg.crs.to_string() if getattr(mg, "crs", None) is not None else None
    except Exception:
        crs = None
    if not crs and getattr(mg, "epsg", None):
        crs = f"EPSG:{mg.epsg}"
    out["crs"] = crs or "not set"
    out["length units"] = getattr(mg, "units", None)
    if header:
        h = {k: v for k, v in header.items() if v not in (None, "")}
        if h:
            out["name-file header"] = ", ".join(f"{k}={v}" for k, v in h.items())
            if crs is None and h.get("crs"):
                out["crs"] = f"{h['crs']} (from name-file header only; modelgrid has none)"
    try:
        top, botm = np.asarray(mg.top, dtype=float), np.asarray(mg.botm, dtype=float)
        act = None
        if getattr(mg, "idomain", None) is not None:
            act = np.asarray(mg.idomain) > 0 if np.asarray(mg.idomain).shape == botm.shape else None
        out["top"] = stats_text(arr_stats(top))
        out["bottom (all layers)"] = stats_text(arr_stats(botm))
        thick = mg.cell_thickness if hasattr(mg, "cell_thickness") else None
        if thick is not None:
            thick = np.asarray(thick, dtype=float)
            if act is not None and act.shape == thick.shape:
                tv = thick[act]
            else:
                tv = thick
            s = arr_stats(tv)
            out["cell thickness (active)"] = stats_text(s)
            if s and s.get("n") and s["min"] <= 0:
                out["warning"] = "zero or negative cell thickness"
        if act is not None:
            per = [int(act[k].sum()) for k in range(act.shape[0])]
            out["active cells"] = f"{int(act.sum())} of {act.size}" + (
                f"  (per layer: {', '.join(map(str, per))})" if 1 < len(per) <= 12 else "")
    except Exception:
        pass
    return {k: v for k, v in out.items() if v not in (None, "")}


def time_info(model, start=None) -> dict:
    out = {}
    try:
        mt = model.modeltime
    except Exception as e:
        return {"error": str(e)}
    perlen = np.asarray(mt.perlen, dtype=float).ravel()
    nstp = np.asarray(mt.nstp).ravel()
    tsm = np.asarray(mt.tsmult, dtype=float).ravel()
    ss = np.asarray(mt.steady_state).ravel() if getattr(mt, "steady_state", None) is not None else None
    units = str(getattr(mt, "time_units", "") or "")
    out["time units"] = units or "not set"
    out["start date/time"] = start or "not set"
    out["stress periods"] = len(perlen)
    out["total time"] = fmt_time(perlen.sum(), units)
    out["time steps"] = int(nstp.sum())
    if ss is not None and ss.size:
        nss = int(ss.astype(bool).sum())
        out["steady / transient periods"] = f"{nss} / {len(ss) - nss}"
    rows = []
    for i in range(len(perlen)):
        rows.append({"period": i + 1, "perlen": float(perlen[i]), "nstp": int(nstp[i]), "tsmult": float(tsm[i]) if i < len(tsm) else None,
                     "state": ("steady" if ss[i] else "transient") if ss is not None and i < len(ss) else ""})
    out["_periods"] = rows
    return {k: v for k, v in out.items() if v not in (None, "")}


def period_rows_text(rows, units, limit=12) -> list:
    out = [[r["period"], g(r["perlen"]), r["nstp"], g(r["tsmult"]), r["state"]] for r in rows[:limit]]
    if len(rows) > limit:
        out.append(["...", f"{len(rows) - limit} more", "", "", ""])
    return out


# --------------------------------------------------------------------------- #
# MODFLOW 6: generic package summary from blocks / datasets
# --------------------------------------------------------------------------- #
def mf6_dataset_summary(ds):
    t = type(ds).__name__
    try:
        if not ds.has_data():
            return None
    except Exception:
        pass
    try:
        d = ds.get_data()
    except Exception as e:
        return t, f"(could not read: {str(e)[:60]})"
    if d is None:
        return None
    try:
        if isinstance(d, dict):
            return t, period_summary(d)
        if isinstance(d, np.ndarray) and d.dtype.names:
            if len(d) <= 3 and len(d.dtype.names) <= 4:
                return t, "; ".join(", ".join(str(x) for x in r) for r in d)
            return t, rec_summary(d)
        if isinstance(d, np.ndarray):
            if d.size == 1:
                return t, g(d.ravel()[0])
            s = arr_stats(d)
            txt = stats_text(s)
            if d.ndim == 3 and 1 < d.shape[0] <= 10 and s and not s.get("constant"):
                per = []
                for k in range(d.shape[0]):
                    sk = arr_stats(d[k])
                    per.append(f"L{k + 1} " + (g(sk["min"]) if sk and sk.get("constant") else
                                               f"{g(sk['min'])}..{g(sk['max'])}" if sk and sk.get("n") else "-"))
                txt += "; " + ", ".join(per)
            return t, txt
        if isinstance(d, (bool, np.bool_)):
            return t, "on" if d else "off"
        if isinstance(d, (list, tuple)):
            return t, ", ".join(map(str, d))[:200]
        return t, str(d)[:200]
    except Exception as e:
        return t, f"(summary failed: {str(e)[:60]})"


def mf6_package(p, full=True) -> dict:
    ptype = p.package_type
    rows = []
    for bname, blk in p.blocks.items():
        for dname, ds in blk.datasets.items():
            r = mf6_dataset_summary(ds)
            if r:
                rows.append({"block": bname, "item": dname, "kind": r[0], "value": r[1]})
    fn = p.filename
    name = p.name[0] if isinstance(p.name, (list, tuple)) else p.name
    out = {"type": ptype, "name": name, "file": fn, "description": pkg_label(ptype), "datasets": rows,
           "is_stress": ptype in STRESS_TYPES}
    if not rows:
        out["note"] = "no data set - all defaults"
    # one-line summary
    keys = [f"{r['item']}={r['value']}" for r in rows if r["block"] in ("options", "dimensions", "nonlinear", "linear")
            and len(r["value"]) < 40][:6]
    per = [r for r in rows if r["block"] == "period"]
    if per:
        keys.append(per[0]["value"][:110])
    elif out["type"] in ("npf", "sto", "ic", "dis", "disv", "disu", "dsp", "adv", "mst", "est", "cnd"):
        keys += [f"{r['item']} {r['value']}" for r in rows if r["block"] == "griddata"][:3]
    out["summary"] = "; ".join(keys)
    return out


def mf6_listfile_status(path: Path) -> dict | None:
    if not path or not path.exists():
        return None
    try:
        txt = path.read_text(errors="replace")
    except OSError:
        return None
    out = {"file": path.name, "size_kB": round(path.stat().st_size / 1024, 1)}
    out["normal termination"] = "yes" if re.search(r"normal termination", txt, re.I) else "NO"
    nf = len(re.findall(r"failed to meet solver convergence|did not converge|failure to converge", txt, re.I))
    out["convergence failures"] = nf
    m = re.search(r"Elapsed run time:\s*(.+)", txt)
    if m:
        out["elapsed run time"] = m.group(1).strip()
    errs = [ln.strip() for ln in txt.splitlines() if re.match(r"\s*(ERROR|\*\*\*\* ERROR)", ln)]
    if errs:
        out["errors"] = errs[:10]
    warns = [ln.strip() for ln in txt.splitlines() if re.match(r"\s*WARNING", ln)]
    if warns:
        out["warnings"] = list(dict.fromkeys(warns))[:15]
    m = re.search(r"MODFLOW 6.*?\n.*?VERSION\s+([^\n]+)", txt, re.S)
    if m:
        out["code version"] = m.group(1).strip()
    return out


def budget_from_df(cum) -> dict | None:
    try:
        last = cum.iloc[-1]
    except Exception:
        return None
    terms = {}
    for k, v in last.items():
        terms[str(k)] = float(v)
    out = {"terms": terms}
    pd_ = terms.get("PERCENT_DISCREPANCY")
    if pd_ is not None:
        out["percent discrepancy"] = pd_
    out["time"] = str(cum.index[-1])
    return out


def head_summary(hf, inactive=()) -> dict:
    out = {}
    try:
        times = hf.get_times()
        out["times saved"] = len(times)
        if times:
            out["time span"] = f"{g(times[0])} .. {g(times[-1])}"
            h = np.asarray(hf.get_data(totim=times[-1]), dtype=float)
            m = np.isfinite(h) & (np.abs(h) < 1e20)
            for v in inactive:
                if v is not None:
                    m &= ~np.isclose(h, v)
            s = arr_stats(h[m])
            if s and s.get("n"):
                out["head at last time"] = f"{g(s['min'])} .. {g(s['max'])} (mean {g(s['mean'])}, {s['n']} active cells)"
            dry = int(np.isclose(h, -1e30).sum())
            if dry:
                out["dry cells (last time)"] = dry
            out["_last"] = h
            out["_last_time"] = times[-1]
    except Exception as e:
        out["error"] = str(e)[:120]
    return out


def cbc_summary(cf) -> dict:
    out = {}
    try:
        recs = [r.decode().strip() if isinstance(r, bytes) else str(r).strip() for r in cf.get_unique_record_names()]
        out["records"] = recs
        out["times saved"] = len(cf.get_times())
    except Exception as e:
        out["error"] = str(e)[:120]
    return out


def run_check(obj) -> dict:
    """FloPy check, run with the model workspace pointed at a temporary folder:
    some checks (e.g. SFR reach connections) write .chk/.csv files into model_ws."""
    import tempfile
    out = {"errors": [], "warnings": [], "passed": 0}
    old_ws = getattr(obj, "model_ws", None)
    tmp = tempfile.mkdtemp(prefix="flo_check_")
    try:
        with quiet() as buf:
            if old_ws is not None and hasattr(obj, "change_model_ws") and not hasattr(obj, "simulation"):
                obj.change_model_ws(tmp, reset_external=False)
                res = obj.check(f=str(Path(tmp) / "check.chk"), verbose=True, level=1)
            else:
                res = obj.check(verbose=True, level=1)
        txt = buf.getvalue()
    except Exception as e:
        return {"errors": [], "warnings": [], "failed": f"check could not run: {str(e)[:120]}"}
    finally:
        try:
            if old_ws is not None and hasattr(obj, "change_model_ws") and not hasattr(obj, "simulation"):
                with quiet():
                    obj.change_model_ws(old_ws, reset_external=False)
        except Exception:
            pass
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)
    chks = res if isinstance(res, list) else [res]
    for c in chks:
        sa = getattr(c, "summary_array", None)
        if sa is None:
            continue
        for r in sa:
            typ = str(r["type"]).strip().lower()
            desc = str(r["desc"]).strip()
            pkg = str(r["package"]).strip()
            loc = ""
            try:
                if "k" in sa.dtype.names and (r["i"] or r["j"] or r["k"]):
                    loc = f" (k={r['k']}, i={r['i']}, j={r['j']})"
            except Exception:
                pass
            (out["errors"] if typ.startswith("error") else out["warnings"]).append(f"{pkg}: {desc}{loc}")
    out["passed"] = len(re.findall(r"^\s{4}\S", txt.split("Checks that passed:")[-1], re.M)) if "Checks that passed" in txt else 0
    # collapse repeated messages
    for k in ("errors", "warnings"):
        cnt = {}
        for m in out[k]:
            base = re.sub(r" \(k=.*\)$", "", m)
            cnt[base] = cnt.get(base, 0) + 1
        out[k] = [f"{b}" + (f"  [{n} cells]" if n > 1 else "") for b, n in cnt.items()]
    return out


def build_mf6(item, only=None, do_check=True) -> list[dict]:
    simpath = item["path"]
    with quiet() as buf:
        sim = flopy.mf6.MFSimulation.load(sim_ws=str(simpath.parent), sim_name=simpath.parent.name,
                                          verbosity_level=0, verify_data=False)
    log = buf.getvalue()
    simd = {"kind": "MODFLOW 6 simulation", "name": simpath.parent.name, "path": str(simpath),
            "folder": str(simpath.parent), "load_log": [ln for ln in log.splitlines() if re.search(r"error|fail|warn", ln, re.I)][:20],
            "packages": [], "models": [], "warnings": []}
    simd["name_file"] = mf6_package(sim.name_file)
    for p in sim.sim_package_list:
        simd["packages"].append(mf6_package(p))
    try:
        tdis = sim.tdis
        simd["tdis"] = {"time units": tdis.time_units.get_data(), "start_date_time": tdis.start_date_time.get_data(),
                        "nper": tdis.nper.get_data()}
    except Exception:
        simd["tdis"] = {}
    simd["exchanges"] = []
    try:
        for ex in getattr(sim, "exchange_files", {}).values() if isinstance(getattr(sim, "exchange_files", None), dict) else []:
            simd["exchanges"].append(mf6_package(ex))
    except Exception:
        pass
    names = sim.model_names
    for mname in names:
        g6 = sim.get_model(mname)
        if only and only.lower() not in (mname.lower(), simpath.parent.name.lower()):
            continue
        if item.get("only_model_nam") and Path(g6.namefile).name.lower() != item["only_model_nam"].lower():
            continue
        md = {"name": mname, "type": g6.model_type, "namefile": g6.namefile,
              "grid_type": str(g6.get_grid_type()).split(".")[-1]}
        md["name_file"] = mf6_package(g6.name_file)
        md["grid"] = grid_info(g6)
        md["time"] = time_info(g6, simd["tdis"].get("start_date_time"))
        md["packages"] = [mf6_package(p) for p in g6.packagelist]
        md["_model"] = g6
        # results
        res = {}
        lst = Path(simpath.parent) / f"{Path(g6.namefile).stem}.lst"
        ls = mf6_listfile_status(lst)
        if ls:
            ls.pop("normal termination", None)      # MF6 writes this only to mfsim.lst
            res["list file"] = ls
        try:
            lb = g6.output.list()
            with quiet():
                inc, cum = lb.get_dataframes(start_datetime=None)
            b = budget_from_df(cum)
            if b:
                res["budget"] = b
        except Exception:
            pass
        for meth, kind in (("head", "head"), ("concentration", "concentration"), ("temperature", "temperature")):
            if hasattr(g6.output, meth):
                try:
                    hf = getattr(g6.output, meth)()
                    if hf is not None:
                        hs = head_summary(hf)
                        hs["file"] = Path(hf.filename).name
                        hs["quantity"] = kind
                        res[kind] = hs
                except Exception:
                    pass
        try:
            cf = g6.output.budget()
            if cf is not None:
                cs = cbc_summary(cf)
                cs["file"] = Path(cf.filename).name
                res["cell budget"] = cs
        except Exception:
            pass
        md["results"] = res
        md["check"] = run_check(g6) if do_check else None
        simd["models"].append(md)
    simd["sim_list"] = mf6_listfile_status(simpath.parent / "mfsim.lst")
    if do_check:
        try:
            with quiet() as buf:
                sim.check(verbose=True, level=0)
            simd["sim_check"] = [ln.strip() for ln in buf.getvalue().splitlines()
                                 if re.search(r"missing|not found|error", ln, re.I)
                                 and not re.match(r"\s*(Checking|No errors)", ln)]
        except Exception:
            simd["sim_check"] = []
    return [simd]


# --------------------------------------------------------------------------- #
# Classic MODFLOW (2005 / NWT / USG / 2000), MT3D, SEAWAT
# --------------------------------------------------------------------------- #
SKIP_ATTRS = {"parent", "name", "file_name", "unit_number", "extension", "url", "allowDuplicates", "acceptable_dtypes",
              "heading", "fn_path", "filenames", "extra", "dtype", "structured", "_generate_heading", "nrow_ncol_nlay_nper",
              "package_type", "_name", "_fn_path", "options"}


def classic_attr_summary(k, v):
    if isinstance(v, (Util2d, Util3d)):
        try:
            a = v.array
            s = arr_stats(a)
            txt = stats_text(s)
            if isinstance(v, Util3d) and a.ndim == 3 and 1 < a.shape[0] <= 10 and s and not s.get("constant"):
                txt += "; " + ", ".join(
                    f"L{i + 1} " + (g(arr_stats(a[i])["min"]) if arr_stats(a[i])["constant"] else
                                    f"{g(arr_stats(a[i])['min'])}..{g(arr_stats(a[i])['max'])}") for i in range(a.shape[0]))
            return "array", txt
        except Exception as e:
            return "array", f"(could not read: {str(e)[:60]})"
    if isinstance(v, (Transient2d, Transient3d)):
        try:
            d = {}
            for kper in range(v._model.nper if hasattr(v, "_model") else 1):
                if kper in getattr(v, "transient_2ds", {}) or kper in getattr(v, "transient_3ds", {}):
                    d[kper] = v.array[kper]
            if not d:
                d = {0: v.array[0]}
            return "transient array", period_summary(d)
        except Exception:
            return "transient array", "(transient array)"
    if isinstance(v, MfList):
        try:
            return "stress-period list", period_summary(v.data)
        except Exception:
            return "stress-period list", "(list)"
    if isinstance(v, np.ndarray):
        if v.dtype.names:
            return "table", rec_summary(v)
        s = arr_stats(v)
        return "array", stats_text(s) if s and v.size > 1 else str(v.ravel()[:1].tolist())
    if isinstance(v, dict):
        if not v:
            return None
        first = next(iter(v.values()))
        if hasattr(first, "dtype") and getattr(first.dtype, "names", None):
            return "per-period table", period_summary(v)
        return "dict", str(v)[:150]
    if isinstance(v, (bool, np.bool_)):
        return "flag", "on" if v else "off"
    if isinstance(v, (int, float, str, np.integer, np.floating)):
        return "value", str(v) if isinstance(v, str) else g(v)
    if isinstance(v, (list, tuple)):
        if not v:
            return None
        return "list", ", ".join(map(str, v))[:150]
    if v is None:
        return None
    return None


def classic_package(p) -> dict:
    ptype = p.name[0].lower() if isinstance(p.name, (list, tuple)) else str(p.name).lower()
    rows = []
    for k, v in vars(p).items():
        if k.startswith("_") or k in SKIP_ATTRS:
            continue
        r = classic_attr_summary(k, v)
        if r:
            rows.append({"block": "", "item": k, "kind": r[0], "value": r[1]})
    opt = getattr(p, "options", None)
    if opt and not hasattr(opt, "__dict__"):
        rows.append({"block": "", "item": "options", "kind": "list", "value": ", ".join(map(str, opt))})
    fn = p.file_name[0] if isinstance(p.file_name, (list, tuple)) else p.file_name
    unit = p.unit_number[0] if isinstance(p.unit_number, (list, tuple)) else p.unit_number
    out = {"type": ptype, "name": ptype.upper(), "file": fn, "unit": unit, "description": pkg_label(ptype),
           "datasets": rows, "is_stress": ptype in STRESS_TYPES}
    scal = [f"{r['item']}={r['value']}" for r in rows if r["kind"] in ("value", "flag") and len(r["value"]) < 20][:7]
    per = [r for r in rows if r["kind"] in ("stress-period list", "transient array", "per-period table")]
    if per:
        scal.append(f"{per[0]['item']}: {per[0]['value'][:110]}")
    out["summary"] = "; ".join(scal)
    return out


def classic_listfile_status(path: Path) -> dict | None:
    if not path or not path.exists():
        return None
    txt = path.read_text(errors="replace")
    out = {"file": path.name, "size_kB": round(path.stat().st_size / 1024, 1),
           "normal termination": "yes" if re.search(r"normal termination", txt, re.I) else "NO",
           "convergence failures": len(re.findall(r"failed to meet solver convergence|failure to converge|"
                                                  r"did not converge", txt, re.I))}
    m = re.search(r"Elapsed run time:\s*(.+)", txt)
    if m:
        out["elapsed run time"] = m.group(1).strip()
    errs = [ln.strip() for ln in txt.splitlines() if "ERROR" in ln.upper() and "PERCENT" not in ln.upper()]
    if errs:
        out["errors"] = list(dict.fromkeys(errs))[:10]
    return out


def load_classic(item, version=None):
    path, ver = item["path"], version or item["version"]
    kw = dict(model_ws=str(path.parent), verbose=True)       # verbose so load failures are reported (captured)
    with quiet() as buf:
        if ver.startswith("mt3d"):
            m = flopy.mt3d.Mt3dms.load(path.name, version=ver, forgive=True, **kw)
        elif ver == "seawat":
            m = flopy.seawat.Seawat.load(path.name, **kw)
        elif ver == "mfusg":
            m = flopy.mfusg.MfUsg.load(path.name, check=False, forgive=True, **kw)
        else:
            m = flopy.modflow.Modflow.load(path.name, version=ver, check=False, forgive=True, **kw)
    return m, buf.getvalue(), ver


VERSION_NAMES = {"mf2005": "MODFLOW-2005", "mfnwt": "MODFLOW-NWT", "mfusg": "MODFLOW-USG", "mf2k": "MODFLOW-2000",
                 "mt3dms": "MT3DMS", "mt3d-usgs": "MT3D-USGS", "seawat": "SEAWAT"}
SOLVERS = {"pcg", "pcgn", "nwt", "gmg", "sip", "de4", "sms", "gcg"}


def build_classic(item, version=None, do_check=True) -> list[dict]:
    path = item["path"]
    m, log, ver = load_classic(item, version)
    header = {}
    try:
        from flopy.utils import mfreadnam
        header = mfreadnam.attribs_from_namfile_header(str(path)) or {}
    except Exception:
        pass
    entries = _nam_entries(path)
    md = {"name": path.stem, "type": VERSION_NAMES.get(ver, ver), "namefile": path.name, "version": ver,
          "version_source": "forced with --version" if version else "detected from the name file",
          "grid_type": "structured" if getattr(m, "structured", True) else "unstructured"}
    md["name_file_entries"] = [{"ftype": e[0], "unit": e[1], "file": e[2], "exists": (path.parent / e[2]).exists()}
                               for e in entries]
    md["grid"] = grid_info(m, header)
    sdt = header.get("start_datetime") or getattr(getattr(m, "dis", None), "start_datetime", None)
    md["time"] = time_info(m, sdt)
    md["packages"] = [classic_package(p) for p in m.packagelist]
    md["_model"] = m
    loaded = {p["type"] for p in md["packages"]}
    issues, failed, lines = [], [], log.splitlines()
    for i, ln in enumerate(lines):
        mm = re.match(r"\s*(\S+)\s+package load\.\.\.failed", ln)
        if mm:
            why = lines[i + 1].strip() if i + 1 < len(lines) else ""
            failed.append(f"{mm.group(1)} ({why})" if why else mm.group(1))
            issues.append(f"{mm.group(1)} package failed to load: {why}")
        elif re.search(r"error|not loaded", ln, re.I):
            issues.append(ln.strip())
    md["load_log"] = issues[:20]
    md["failed_packages"] = failed
    md["header"] = header
    # name-file entries not loaded / missing
    notes = []
    for e in md["name_file_entries"]:
        if not e["exists"] and not e["ftype"].startswith("DATA") and e["ftype"] != "LIST":
            notes.append(f"{e['ftype']} file '{e['file']}' listed in the name file but missing")
    if ver not in ("mt3dms", "mt3d-usgs", "seawat") and not (loaded & SOLVERS):
        stray = sorted({f.suffix[1:].upper() for f in path.parent.glob(f"{path.stem}.*")
                        if f.suffix[1:].lower() in SOLVERS})
        notes.append("no solver package in the name file" + (f" ({', '.join(stray)} file exists in the folder but is not listed)" if stray else ""))
    if failed:
        notes.insert(0, f"FloPy could not load {len(failed)} package(s) listed in the name file: {'; '.join(failed)}"
                        + (" - wrong MODFLOW version? try --version mfnwt" if ver == "mf2005" else ""))
    md["notes"] = notes
    # results
    res = {}
    lst = next((path.parent / e[2] for e in entries if e[0] == "LIST"), None)
    ls = classic_listfile_status(lst) if lst else None
    if ls:
        res["list file"] = ls
        try:
            lbcls = flopy.utils.MfusgListBudget if ver == "mfusg" else (
                flopy.utils.MtListBudget if ver.startswith("mt3d") else flopy.utils.MfListBudget)
            with quiet():
                lb = lbcls(str(lst))
                if ver.startswith("mt3d"):
                    raise ValueError
                inc, cum = lb.get_dataframes(start_datetime=None)
            b = budget_from_df(cum)
            if b:
                res["budget"] = b
        except Exception:
            pass
    hnoflo = getattr(getattr(m, "bas6", None), "hnoflo", None)
    hdry = None
    for pk in ("lpf", "upw", "bcf6"):
        if getattr(m, pk, None) is not None:
            hdry = getattr(getattr(m, pk), "hdry", None)
    for e in entries:
        if not e[0].startswith("DATA"):
            continue
        f = path.parent / e[2]
        if not f.exists():
            continue
        try:
            with quiet():
                hf = flopy.utils.HeadFile(str(f))
            hs = head_summary(hf, (hnoflo, hdry))
            if hs.get("times saved"):
                hs["file"] = f.name
                hs["quantity"] = "head"
                res.setdefault("head", hs)
                continue
        except Exception:
            pass
        try:
            with quiet():
                cf = flopy.utils.CellBudgetFile(str(f))
            cs = cbc_summary(cf)
            if cs.get("records"):
                cs["file"] = f.name
                res.setdefault("cell budget", cs)
        except Exception:
            pass
    for ext in ("ucn",):
        for f in sorted(path.parent.glob(f"*.{ext}"))[:1]:
            try:
                with quiet():
                    uf = flopy.utils.UcnFile(str(f))
                hs = head_summary(uf)
                hs["file"] = f.name
                hs["quantity"] = "concentration"
                res["concentration"] = hs
            except Exception:
                pass
    md["results"] = res
    md["check"] = run_check(m) if do_check and hasattr(m, "check") else None
    return [{"kind": VERSION_NAMES.get(ver, ver) + " model", "name": path.stem, "path": str(path),
             "folder": str(path.parent), "packages": [], "models": [md], "warnings": [], "load_log": md["load_log"]}]


# --------------------------------------------------------------------------- #
# Warnings & summary
# --------------------------------------------------------------------------- #
def finalize(entry):
    w = entry["warnings"]
    for m in entry["models"]:
        tag = f"{m['name']}: " if len(entry["models"]) > 1 else ""
        for n in m.get("notes", []):
            w.append(tag + n)
        if m["grid"].get("warning"):
            w.append(tag + m["grid"]["warning"])
        if m["grid"].get("crs") == "not set":
            pass
        r = m["results"]
        b = r.get("budget")
        if b and b.get("percent discrepancy") is not None and abs(b["percent discrepancy"]) > 1:
            w.append(tag + f"water budget discrepancy {g(b['percent discrepancy'])} % (> 1 %) - check solver closure / convergence")
        ls = r.get("list file")
        if ls:
            if ls.get("normal termination", "yes") != "yes":
                w.append(tag + f"{ls['file']}: no 'normal termination' found")
            if ls["convergence failures"]:
                w.append(tag + f"{ls['file']}: {ls['convergence failures']} convergence failure message(s)")
        c = m.get("check")
        if c and c.get("errors"):
            w.append(tag + f"FloPy check: {len(c['errors'])} error(s) - " + "; ".join(c["errors"][:3]))
        m["summary"] = model_summary(entry, m)
    sl = entry.get("sim_list")
    if sl and sl["normal termination"] != "yes":
        w.append("mfsim.lst: no 'normal termination' found")
    for x in entry.get("sim_check", []) or []:
        w.append("simulation check: " + x)


def model_summary(entry, m) -> dict:
    gr, t, r = m["grid"], m["time"], m["results"]
    pk = m["packages"]
    out = {
        "Model type": m["type"] + (f" ({m['version_source']})" if m.get("version_source") else ""),
        "Grid": ", ".join(x for x in (gr.get("type"), f"{gr.get('nlay')} layer(s)",
                                      gr.get("rows x columns") and f"{gr['rows x columns']} rows x cols",
                                      gr.get("cells per layer") and f"{gr['cells per layer']} cells/layer") if x),
        "Active cells": gr.get("active cells"),
        "Extent": f"{gr.get('extent x', '')}; {gr.get('extent y', '')}" if gr.get("extent x") else None,
        "Length / time units": f"{gr.get('length units') or '?'} / {t.get('time units') or '?'}",
        "CRS": gr.get("crs"),
        "Stress periods": f"{t.get('stress periods')} ({t.get('steady / transient periods', '-')} steady/transient), "
                          f"{t.get('time steps')} time steps, total {t.get('total time')}",
        "Start date": t.get("start date/time"),
        "Flow / property packages": ", ".join(p["type"].upper() for p in pk if not p["is_stress"]
                                              and p["type"] not in ("oc", "list") and not p["type"].startswith("obs")),
        "Boundary / stress packages": ", ".join(pkg_tag(p) for p in pk if p["is_stress"]) or None,
        "Solver": ", ".join(p["type"].upper() for p in entry["packages"] if p["type"] in ("ims", "ems", "pts"))
                  or ", ".join(p["type"].upper() for p in pk if p["type"] in SOLVERS) or "none",
    }
    ls = entry.get("sim_list") or r.get("list file")
    if ls and r.get("list file"):
        out["Run"] = f"normal termination: {ls.get('normal termination', '?')}, convergence failures: {r['list file']['convergence failures']}" + (
            f", {ls['elapsed run time']}" if ls.get("elapsed run time") else "")
    else:
        out["Run"] = "no listing file - not run here"
    if r.get("budget") and r["budget"].get("percent discrepancy") is not None:
        out["Budget discrepancy"] = f"{g(r['budget']['percent discrepancy'])} %"
    h = r.get("head")
    if h and h.get("head at last time"):
        out["Head (last time)"] = h["head at last time"]
    if m.get("check"):
        c = m["check"]
        out["FloPy check"] = f"{len(c['errors'])} error(s), {len(c['warnings'])} warning(s), {c.get('passed', 0)} checks passed"
    return {k: v for k, v in out.items() if v not in (None, "")}


def build_config(arg, only=None, version=None, do_check=True) -> dict:
    items = discover(arg)
    if only:
        sel = [i for i in items if only.lower() in (_label(i).lower(), i["path"].stem.lower())]
        if not sel:  # maybe a model name inside an MF6 simulation
            sel = [i for i in items if i["kind"] == "mf6"]
        items = sel
    entries = []
    for it in items:
        try:
            ents = build_mf6(it, only if it["kind"] == "mf6" and only and only.lower() != _label(it).lower() else None,
                             do_check) if it["kind"] == "mf6" else build_classic(it, version, do_check)
        except Exception as e:
            ents = [{"kind": "load failed", "name": _label(it), "path": str(it["path"]), "folder": str(it["path"].parent),
                     "error": f"{type(e).__name__}: {e}", "models": [], "packages": [], "warnings": [f"FloPy could not load it: {e}"]}]
        for e in ents:
            if e["models"]:
                finalize(e)
            entries.append(e)
    entries = [e for e in entries if e["models"] or e.get("error")]
    if not entries:
        sys.exit("Nothing loaded (check --model).")
    root = Path(arg).expanduser().resolve()
    return {"root": str(root if root.is_dir() else root.parent), "flopy_version": flopy.__version__,
            "name": root.name if root.is_dir() else root.parent.name, "entries": entries}


# --------------------------------------------------------------------------- #
# HTML template & helpers (same layout as ras_config.py / hms_config.py)
# --------------------------------------------------------------------------- #
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<script>(function(){var t;try{t=localStorage.getItem('modelvue-theme')}catch(e){}if(t!=='dark')t='light';document.documentElement.dataset.theme=t;})();</script>
<title>MODFLOW · @@TITLE@@</title>
<style>
:root{--bg:#f6f7f9;--panel:#fff;--ink:#1d2433;--muted:#667085;--line:#e4e7ec;--acc:#2f6b2f;--accs:#e8f2e3;--warn:#b42318}
:root[data-theme=dark]{color-scheme:dark;--bg:#0f1318;--panel:#171c23;--ink:#e6e9ef;--muted:#98a2b3;--line:#2a313b;--acc:#9fd49a;--accs:#1b2e1a;--warn:#f97066}
*{box-sizing:border-box} body{margin:0;font:14px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,"Noto Sans TC",sans-serif;background:var(--bg);color:var(--ink)}
header{position:sticky;top:0;z-index:5;background:var(--acc);color:#fff;padding:14px 24px;display:flex;gap:16px;align-items:baseline;flex-wrap:wrap}
header h1{margin:0;font-size:20px} header span{opacity:.85}
[data-theme=dark] header{background:#1f3d1f}
header .theme{margin-left:auto;font:inherit;font-size:12px;padding:3px 12px;border:1px solid rgba(255,255,255,.55);border-radius:14px;background:transparent;color:#fff;cursor:pointer} header .theme:hover{background:rgba(255,255,255,.15)}

.wrap{display:grid;grid-template-columns:240px 1fr;gap:0}
nav{position:sticky;top:56px;align-self:start;height:calc(100vh - 56px);overflow:auto;padding:16px 10px;border-right:1px solid var(--line)}
nav a{display:block;padding:5px 10px;border-radius:6px;color:var(--ink);text-decoration:none;font-size:13px}
nav a.sub{padding-left:22px;font-size:12.5px;color:var(--muted)}
nav a:hover,nav a.on{background:var(--accs);color:var(--acc)}
nav .grp{margin:14px 10px 4px;font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}
.tools{display:flex;gap:6px;padding:0 6px 10px;border-bottom:1px solid var(--line);margin-bottom:6px}
.tools button{flex:1;font:inherit;font-size:12px;padding:4px 6px;border:1px solid var(--line);border-radius:6px;background:var(--panel);color:var(--ink);cursor:pointer}
.tools button:hover{background:var(--accs);color:var(--acc)}
main{padding:20px 28px 60px;min-width:0} main a{color:var(--acc)}
section{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:18px 22px;margin-bottom:22px;scroll-margin-top:70px}
h2{margin:0 0 10px;font-size:19px;color:var(--acc);display:flex;gap:10px;align-items:center;flex-wrap:wrap}
h4{margin:16px 0 6px;font-size:14px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
.muted{color:var(--muted)} .warn{color:var(--warn)} code{font-size:12px;background:var(--accs);padding:1px 5px;border-radius:4px}
.badge{font-size:11px;background:var(--acc);color:#fff;border-radius:10px;padding:0 8px;font-weight:500}
.badge.o{background:var(--accs);color:var(--acc)} .badge.w{background:var(--warn)}
.chip{display:inline-block;background:var(--accs);color:var(--acc);border-radius:4px;padding:0 6px;margin:1px 3px 1px 0;font-size:12px}
.tw{overflow-x:auto} table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
td{white-space:nowrap} td.wrap,td .wrap{white-space:normal}
thead th{background:var(--accs);color:var(--acc);cursor:pointer;user-select:none;white-space:nowrap}
thead th:hover{text-decoration:underline} tbody tr:hover{background:var(--accs)}
table.kv{width:auto;min-width:50%} table.kv th{color:var(--muted);font-weight:500;background:none;cursor:auto;padding-right:24px;white-space:nowrap}
table.kv td{white-space:normal;overflow-wrap:anywhere}
td.msg{white-space:normal;width:88%}
td.desc{white-space:normal;color:var(--muted);font-size:12px;min-width:260px}
.stats{display:flex;flex-wrap:wrap;gap:10px;margin-bottom:12px} .stat{border:1px solid var(--line);border-radius:8px;padding:8px 14px;min-width:100px}
.stat b{display:block;font-size:22px;color:var(--acc)} .stat span{color:var(--muted);font-size:12px}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:0 24px;align-items:start}
details.blk{border:1px solid var(--line);border-radius:8px;margin:12px 0;background:var(--panel);scroll-margin-top:70px}
details.blk>summary{cursor:pointer;padding:8px 12px;font-weight:600;color:var(--acc);list-style:none}
details.blk>summary::-webkit-details-marker{display:none}
details.blk>summary:before{content:"\\25B8";display:inline-block;width:16px;color:var(--muted)}
details.blk[open]>summary:before{content:"\\25BE"}
details.blk>summary small{color:var(--muted);font-weight:400;margin-left:6px}
details.blk>.cb{padding:0 12px 12px}
svg.plot{width:100%;max-width:980px;height:auto;display:block;margin:8px 0;font-size:11px}
svg.plot text{fill:var(--muted)} svg.plot .lab{fill:var(--ink)} svg.plot .gl{stroke:var(--line)} svg.plot .frame{fill:none;stroke:var(--line)}
svg.plan{width:100%;max-width:980px;height:auto;display:block;margin:8px 0;background:var(--bg);border:1px solid var(--line);border-radius:6px}
svg.plan .cell{stroke:rgba(0,0,0,.18);stroke-width:.4} svg.plan .cell:hover{stroke:#000;stroke-width:1.2}
svg.plan .chan{fill:none;stroke:#1565c0;stroke-width:3} svg.plan .pipe{fill:none;stroke:#e65100;stroke-width:2.4;stroke-dasharray:6 3}
svg.plan .lnk{stroke:#7b1fa2;stroke-width:1;opacity:.8} svg.plan .bnd{fill:none;stroke:#d32f2f;stroke-width:4}
svg.plan .m1{fill:#0d47a1} svg.plan .stor{fill:#e65100;stroke:#fff} svg.plan .nn{fill:#1565c0;stroke:#fff}
.legend{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:12px;color:var(--muted);margin-bottom:8px}
.legend i{display:inline-block;vertical-align:middle;margin-right:5px}
.legend .sw{width:90px;height:10px;border-radius:2px} .legend .ln{width:22px;height:0;border-top:3px solid}
.legend .ln.chan{border-color:#1565c0} .legend .ln.pipe{border-color:#e65100;border-top-style:dashed}
.legend .ln.lnk{border-color:#7b1fa2;border-top-width:1px} .legend .ln.bnd{border-color:#d32f2f}
.legend .bx{width:9px;height:9px} .legend .bx.stor{background:#e65100} .legend .bx.nn{background:#1565c0}
ul.warns{margin:6px 0;padding-left:20px;color:var(--warn)}
@media (max-width:1000px){.wrap{grid-template-columns:1fr} nav{display:none}}
footer{color:var(--muted);font-size:12px;text-align:center;padding:10px}
@media print{nav,header{display:none} .wrap{display:block} details>*{display:block}}
</style></head>
<body>
<header><h1>MODFLOW · @@TITLE@@</h1><span>@@SUB@@</span><button class="theme" type="button" title="Switch between the light and dark theme">☾ Dark</button></header>
<div class="wrap"><nav><div class="tools"><button data-x="1">Expand all</button><button data-x="0">Collapse all</button></div>@@NAV@@</nav>
<main>@@BODY@@
<footer>Generated by modelvue.flovue on @@GEN@@ from <code>@@PATH@@</code></footer></main></div>
<script>
document.querySelectorAll('.tools button').forEach(b=>b.addEventListener('click',()=>{
  const o=b.dataset.x==='1'; document.querySelectorAll('main details').forEach(d=>d.open=o);}));
document.querySelectorAll('nav a').forEach(a=>a.addEventListener('click',()=>{
  const t=document.querySelector(a.getAttribute('href')); if(t&&t.tagName==='DETAILS')t.open=true;}));
document.querySelectorAll('table.sortable thead th').forEach((th,i)=>th.addEventListener('click',()=>{
  const tb=th.closest('table').tBodies[0], rows=[...tb.rows], asc=th.dataset.asc!=='1';
  rows.sort((a,b)=>{const x=a.cells[i].innerText,y=b.cells[i].innerText,nx=parseFloat(x),ny=parseFloat(y);
    const r=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y); return asc?r:-r;});
  th.dataset.asc=asc?'1':'0'; rows.forEach(r=>tb.appendChild(r));}));
// MDU parameter filter
document.querySelectorAll('input.pfilter').forEach(inp=>inp.addEventListener('input',()=>{
  const q=inp.value.toLowerCase(), box=inp.closest('.cb');
  box.querySelectorAll('tbody tr').forEach(r=>r.style.display=r.innerText.toLowerCase().includes(q)?'':'none');
  if(q) box.querySelectorAll('details').forEach(d=>d.open=true);}));
const links=[...document.querySelectorAll('nav a')];
const obs=new IntersectionObserver(es=>es.forEach(e=>{if(e.isIntersecting){links.forEach(l=>l.classList.toggle('on',l.getAttribute('href')==='#'+e.target.id));}}),{rootMargin:'-30% 0px -60% 0px'});
document.querySelectorAll('main section, main details.blk[id]').forEach(s=>obs.observe(s));
const themeBtn=document.querySelector('header .theme');
const themeLabel=()=>{themeBtn.textContent=document.documentElement.dataset.theme==='dark'?'☀ Light':'☾ Dark';};
themeLabel();
themeBtn.addEventListener('click',()=>{const t=document.documentElement.dataset.theme==='dark'?'light':'dark';
  document.documentElement.dataset.theme=t; themeLabel(); try{localStorage.setItem('modelvue-theme',t)}catch(e){}});
</script></body></html>"""


# --------------------------------------------------------------------------- #
# Terminal output
# --------------------------------------------------------------------------- #


def text_table(headers, rows, indent=2, maxw=60):
    if not rows:
        return
    w = [min(maxw, max(len(str(h)), *(len(str(r[i])) for r in rows))) for i, h in enumerate(headers)]
    pad = " " * indent
    print(pad + "  ".join(str(h)[:w[i]].ljust(w[i]) for i, h in enumerate(headers)))
    print(pad + "  ".join("-" * x for x in w))
    for r in rows:
        print(pad + "  ".join(str("" if c is None else c)[:w[i]].ljust(w[i]) for i, c in enumerate(r)))


def show_packages(pkgs, full):
    for p in pkgs:
        lab = pkg_tag(p)
        print(f"  {lab:<10} {p['description'][:34]:<34} {p['file'] or ''}")
        if full:
            for r in p["datasets"]:
                blk = f"[{r['block']}] " if r["block"] else ""
                print(f"      {blk}{r['item']}: {r['value']}")
            if p.get("note"):
                print(f"      ({p['note']})")
        elif p.get("summary"):
            print(f"      {p['summary'][:150]}")
        elif p.get("note"):
            print(f"      ({p['note']})")


def show_model(entry, m, full):
    h2(f"Model {m['name']} ({m['type']})")
    for k, v in m["summary"].items():
        row(k, v)
    h2("Grid")
    for k, v in m["grid"].items():
        if k != "warning":
            row(k, v)
    h2("Time discretisation")
    t = m["time"]
    for k, v in t.items():
        if not k.startswith("_"):
            row(k, v)
    if t.get("_periods"):
        text_table(["Period", "Length", "Steps", "Mult", "State"], period_rows_text(t["_periods"], t.get("time units")), indent=4)
    if m.get("name_file_entries"):
        h2(f"Name file {m['namefile']}")
        text_table(["Ftype", "Unit", "File", "Exists"],
                   [[e["ftype"], e["unit"], e["file"], "yes" if e["exists"] else "NO"] for e in m["name_file_entries"]])
    h2(f"Packages ({len(m['packages'])})" + ("" if full else " - use --all for every dataset"))
    show_packages(m["packages"], full)
    c = m.get("check")
    if c:
        h2("FloPy model check")
        if c.get("failed"):
            print("  " + c["failed"])
        print(f"  {len(c['errors'])} error(s), {len(c['warnings'])} warning(s), {c.get('passed', 0)} checks passed")
        for x in c["errors"]:
            print(f"    ERROR   {x}")
        for x in c["warnings"]:
            print(f"    warning {x}")
    h2("Results")
    r = m["results"]
    if not r:
        print("  No listing / head / budget output found (model not run here).")
    for k, v in r.items():
        if k == "budget":
            print(f"  Water budget (cumulative, end of run{', time ' + v['time'] if v.get('time') else ''})")
            terms = v["terms"]
            ins = [(a, b) for a, b in terms.items() if a.endswith("_IN") and a != "TOTAL_IN" and b]
            outs = [(a, b) for a, b in terms.items() if a.endswith("_OUT") and a != "TOTAL_OUT" and b]
            for a, b in ins + outs:
                print(f"      {a:<28} {g(b, 6):>14}")
            for a in ("TOTAL_IN", "TOTAL_OUT", "IN-OUT", "PERCENT_DISCREPANCY"):
                if a in terms:
                    print(f"      {a:<28} {g(terms[a], 6):>14}")
            continue
        print(f"  {k}: {v.get('file', '')}")
        for kk, vv in v.items():
            if kk not in ("file",) and not kk.startswith("_"):
                row(kk, vv, indent=6)


def show(cfg, full=False):
    for e in cfg["entries"]:
        h1(f"{e['kind'].upper()}: {e['name']}")
        row("File", e["path"])
        if e.get("error"):
            print(f"  ERROR: {e['error']}")
            continue
        if e["kind"].startswith("MODFLOW 6"):
            row("Models", ", ".join(f"{m['name']} ({m['type']})" for m in e["models"]))
            td = e.get("tdis", {})
            row("TDIS", f"nper {td.get('nper')}, time units {td.get('time units') or 'not set'}, start {td.get('start_date_time') or 'not set'}")
            h2("Simulation packages")
            show_packages([e["name_file"]] + e["packages"] + e.get("exchanges", []), full or True)
            if e.get("sim_list"):
                row("mfsim.lst", {k: v for k, v in e["sim_list"].items() if k not in ("file", "warnings", "errors")})
        if e.get("load_log"):
            h2("Load messages")
            for ln in e["load_log"]:
                print("  " + ln)
        for m in e["models"]:
            show_model(e, m, full)
        if e["warnings"]:
            h2("Warnings")
            for w in e["warnings"]:
                print(f"  ! {w}")
    print(f"\n(read with FloPy {cfg['flopy_version']})\n")


# --------------------------------------------------------------------------- #
# Figures (FloPy plotting + matplotlib, optional)
# --------------------------------------------------------------------------- #
AREAL = {"ssm", "api", "uzf", "rch", "rcha", "evt", "evta", "ets", "fmi"}
BC_COLORS = ["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]


def _mpl():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except Exception:
        return None


def _png(fig, plt) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _stress_objs(m):
    mod = m["_model"]
    out = []
    for p in m["packages"]:
        if not p["is_stress"] or p["type"] in AREAL:      # areal packages would cover the whole map
            continue
        if hasattr(mod, "simulation"):  # MF6
            obj = mod.get_package(p["name"])
            out.append((p["name"].upper(), {"package": obj}))
        else:
            out.append((p["type"].upper(), {"ftype": p["type"].upper()}))
    return out


def _plot_bcs(view, m, ax):
    from matplotlib.patches import Patch
    handles = []
    for i, (lab, kw) in enumerate(_stress_objs(m)):
        col = BC_COLORS[i % len(BC_COLORS)]
        try:
            with quiet():
                if "package" in kw:
                    view.plot_bc(package=kw["package"], color=col)
                else:
                    view.plot_bc(kw["ftype"], color=col)
            handles.append(Patch(color=col, label=lab))
        except Exception:
            pass
    return handles


def fig_grid(m) -> list[tuple[str, str]]:
    plt = _mpl()
    if plt is None:
        return []
    mod = m["_model"]
    mg = mod.modelgrid
    figs = []
    xsec = getattr(mg, "grid_type", "") == "structured" and mg.nrow == 1
    try:
        if not xsec:
            fig, ax = plt.subplots(figsize=(8, 6.5))
            pmv = flopy.plot.PlotMapView(model=mod, ax=ax, layer=0)
            top = np.asarray(mg.top, dtype=float)
            idom = np.asarray(mg.idomain) if mg.idomain is not None else None
            if idom is not None and idom.ndim >= 2:
                top = np.ma.masked_where(idom[0] <= 0, top)
            pa = pmv.plot_array(top, cmap="terrain")
            fig.colorbar(pa, ax=ax, shrink=0.8, label="top elevation (layer 1)")
            try:
                pmv.plot_inactive(color_noflow="#d9d9d9")
            except Exception:
                pass
            if mg.nnodes / max(mg.nlay, 1) <= 6000:
                pmv.plot_grid(lw=0.3, color="0.35")
            hs = _plot_bcs(pmv, m, ax)
            if hs:
                ax.legend(handles=hs, loc="upper right", fontsize=8, framealpha=0.9)
            ax.set_title(f"{m['name']}: plan view, layer 1 - top, inactive cells and boundary packages", fontsize=10)
            ax.set_aspect("equal")
            figs.append(("Plan view", _png(fig, plt)))
        if getattr(mg, "grid_type", "") == "structured" and (mg.nlay > 1 or xsec):
            line = {"row": 0 if xsec else mg.nrow // 2}
            fig, ax = plt.subplots(figsize=(10, 4))
            pxs = flopy.plot.PlotCrossSection(model=mod, ax=ax, line=line)
            try:
                pxs.plot_inactive(color_noflow="#d9d9d9")
            except Exception:
                pass
            pxs.plot_grid(lw=0.25, color="0.4")
            hs = _plot_bcs(pxs, m, ax)
            if hs:
                ax.legend(handles=hs, loc="upper right", fontsize=8, framealpha=0.9)
            ax.set_title(f"Cross-section along row {line['row'] + 1}: layers, inactive cells and boundary packages", fontsize=10)
            figs.append(("Cross-section", _png(fig, plt)))
    except Exception as e:
        figs.append(("Figure failed", f"ERROR:{type(e).__name__}: {e}"))
    return figs


def fig_head(m) -> list[tuple[str, str]]:
    plt = _mpl()
    r = m["results"]
    h = next((r[k] for k in ("head", "concentration", "temperature") if k in r and r[k].get("_last") is not None), None)
    if plt is None or h is None:
        return []
    mod = m["_model"]
    mg = mod.modelgrid
    arr = np.ma.masked_where(~np.isfinite(h["_last"]) | (np.abs(h["_last"]) >= 1e20) | np.isclose(h["_last"], -999.99)
                             | np.isclose(h["_last"], -999.0), h["_last"])
    q = h.get("quantity", "head")
    figs = []
    try:
        xsec = getattr(mg, "grid_type", "") == "structured" and mg.nrow == 1
        if xsec:
            fig, ax = plt.subplots(figsize=(10, 4))
            v = flopy.plot.PlotCrossSection(model=mod, ax=ax, line={"row": 0})
            pa = v.plot_array(arr, cmap="viridis")
            try:
                cs = v.contour_array(arr, colors="white", linewidths=0.8)
                ax.clabel(cs, fmt="%.3g", fontsize=7)
            except Exception:
                pass
            v.plot_grid(lw=0.15, color="0.6")
        else:
            fig, ax = plt.subplots(figsize=(8, 6.5))
            v = flopy.plot.PlotMapView(model=mod, ax=ax, layer=0)
            pa = v.plot_array(arr, cmap="viridis")
            try:
                cs = v.contour_array(arr, colors="white", linewidths=0.8)
                ax.clabel(cs, fmt="%.3g", fontsize=7)
            except Exception:
                pass
            try:
                v.plot_inactive(color_noflow="#d9d9d9")
            except Exception:
                pass
            ax.set_aspect("equal")
        fig.colorbar(pa, ax=ax, shrink=0.8, label=q)
        ax.set_title(f"{q.capitalize()} at last saved time ({g(h.get('_last_time'))})" + ("" if xsec else ", layer 1"), fontsize=10)
        figs.append((q.capitalize(), _png(fig, plt)))
    except Exception as e:
        figs.append(("Figure failed", f"ERROR:{type(e).__name__}: {e}"))
    return figs


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #
EXTRA_CSS = """
figure{margin:10px 0} figure img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px;background:#fff}
figure figcaption{font-size:12px;color:var(--muted)}
.figs{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:12px}
td.msg{white-space:normal} .pk{font-weight:600}
"""


def _figs_html(figs):
    out = []
    for cap, data in figs:
        if data.startswith("ERROR:"):
            out.append(f'<p class="muted">{_e(cap)}: {_e(data[6:])}</p>')
        else:
            out.append(f'<figure><img alt="{_e(cap)}" src="data:image/png;base64,{data}"><figcaption>{_e(cap)}</figcaption></figure>')
    return f'<div class="figs">{"".join(out)}</div>' if out else ""


def _pkg_table(pkgs):
    rows = []
    for p in pkgs:
        lab = pkg_tag(p)
        body = _table(["Block", "Item", "Kind", "Value"],
                      [[_e(r["block"]), f"<b>{_e(r['item'])}</b>", _e(r["kind"]), (_e(r["value"]), "msg")]
                       for r in p["datasets"]], sortable=False) or f'<p class="muted">{_e(p.get("note", "no data"))}</p>'
        rows.append(_blk(f'<span class="pk">{_e(lab)}</span> <span class="muted" style="font-weight:400">'
                         f'{_e(p["description"])} · {_e(p["file"])}</span>',
                         body, (p.get("summary") or "")[:140], open_=False))
    return "".join(rows)


def _results_html(m):
    r = m["results"]
    if not r:
        return '<p class="muted">No listing, head or budget output found (model not run here).</p>'
    body = ""
    for k, v in r.items():
        if k == "budget":
            t = v["terms"]
            ins = [(a, b) for a, b in t.items() if a.endswith("_IN") and a != "TOTAL_IN"]
            outs = [(a, b) for a, b in t.items() if a.endswith("_OUT") and a != "TOTAL_OUT"]
            rows = [[_e(a), _e(g(b, 6)), ""] for a, b in ins] + [[_e(a), "", _e(g(b, 6))] for a, b in outs]
            rows = [r + [""] for r in rows]
            rows += [[f"<b>{_e(a)}</b>", _e(g(t[a], 6)) if a == "TOTAL_IN" else "",
                      _e(g(t[a], 6)) if a == "TOTAL_OUT" else "",
                      f"<b>{_e(g(t[a], 6))}</b>" if a in ("IN-OUT", "PERCENT_DISCREPANCY") else ""]
                     for a in ("TOTAL_IN", "TOTAL_OUT", "IN-OUT", "PERCENT_DISCREPANCY") if a in t]
            pdisc = v.get("percent discrepancy")
            flag = f' <span class="badge{" w" if pdisc is not None and abs(pdisc) > 1 else ""}">discrepancy {g(pdisc)} %</span>' if pdisc is not None else ""
            body += f"<h4>Water budget, cumulative at end of run{flag}</h4>" + _table(["Term", "In", "Out", "Balance"], rows, sortable=False)
            continue
        pairs = [(kk, vv) for kk, vv in v.items() if not kk.startswith("_") and kk not in ("errors", "warnings")]
        body += f"<h4>{_e(k)}</h4>" + _kv(pairs)
        for lvl in ("errors", "warnings"):
            if v.get(lvl):
                body += _table([lvl.capitalize()], [[(_e(x), "msg")] for x in v[lvl]], sortable=False)
    return body + _figs_html(fig_head(m))


def _model_html(entry, m, mi, ei):
    pid = lambda k: f"e{ei}m{mi}-{k}"
    present, out = [], []

    def add(key, title, body, note="", open_=True):
        if body:
            present.append((pid(key), title))
            out.append(_blk(_e(title), body, note, open_, pid(key)))

    out.append(f'<h3 style="margin:18px 0 4px;color:var(--acc)">Model {_e(m["name"])} '
               f'<span class="badge o">{_e(m["type"])}</span></h3>')
    s = list(m["summary"].items())
    add("sum", "Summary", _two(s))
    gr = {k: v for k, v in m["grid"].items() if k != "warning"}
    add("grid", "Grid", _two(list(gr.items())) + _figs_html(fig_grid(m)), m.get("grid_type", ""))
    t = m["time"]
    body = _two([(k, v) for k, v in t.items() if not k.startswith("_")])
    if t.get("_periods"):
        body += _table(["Period", "Length", "Steps", "Multiplier", "State"],
                       [[_e(x) for x in r] for r in period_rows_text(t["_periods"], t.get("time units"), limit=60)])
    add("time", "Time discretisation", body, f"{t.get('stress periods')} periods")
    if m.get("name_file_entries"):
        add("nam", f"Name file · {m['namefile']}", _table(["Ftype", "Unit", "File", "Exists"],
            [[_e(e["ftype"]), _e(e["unit"]), _e(e["file"]), "yes" if e["exists"] else '<span class="warn">NO</span>']
             for e in m["name_file_entries"]]), open_=False)
    elif m.get("name_file"):
        add("nam", f"Name file · {Path(m['namefile']).name}", _pkg_table([m["name_file"]]), open_=False)
    add("pkgs", f"Packages ({len(m['packages'])})", _pkg_table(m["packages"]), "click a package for every dataset")
    c = m.get("check")
    if c:
        body = (f'<p>{len(c["errors"])} error(s), {len(c["warnings"])} warning(s), {c.get("passed", 0)} checks passed</p>'
                + (f'<p class="warn">{_e(c["failed"])}</p>' if c.get("failed") else "")
                + _table(["Level", "Message"], [['<span class="warn">ERROR</span>', (_e(x), "msg")] for x in c["errors"]]
                         + [["warning", (_e(x), "msg")] for x in c["warnings"]], sortable=False))
        add("check", "FloPy model check", body, f"{len(c['errors'])} errors · {len(c['warnings'])} warnings")
    add("res", "Results", _results_html(m))
    if m.get("load_log"):
        add("log", "Load messages", _table(["Message"], [[(_e(x), "msg")] for x in m["load_log"]]), open_=False)
    return "".join(out), present


def render_html(cfg) -> str:
    nav, body = [], []
    for ei, e in enumerate(cfg["entries"]):
        sec = [f'<section id="e{ei}"><h2>{_e(e["name"])} <span class="badge o">{_e(e["kind"])}</span>']
        if e.get("error"):
            sec.append(f'</h2><p class="warn">{_e(e["error"])}</p></section>')
            body.append("".join(sec))
            nav.append(f'<a href="#e{ei}"><b>{_e(e["name"])}</b></a>')
            continue
        bad = [m for m in e["models"] if (e.get("sim_list") or m["results"].get("list file", {})).get("normal termination") == "yes"]
        if any(m["results"] for m in e["models"]):
            sec.append(f' <span class="badge{"" if len(bad) == len(e["models"]) else " w"}">'
                       f'{"run finished" if len(bad) == len(e["models"]) else "run problem"}</span>')
        else:
            sec.append(' <span class="badge o">not run</span>')
        if e["warnings"]:
            sec.append(f' <span class="badge w">{len(e["warnings"])} warning(s)</span>')
        sec.append(f'</h2><p class="muted"><code>{_e(e["path"])}</code></p>')
        if e["warnings"]:
            sec.append('<ul class="warns">' + "".join(f"<li>{_e(w)}</li>" for w in e["warnings"]) + "</ul>")
        nav.append(f'<a href="#e{ei}"><b>{_e(e["name"])}</b><br><small class="muted">{_e(e["kind"])}</small></a>')
        if e["kind"].startswith("MODFLOW 6"):
            td = e.get("tdis", {})
            sb = _kv([("models", ", ".join(f"{m['name']} ({m['type']})" for m in e["models"])),
                      ("time units", td.get("time units") or "not set"), ("start date/time", td.get("start_date_time") or "not set"),
                      ("stress periods", td.get("nper"))])
            if e.get("sim_list"):
                sb += "<h4>mfsim.lst</h4>" + _kv([(k, v) for k, v in e["sim_list"].items() if k not in ("warnings", "errors")])
            sb += _pkg_table([e["name_file"]] + e["packages"] + e.get("exchanges", []))
            sec.append(_blk("Simulation: name file, TDIS, solver", sb, "mfsim.nam", True, f"e{ei}-sim"))
            nav.append(f'<a class="sub" href="#e{ei}-sim">Simulation</a>')
        for mi, m in enumerate(e["models"]):
            html_, present = _model_html(e, m, mi, ei)
            sec.append(html_)
            if len(e["models"]) > 1:
                nav.append(f'<a class="sub" href="#{present[0][0]}"><i>{_e(m["name"])}</i></a>')
            nav += [f'<a class="sub" href="#{b}">{_e(t)}</a>' for b, t in present]
        sec.append("</section>")
        body.append("".join(sec))
    sub = f'{len(cfg["entries"])} model set(s) · FloPy {cfg["flopy_version"]}'
    rep = {"@@TITLE@@": _e(cfg["name"]), "@@SUB@@": _e(sub), "@@NAV@@": "".join(nav), "@@BODY@@": "".join(body),
           "@@GEN@@": _dt.datetime.now().strftime("%Y-%m-%d %H:%M"), "@@PATH@@": _e(cfg["root"])}
    out = HTML_TEMPLATE.replace("</style>", EXTRA_CSS + "</style>")
    for k, v in rep.items():
        out = out.replace(k, v)
    return out


def write_html(cfg, out: Path, open_browser=True):
    out.write_text(render_html(cfg), encoding="utf-8")
    print(f"Webpage written to {out.resolve()}")
    if open_browser:
        import webbrowser
        webbrowser.open(out.resolve().as_uri())


def strip_private(o):
    if isinstance(o, dict):
        return {k: strip_private(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, list):
        return [strip_private(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    return o


# --------------------------------------------------------------------------- #
def main(argv=None, prog=None):
    ap = argparse.ArgumentParser(prog=prog, description="Show the configuration of a MODFLOW model, read with FloPy.")
    ap.add_argument("model", help="mfsim.nam, a classic .nam file, or a folder containing models")
    ap.add_argument("--model", dest="only", metavar="NAME", help="only this simulation folder / model / name file")
    ap.add_argument("--version", choices=sorted(VERSION_NAMES), help="force the version of a classic name-file model")
    ap.add_argument("--all", action="store_true", help="terminal: show every dataset of every package")
    ap.add_argument("--no-check", action="store_true", help="skip FloPy's model check")
    ap.add_argument("--json", metavar="FILE", help="also write the configuration to a JSON file")
    ap.add_argument("--html", nargs="?", const="", metavar="FILE",
                    help="write the configuration as a webpage and open it in the browser, instead of printing it "
                         "(default file: <name>_config.html in the model folder)")
    ap.add_argument("--no-browser", action="store_true", help="with --html, write the page but don't open it")
    a = ap.parse_args(argv)

    cfg = build_config(a.model, a.only, a.version, not a.no_check)
    if a.html is None:
        show(cfg, a.all)
    else:
        e0 = cfg["entries"][0]
        folder = Path(e0["folder"]) if len(cfg["entries"]) == 1 else Path(cfg["root"])
        name = e0["name"] if len(cfg["entries"]) == 1 else cfg["name"]
        out = Path(a.html) if a.html else folder / f"{name}_config.html"
        write_html(cfg, out, open_browser=not a.no_browser)
    if a.json:
        Path(a.json).write_text(json.dumps(strip_private(cfg), indent=2, default=str), encoding="utf-8")
        print(f"Configuration written to {a.json}")


if __name__ == "__main__":
    main()
