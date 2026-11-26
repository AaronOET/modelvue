"""
modelvue.xbvue - Read an XBeach model setup and print its configuration.

Parses the plain-text XBeach inputs - params.txt, grid files (x.grd / y.grd),
bathymetry (bed.dep), wave boundary files (jonswap / jons_table / FILELIST),
water-level files (zs0file), vegetation files - and, for case sets, the
Deltares OET test-bed '.config' file that lists several runs.

Usage
-----
    xbvue  path/to/params.txt                 # one case
    xbvue  path/to/case_folder                # finds params.txt
    xbvue  path/to/case_set                   # .config or sub-folders with params.txt
    xbvue  Holland1976 --case Raai1085        # one case of a set
    xbvue  Holland1976 --params params.txt    # override the params file name
    xbvue  Holland1976 --details              # terminal: full detail for every case
    xbvue  case_folder --json config.json     # also dump to JSON
    xbvue  case_folder --html                 # webpage only, opens in browser
    xbvue  case_folder --html report.html --no-browser

Requires Python 3.8+. Optional: numpy (faster grid reading), and netCDF4,
h5py or scipy to summarise an existing xboutput.nc.
"""
from __future__ import annotations

import argparse
import configparser
import html as _h
import json
import math
import re
import sys
from pathlib import Path

try:
    import numpy as np  # optional
except ImportError:  # pragma: no cover
    np = None

from ._common import _e, h1, h2, PALETTE, strip_private, wrap_pairs

ENC = "latin-1"  # never fails to decode


class Placeholder(Exception):
    """File exists but holds only NUL bytes (e.g. a Dropbox/OneDrive online-only stub)."""


# --------------------------------------------------------------------------- #
# Generic helpers
# --------------------------------------------------------------------------- #
def read_text(path: Path) -> str:
    data = path.read_bytes()
    if data and not data.strip(b"\0"):
        raise Placeholder(f"{path.name} contains only NUL bytes - cloud placeholder not downloaded?")
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode(ENC, "replace")


def find_ci(folder: Path, name: str) -> Path | None:
    """Case-insensitive lookup of a (possibly relative) file name inside folder."""
    if not name:
        return None
    p = folder / name
    if p.exists():
        return p
    parent = p.parent
    if parent.is_dir():
        low = p.name.lower()
        for q in parent.iterdir():
            if q.name.lower() == low:
                return q
    return None


def to_float(v, default=None):
    try:
        return float(str(v).replace("D", "E").replace("d", "e"))
    except (TypeError, ValueError):
        return default


def to_int(v, default=None):
    f = to_float(v)
    return int(f) if f is not None else default


def read_numbers(path: Path) -> list[float]:
    txt = read_text(path).replace("D", "E").replace("d", "e")
    toks = txt.split()
    if np is not None:
        try:
            return np.array(toks, dtype=float)
        except ValueError:
            pass
    return [x for x in (to_float(t) for t in toks) if x is not None]


def read_rows(path: Path) -> list[list[float]]:
    rows = []
    for ln in read_text(path).splitlines():
        ln = ln.split("%", 1)[0].strip()
        if not ln:
            continue
        vals = [to_float(t) for t in ln.replace(",", " ").split()]
        if vals and all(v is not None for v in vals):
            rows.append(vals)
    return rows


def g(v, nd=4):
    """Compact number formatting."""
    if v is None:
        return None
    if isinstance(v, (int,)) or (isinstance(v, float) and v.is_integer() and abs(v) < 1e9):
        return str(int(v))
    return f"{v:.{nd}g}"


def fmt_dur(sec) -> str | None:
    s = to_float(sec)
    if s is None:
        return None
    if s >= 86400:
        return f"{g(s)} s ({s / 3600:.1f} h = {s / 86400:.2f} d)"
    if s >= 3600:
        return f"{g(s)} s ({s / 3600:.2f} h)"
    if s >= 60:
        return f"{g(s)} s ({s / 60:.1f} min)"
    return f"{g(s)} s"


def stats(a) -> dict:
    if a is None or len(a) == 0:
        return {}
    if np is not None:
        a = np.asarray(a, dtype=float)
        return {"min": float(a.min()), "max": float(a.max()), "mean": float(a.mean())}
    return {"min": min(a), "max": max(a), "mean": sum(a) / len(a)}


# --------------------------------------------------------------------------- #
# Keyword catalogue (used to group parameters; unknown keys keep their file section)
# --------------------------------------------------------------------------- #
CATEGORIES = {
    "Physical processes": "wavemodel swave lwave flow sedtrans morphology avalanching nonh roller "
                          "gwflow ships vegetation snells single_dir wind viscosity advection q3d "
                          "swrunup setbathy bedupdate cyclic wci",
    "Grid": "gridform depfile posdwn nx ny alfa vardx dx dy xfile yfile xyfile xori yori thetamin "
            "thetamax dtheta dtheta_s thetanaut projection",
    "Model time": "tstop cfl dtset maxdtfac tunits",
    "Initial conditions": "zs0 zsinitfile hotstartflow hotstartfileno",
    "Wave boundary": "wbctype instat bcfile random rt dtbc taper hrms tm01 trep tlong dir0 m nmax "
                     "fcutoff sprdthr nspr wbcevarreduce wbcremovestokes wbcscaleenergy wbcversion "
                     "nspectrumloc correcthm0 tm01switch oldwbc nonhspectrum_dtheta_s",
    "Flow boundary": "front back left right arc order epsi tidetype freewave lateralwave carspan "
                     "nc highcomp",
    "Water level / tide": "zs0file tideloc tidelen paulrevere",
    "Wind": "windfile windv windth rhoa cd",
    "Flow parameters": "bedfriction bedfriccoef bedfricfile nuh nuhfac nuhv smag umin hmin eps "
                       "c cf fw fwcutoff fwfile",
    "Wave breaking": "break gamma gamma2 alpha n gammax delta breakerdelay",
    "Roller": "beta rfb",
    "Wave numerics": "scheme wavint maxiter maxerror",
    "Non-hydrostatic": "solver kdmin nhbreaker maxbrsteep secbrsteep reformsteep topt dispc nhq3d "
                       "nhlay",
    "Sediment": "d50 d90 rhos dzg dzg1 dzg2 dzg3 ngd nd por sedcal ucrcal struct ne_layer form "
                "waveform facsk facas facua lws sws turb tsmin smax cmax lsgrad bdslpeffmag "
                "bdslpeffdir bdslpeffini sus bed bulk fallvelred",
    "Morphology": "morfac morfacopt morstart morstop wetslp dryslp hswitch dzmax deltahmin "
                  "dryslpfac",
    "Vegetation": "veggiefile veggiemapfile vegnonlin vegcanflo vegstokes vegbreaking",
    "Physical constants": "rho g depthscale lat",
    "Groundwater": "kx ky kz gwdelay aquiferbotfile gw0 gw0file gwhorinfl gwscheme gwnonh",
    "Output": "outputformat tstart tintg tintm tintp tintc tsglobal tsmean tspoints ncfilename "
              "outputprecision nglobalvar nmeanvar npointvar npoints nrugauge nrugdepth rugdepth "
              "remdryoutput timings",
    "Parallel": "mpiboundary mmpi nmpi",
}
KEY2CAT = {k: cat for cat, keys in CATEGORIES.items() for k in keys.split()}
CAT_ORDER = list(CATEGORIES) + ["Other"]

# Values XBeach uses when a key is NOT in params.txt (only those relied on below)
DEFAULTS = {"posdwn": "1", "morfac": "1", "morphology": "1", "sedtrans": "1", "zs0": "0",
            "outputformat": "fortran", "wavemodel": "surfbeat", "nonh": "0", "vardx": "0",
            "cfl": "0.7"}

LIST_KEYS = {"nglobalvar", "nmeanvar", "npointvar", "npoints", "nrugauge"}
SEC_RE = re.compile(r"^%{2,}\s*([A-Za-z][^%:]*?)\s*%{3,}\s*$")
HDR_RE = re.compile(r"^%+\s*(date|function)\s*:\s*(.*?)\s*%*\s*$", re.I)


# --------------------------------------------------------------------------- #
# params.txt
# --------------------------------------------------------------------------- #
def parse_params(path: Path) -> dict:
    lines = read_text(path).splitlines()
    params, lists, header, dups, unparsed = {}, {}, {}, [], []
    section, i = None, 0
    while i < len(lines):
        s = lines[i].strip()
        i += 1
        if not s:
            continue
        if s.startswith("%"):
            m = HDR_RE.match(s)
            if m:
                header[m.group(1).lower()] = m.group(2)
                continue
            m = SEC_RE.match(s)
            if m and not m.group(1).lower().startswith("xbeach parameter"):
                section = m.group(1)
            continue
        s = s.split("%", 1)[0].strip()
        if "=" not in s:
            unparsed.append(s)
            continue
        k, v = (x.strip() for x in s.split("=", 1))
        kl = k.lower()
        if kl in params:
            dups.append(k)
        params[kl] = {"key": k, "value": v, "section": section}
        if kl in LIST_KEYS:
            n, items = to_int(v, 0), []
            while len(items) < n and i < len(lines):
                t = lines[i].split("%", 1)[0].strip()
                if t and "=" in t:
                    break
                i += 1
                if t:
                    items.append(t)
            lists[kl] = items
    return {"file": path.name, "header": header, "params": params, "lists": lists,
            "duplicates": dups, "unparsed": unparsed}


class P:
    """Case-insensitive accessor for parsed parameters."""

    def __init__(self, parsed):
        self.d = parsed["params"]

    def __call__(self, key, default=None):
        e = self.d.get(key.lower())
        return e["value"] if e is not None and e["value"] != "" else default

    def has(self, key):
        return key.lower() in self.d

    def eff(self, key):
        """(value, is_default)"""
        v = self(key)
        return (v, False) if v is not None else (DEFAULTS.get(key.lower()), True)


def grouped_params(parsed) -> dict:
    out = {c: [] for c in CAT_ORDER}
    for kl, e in parsed["params"].items():
        cat = KEY2CAT.get(kl) or "Other"
        out[cat].append({"key": e["key"], "value": e["value"], "file_section": e["section"]})
    return {c: v for c, v in out.items() if v}


# --------------------------------------------------------------------------- #
# Grid and bathymetry
# --------------------------------------------------------------------------- #
def _diffs(a, nxp, nyp, axis):
    if np is not None:
        A = np.asarray(a, dtype=float).reshape(nyp, nxp)
        d = np.abs(np.diff(A, axis=axis))
        return d.ravel() if d.size else None
    rows = [a[j * nxp:(j + 1) * nxp] for j in range(nyp)]
    if axis == 1:
        return [abs(r[k + 1] - r[k]) for r in rows for k in range(nxp - 1)] or None
    return [abs(rows[j + 1][k] - rows[j][k]) for j in range(nyp - 1) for k in range(nxp)] or None


def read_grid(folder: Path, p: P, warnings: list) -> dict:
    nx, ny = to_int(p("nx"), 0), to_int(p("ny"), 0)
    nxp, nyp = nx + 1, ny + 1
    npts = nxp * nyp
    posdwn, posdwn_def = p.eff("posdwn")
    posdwn = to_float(posdwn, 1)
    vardx, _ = p.eff("vardx")
    out = {"nx": nx, "ny": ny, "dimension": "1D (cross-shore profile)" if ny == 0 else "2D",
           "grid_points": f"{nxp} x {nyp} = {npts}",
           "vardx": "variable (xfile/yfile)" if vardx == "1" else "uniform (dx/dy)",
           "origin_xori_yori": f"{p('xori', '0')}, {p('yori', '0')}",
           "rotation_alfa": p("alfa", "0"),
           "posdwn": f"{g(posdwn)}{' (default)' if posdwn_def else ''} -> "
                     + ("depth positive DOWN in depfile" if posdwn == 1 else "bed level positive UP in depfile"),
           "_plot": {}}

    def load(key):
        fn = p(key)
        if not fn:
            return None
        fp = find_ci(folder, fn)
        if fp is None:
            warnings.append(f"{key} '{fn}' not found")
            return None
        try:
            a = read_numbers(fp)
        except Placeholder as e:
            warnings.append(str(e))
            return None
        if len(a) != npts:
            warnings.append(f"{fn}: {len(a)} values, expected (nx+1)*(ny+1) = {npts}")
            if len(a) < npts:
                return None
            a = a[:npts]
        return a

    x = y = None
    if vardx == "1":
        x, y = load("xfile"), load("yfile")
    else:
        dx, dy = to_float(p("dx")), to_float(p("dy"))
        if dx:
            x = [i * dx for _ in range(nyp) for i in range(nxp)]
            out["dx_dy"] = f"{g(dx)} / {g(dy) if dy else '-'}"
        if dy and ny:
            y = [j * dy for j in range(nyp) for _ in range(nxp)]
    if x is not None:
        sx = stats(x)
        out["x_range"] = f"{g(sx['min'], 6)} .. {g(sx['max'], 6)}  (length {g(sx['max'] - sx['min'], 6)})"
        d = _diffs(x, nxp, nyp, 1)
        if d is not None:
            sd = stats(d)
            out["dx_min_max"] = f"{g(sd['min'])} .. {g(sd['max'])}"
    if y is not None and ny:
        sy = stats(y)
        out["y_range"] = f"{g(sy['min'], 6)} .. {g(sy['max'], 6)}  (width {g(sy['max'] - sy['min'], 6)})"
        d = _diffs(y, nxp, nyp, 0)
        if d is not None:
            sd = stats(d)
            out["dy_min_max"] = f"{g(sd['min'])} .. {g(sd['max'])}"

    dep = load("depfile")
    if dep is not None:
        zb = [-v for v in dep] if posdwn == 1 else list(dep)
        if np is not None:
            zb = np.asarray(zb, dtype=float)
        s = stats(zb)
        out["zb_range"] = f"{g(s['min'])} .. {g(s['max'])} m (mean {g(s['mean'])})"
        out["zb_min"], out["zb_max"] = s["min"], s["max"]
        out["_zb"] = zb
        if ny == 0 and x is not None:
            out["_plot"]["profile"] = (list(map(float, x)), list(map(float, zb)))
        elif ny > 0:
            sj, si = max(1, math.ceil(nyp / 220)), max(1, math.ceil(nxp / 220))
            grid = [[float(zb[j * nxp + i]) for i in range(0, nxp, si)] for j in range(0, nyp, sj)]
            out["_plot"]["map"] = grid
    elif not p("depfile"):
        warnings.append("no depfile set")
    return out


# --------------------------------------------------------------------------- #
# Boundary conditions
# --------------------------------------------------------------------------- #
JONS_TABLE_COLS = ["Hm0", "Tp", "mainang", "gammajsp", "s", "duration", "dtbc"]


def parse_wave_bc(folder: Path, p: P, warnings: list) -> dict:
    wbctype = p("wbctype") or p("instat")
    bcf = p("bcfile")
    out = {"type": wbctype or "not set", "type_key": "wbctype" if p("wbctype") else
           ("instat" if p("instat") else None), "bcfile": bcf,
           "random": p("random"), "rt": p("rt"), "dtbc": p("dtbc")}
    out["params_in_params_txt"] = {k: p(k) for k in ("Hrms", "Tm01", "Trep", "Tlong", "dir0", "m", "taper")
                                   if p(k) is not None}
    if not bcf:
        return out
    fp = find_ci(folder, bcf)
    if fp is None:
        warnings.append(f"bcfile '{bcf}' not found")
        return out
    try:
        txt = read_text(fp)
    except Placeholder as e:
        warnings.append(str(e))
        return out
    lines = [ln.split("%", 1)[0].strip() for ln in txt.splitlines()]
    lines = [ln for ln in lines if ln]
    if lines and lines[0].upper().startswith("FILELIST"):
        files = []
        for ln in lines[1:]:
            f = ln.split()
            if len(f) >= 3:
                files.append({"duration": f[0], "dtbc": f[1], "file": f[2],
                              "exists": find_ci(folder, f[2]) is not None})
        out["format"] = "FILELIST"
        out["filelist"] = files
        out["total_duration"] = fmt_dur(sum(to_float(f["duration"], 0) for f in files))
        miss = [f["file"] for f in files if not f["exists"]]
        if miss:
            warnings.append(f"FILELIST entries missing: {', '.join(miss[:5])}{' ...' if len(miss) > 5 else ''}")
    elif any("=" in ln for ln in lines):
        out["format"] = "parametric (key = value)"
        out["spectrum"] = dict(tuple(x.strip() for x in ln.split("=", 1)) for ln in lines if "=" in ln)
    else:
        rows = read_rows(fp)
        out["format"] = f"table, {len(rows)} row(s) x {len(rows[0]) if rows else 0} col(s)"
        is_jt = (wbctype or "").lower() in ("jons_table", "jonstable") and rows and len(rows[0]) >= 6
        cols = JONS_TABLE_COLS[:len(rows[0])] if is_jt else [f"col{i + 1}" for i in range(len(rows[0]) if rows else 0)]
        out["columns"] = cols
        out["rows"] = rows
        if is_jt:
            hm0 = [r[0] for r in rows]
            tp = [r[1] for r in rows]
            dur = [r[5] for r in rows]
            out["Hm0_range"] = f"{g(min(hm0))} .. {g(max(hm0))} m"
            out["Tp_range"] = f"{g(min(tp))} .. {g(max(tp))} s"
            out["direction"] = ", ".join(sorted({g(r[2]) for r in rows}))
            out["total_duration"] = fmt_dur(sum(dur))
            t, ts, h, tps = 0.0, [], [], []
            for r in rows:
                ts += [t / 3600, (t + r[5]) / 3600]
                h += [r[0], r[0]]
                tps += [r[1], r[1]]
                t += r[5]
            out["_plot"] = {"t_h": ts, "Hm0": h, "Tp": tps}
            tstop = to_float(p("tstop"))
            if tstop and t < tstop:
                warnings.append(f"wave table covers {g(t)} s < tstop {g(tstop)} s (table is recycled)")
    return out


def parse_water_level(folder: Path, p: P, warnings: list) -> dict:
    zs0, zs0_def = p.eff("zs0")
    out = {"zs0": f"{zs0}{' (default)' if zs0_def else ''}", "tidetype": p("tidetype")}
    f = p("zs0file")
    if not f:
        out["mode"] = "constant (zs0)"
        return out
    out["mode"] = "time-varying (zs0file)"
    out["zs0file"] = f
    fp = find_ci(folder, f)
    if fp is None:
        warnings.append(f"zs0file '{f}' not found")
        return out
    try:
        rows = read_rows(fp)
    except Placeholder as e:
        warnings.append(str(e))
        return out
    if not rows:
        return out
    ncol = len(rows[0])
    tideloc = to_int(p("tideloc"), ncol - 1)
    out["tideloc"] = f"{tideloc}{'' if p('tideloc') else ' (from file columns)'}"
    if p("tidelen"):
        out["tidelen"] = p("tidelen")
    t = [r[0] for r in rows]
    out["records"] = len(rows)
    out["time_span"] = f"{g(t[0])} .. {g(t[-1])} s ({(t[-1] - t[0]) / 3600:.2f} h), dt ~ {g((t[-1] - t[0]) / max(1, len(t) - 1))} s"
    names = {1: ["offshore"], 2: ["offshore (front)", "landward (back)"],
             4: ["offshore left", "offshore right", "back right", "back left"]}.get(tideloc, [])
    locs, series = [], []
    for c in range(1, ncol):
        v = [r[c] for r in rows if len(r) > c]
        nm = names[c - 1] if c - 1 < len(names) else f"loc {c}"
        locs.append({"column": c + 1, "location": nm, "min": g(min(v)), "max": g(max(v))})
        series.append((nm, v))
    out["locations"] = locs
    out["_plot"] = {"t_h": [x / 3600 for x in t], "series": series}
    tstop = to_float(p("tstop"))
    if tstop and t[-1] < tstop:
        warnings.append(f"zs0file ends at {g(t[-1])} s, before tstop = {g(tstop)} s")
    return out


def parse_vegetation(folder: Path, p: P, warnings: list) -> dict | None:
    if p("vegetation", "0") != "1" and not p("veggiefile"):
        return None
    out = {"vegetation": p("vegetation"), "veggiefile": p("veggiefile"),
           "veggiemapfile": p("veggiemapfile"), "species": []}
    vf = find_ci(folder, p("veggiefile") or "")
    if p("veggiefile") and vf is None:
        warnings.append(f"veggiefile '{p('veggiefile')}' not found")
    if vf:
        try:
            names = [ln.strip() for ln in read_text(vf).splitlines() if ln.strip()]
        except Placeholder as e:
            warnings.append(str(e))
            names = []
        for nm in names:
            sp = {"file": nm}
            sf = find_ci(folder, nm)
            if sf is None:
                sp["missing"] = True
                warnings.append(f"vegetation species file '{nm}' not found")
            else:
                try:
                    for ln in read_text(sf).splitlines():
                        ln = ln.split("%", 1)[0]
                        if "=" in ln:
                            k, v = ln.split("=", 1)
                            sp[k.strip()] = v.strip()
                except Placeholder as e:
                    warnings.append(str(e))
            out["species"].append(sp)
    mf = find_ci(folder, p("veggiemapfile") or "")
    if p("veggiemapfile") and mf is None:
        warnings.append(f"veggiemapfile '{p('veggiemapfile')}' not found")
    if mf:
        try:
            a = [int(round(v)) for v in read_numbers(mf)]
            cnt = {}
            for v in a:
                cnt[v] = cnt.get(v, 0) + 1
            out["map_cells_per_species"] = {("none" if k == 0 else f"species {k}"): n for k, n in sorted(cnt.items())}
            out["_map"] = a
        except Placeholder as e:
            warnings.append(str(e))
    return out


# --------------------------------------------------------------------------- #
# Referenced files, output request, results
# --------------------------------------------------------------------------- #
FILE_KEYS = {"ne_layer", "bcfile"}


def referenced_files(folder: Path, parsed) -> list[dict]:
    out = []
    for kl, e in parsed["params"].items():
        if (kl.endswith("file") or kl in FILE_KEYS) and e["value"] and kl not in ("ncfilename",):
            out.append({"key": e["key"], "file": e["value"],
                        "exists": find_ci(folder, e["value"]) is not None})
    return out


def output_request(p: P, parsed) -> dict:
    fmt, d = p.eff("outputformat")
    L = parsed["lists"]
    out = {"outputformat": f"{fmt}{' (default)' if d else ''}", "ncfilename": p("ncfilename"),
           "tstart": fmt_dur(p("tstart")),
           "tintg (global)": fmt_dur(p("tintg")), "tintm (mean)": fmt_dur(p("tintm")),
           "tintp (points)": fmt_dur(p("tintp")), "tintc": fmt_dur(p("tintc")),
           "global_vars": L.get("nglobalvar", []), "mean_vars": L.get("nmeanvar", []),
           "point_vars": L.get("npointvar", []), "points": L.get("npoints", []),
           "runup_gauges": L.get("nrugauge", [])}
    if p("tstop") and p("tstart"):
        tg = to_float(p("tintg"))
        if tg:
            out["global_snapshots"] = int((to_float(p("tstop")) - to_float(p("tstart"))) / tg) + 1
    return out


def _nc_summary(fp: Path) -> dict:
    out = {}
    try:
        import netCDF4  # type: ignore
        with netCDF4.Dataset(fp) as ds:
            out["dimensions"] = {k: len(v) for k, v in ds.dimensions.items()}
            out["variables"] = [k for k in ds.variables if k not in ds.dimensions]
            for tn in ("globaltime", "time"):
                if tn in ds.variables and len(ds.variables[tn]):
                    out["last_global_time"] = float(ds.variables[tn][-1])
                    break
        return out
    except ImportError:
        pass
    except Exception as e:
        out["error"] = str(e)
        return out
    try:
        from scipy.io import netcdf_file  # netCDF3 classic
        with netcdf_file(fp, "r", mmap=False) as ds:
            out["dimensions"] = {k: (v if v is not None else "unlimited") for k, v in ds.dimensions.items()}
            out["variables"] = [k for k in ds.variables if k not in ds.dimensions]
            if "globaltime" in ds.variables:
                t = ds.variables["globaltime"].data
                out["dimensions"]["globaltime"] = len(t)
                if len(t):
                    out["last_global_time"] = float(t[-1])
        return out
    except ImportError:
        pass
    except Exception:
        pass
    try:
        import h5py  # netCDF4/HDF5 files only
        with h5py.File(fp, "r") as f:
            out["variables"] = list(f.keys())
            if "globaltime" in f and len(f["globaltime"]):
                out["last_global_time"] = float(f["globaltime"][-1])
    except Exception:
        out["note"] = "install netCDF4 (or scipy) for variables / times"
    return out


def results_status(folder: Path, p: P) -> dict | None:
    out = {}
    nc = find_ci(folder, p("ncfilename") or "xboutput.nc")
    if nc:
        out["netcdf"] = f"{nc.name} ({nc.stat().st_size / 1e6:.1f} MB)"
        info = _nc_summary(nc)
        if "last_global_time" in info:
            t = info.pop("last_global_time")
            out["last_output_time"] = fmt_dur(t)
            tstop = to_float(p("tstop"))
            if tstop:
                out["completed"] = "yes" if t >= tstop - 1e-6 else f"NO - stopped at {100 * t / tstop:.0f} % of tstop"
        out.update(info)
    dats = sorted(x.name for x in folder.glob("*.dat"))
    if dats:
        out["fortran_dat_files"] = f"{len(dats)} (.dat): {', '.join(dats[:8])}{' ...' if len(dats) > 8 else ''}"
    for nm in ("XBlog.txt", "XBerror.txt", "XBwarning.txt"):
        f = find_ci(folder, nm)
        if f is None:
            continue
        try:
            tail = [ln.strip() for ln in read_text(f).splitlines() if ln.strip()]
        except Placeholder:
            tail = ["(cloud placeholder)"]
        if nm == "XBlog.txt":
            out["XBlog_tail"] = tail[-3:]
        elif tail:
            out[nm] = tail[:3]
    return out or None


# --------------------------------------------------------------------------- #
# Model summary
# --------------------------------------------------------------------------- #
def model_summary(p: P, grid: dict, waves: dict, wl: dict, veg) -> dict:
    wm, wm_d = p.eff("wavemodel")
    if wm_d and p("nonh") == "1":
        wm, note = "nonh", " (from nonh = 1)"
    else:
        note = " (default)" if wm_d else ""
    morph, md = p.eff("morphology")
    sed, sd = p.eff("sedtrans")
    mf, _ = p.eff("morfac")
    cfl, cd = p.eff("cfl")
    morph_txt = ("ON" if morph == "1" else "off") + (" (default)" if md else "")
    if morph == "1":
        morph_txt += f", morfac {mf}" + (f", morstart {p('morstart')}" if p("morstart") else "")
    procs = [k for k in ("roller", "avalanching", "gwflow", "wind", "vegetation", "ships", "nonh", "snells", "single_dir")
             if p(k) == "1"]
    out = {
        "Wave model": wm + note,
        "Dimension": f"{grid.get('dimension')}, nx = {grid.get('nx')}, ny = {grid.get('ny')}",
        "Simulation time (tstop)": fmt_dur(p("tstop")),
        "CFL": f"{cfl}{' (default)' if cd else ''}",
        "Wave forcing": f"{waves['type']}" + (f" <- {waves['bcfile']}" if waves.get("bcfile") else ""),
        "Water level": wl["mode"] + (f" <- {wl['zs0file']}" if wl.get("zs0file") else f" = {wl['zs0']}"),
        "Boundaries front/back/l/r": (" / ".join(p(k, "-") for k in ("front", "back", "left", "right"))
                                      if any(p(k) for k in ("front", "back", "left", "right")) else None),
        "Sediment transport": ("ON" if sed == "1" else "off") + (" (default)" if sd else ""),
        "Morphology": morph_txt,
        "Breaking": " ".join(x for x in (p("break"), p("gamma") and f"gamma={p('gamma')}",
                                         p("alpha") and f"alpha={p('alpha')}") if x) or None,
        "Bed friction": " ".join(x for x in (p("bedfriction"), p("bedfriccoef")) if x) or None,
        "Grain size D50 / D90": " / ".join(x for x in (p("D50"), p("D90")) if x) or None,
        "Switched-on processes": ", ".join(procs) or None,
        "Vegetation species": ", ".join(s["file"] for s in veg["species"]) if veg else None,
    }
    return {k: v for k, v in out.items() if v not in (None, "")}


# --------------------------------------------------------------------------- #
# Case / case-set discovery
# --------------------------------------------------------------------------- #
def parse_oet_config(path: Path) -> dict | None:
    cp = configparser.ConfigParser(interpolation=None, strict=False)
    try:
        cp.read_string(read_text(path))
    except (configparser.Error, Placeholder):
        return None
    gen = dict(cp["general"]) if cp.has_section("general") else {}
    names = [r.strip() for r in gen.get("runs", "").split(",") if r.strip()]
    sec_runs = [s[4:] for s in cp.sections() if s.startswith("run_")]
    names += [r for r in sec_runs if r not in names]
    runs = [{"name": r, **(dict(cp["run_" + r]) if cp.has_section("run_" + r) else {})} for r in names]
    return {"file": path.name, "general": gen, "runs": runs,
            "categories": dict(cp["categories"]) if cp.has_section("categories") else {}}


def discover(arg: str, params_name: str | None) -> dict:
    p = Path(arg).expanduser().resolve()
    if not p.exists():
        sys.exit(f"Path not found: {p}")
    if p.is_file() and p.name != ".config":
        return {"root": p.parent, "name": p.parent.name, "config": None,
                "cases": [{"name": p.parent.name, "folder": p.parent, "params": p.name, "enabled": True}]}
    root = p if p.is_dir() else p.parent
    pname = params_name or "params.txt"
    cfg = parse_oet_config(root / ".config") if (root / ".config").exists() else None
    cases = []
    if cfg and cfg["runs"]:
        for r in cfg["runs"]:
            sub = r.get("path", ".")
            folder = (root / sub).resolve()
            cases.append({"name": r["name"] if sub not in (".", "") else root.name,
                          "folder": folder, "params": params_name or r.get("params", "params.txt"),
                          "enabled": r.get("enable", "1") == "1",
                          "config_run": {k: v for k, v in r.items() if k not in ("name", "path", "params")}})
    elif find_ci(root, pname):
        cases.append({"name": root.name, "folder": root, "params": pname, "enabled": True})
    else:
        for f in sorted(root.rglob("*")):
            if f.is_file() and f.name.lower() == pname.lower() and len(f.relative_to(root).parts) <= 4 \
                    and not any(x.startswith(".") for x in f.relative_to(root).parts):
                cases.append({"name": str(f.parent.relative_to(root)).replace("\\", "/"),
                              "folder": f.parent, "params": f.name, "enabled": True})
    if not cases:
        sys.exit(f"No XBeach '{pname}' or OET .config found in {root}")
    return {"root": root, "name": root.name, "config": cfg, "cases": cases}


def readme_text(root: Path) -> str | None:
    for nm in ("README.md", "README.txt", "readme.md", "README"):
        f = root / nm
        if f.exists():
            try:
                return read_text(f)[:20000]
            except Placeholder:
                return None
    return None


def build_case(c: dict) -> dict:
    folder = c["folder"]
    warnings = []
    case = {"name": c["name"], "folder": str(folder), "params_file": c["params"],
            "enabled": c["enabled"], "config_run": c.get("config_run")}
    pf = find_ci(folder, c["params"])
    if pf is None:
        case["error"] = f"{c['params']} not found in {folder}"
        return case
    try:
        parsed = parse_params(pf)
    except Placeholder as e:
        case["error"] = str(e)
        return case
    p = P(parsed)
    case["header"] = parsed["header"]
    case["other_params_files"] = sorted(x.name for x in folder.glob("params*.txt") if x.name != pf.name)
    case["n_parameters"] = len(parsed["params"])
    case["params"] = {e["key"]: e["value"] for e in parsed["params"].values()}
    case["parameters_by_category"] = grouped_params(parsed)
    case["grid"] = read_grid(folder, p, warnings)
    case["waves"] = parse_wave_bc(folder, p, warnings)
    case["water_level"] = parse_water_level(folder, p, warnings)
    zb, wl = case["grid"].pop("_zb", None), case["water_level"]
    if zb is not None:
        pl = wl.get("_plot")
        z0, lab = (pl["series"][0][1][0], "initial offshore zs") if pl and pl["series"] \
            else (to_float(p.eff("zs0")[0], 0.0), "zs0")
        wet = sum(1 for v in zb if v < z0) if np is None else int((zb < z0).sum())
        case["grid"]["wet_initially"] = f"{wet} of {len(zb)} points ({100 * wet / len(zb):.0f} %) below {lab} = {g(z0)}"
    case["vegetation"] = parse_vegetation(folder, p, warnings)
    case["output"] = output_request(p, parsed)
    case["files"] = referenced_files(folder, parsed)
    case["results"] = results_status(folder, p)
    case["summary"] = model_summary(p, case["grid"], case["waves"], case["water_level"], case["vegetation"])
    if parsed["lists"].get("npointvar") and not parsed["lists"].get("npoints") and not parsed["lists"].get("nrugauge"):
        warnings.append("npointvar requested but no npoints / nrugauge locations given")
    if parsed["duplicates"]:
        warnings.append(f"keys set more than once (last wins): {', '.join(parsed['duplicates'])}")
    if parsed["unparsed"]:
        warnings.append(f"lines not understood: {parsed['unparsed'][:5]}")
    for f in case["files"]:
        if not f["exists"] and not any(f["file"] in w for w in warnings):
            warnings.append(f"{f['key']} '{f['file']}' not found")
    case["warnings"] = warnings
    return case


def differing_params(cases: list[dict]) -> list[str]:
    ok = [c for c in cases if "params" in c]
    if len(ok) < 2:
        return []
    keys = []
    for c in ok:
        keys += [k for k in c["params"] if k.lower() not in [x.lower() for x in keys]]

    def val(c, k):
        for kk, v in c["params"].items():
            if kk.lower() == k.lower():
                f = to_float(v)
                return f if f is not None else v.lower()
        return None

    return [k for k in keys if len({val(c, k) for c in ok}) > 1]


def param_ci(c, k):
    for kk, v in c.get("params", {}).items():
        if kk.lower() == k.lower():
            return v
    return None


def build_config(arg: str, params_name=None, only_case=None) -> dict:
    d = discover(arg, params_name)
    sel = d["cases"]
    if only_case:
        sel = [c for c in sel if c["name"].lower() == only_case.lower()]
        if not sel:
            sys.exit(f"Case '{only_case}' not found. Cases: {', '.join(c['name'] for c in d['cases'])}")
    cases = [build_case(c) for c in sel]
    return {"name": d["name"], "root": str(d["root"]), "oet_config": d["config"],
            "readme": readme_text(d["root"]) if len(d["cases"]) > 1 else None,
            "n_cases_total": len(d["cases"]), "cases": cases,
            "differing_parameters": differing_params(cases)}


# --------------------------------------------------------------------------- #
# Terminal output
# --------------------------------------------------------------------------- #


def row(label, value, indent=2):
    if value in (None, "", [], {}):
        return
    if isinstance(value, dict):
        value = ", ".join(f"{k}: {v}" for k, v in value.items())
    elif isinstance(value, list):
        value = ", ".join(map(str, value))
    print(f"{' ' * indent}{label:<28}: {value}")


def text_table(headers, rows):
    w = [max(len(str(h)), *(len(str(r[i])) for r in rows)) for i, h in enumerate(headers)]
    print("  " + "  ".join(str(h).ljust(w[i]) for i, h in enumerate(headers)))
    print("  " + "  ".join("-" * x for x in w))
    for r in rows:
        print("  " + "  ".join(str(c).ljust(w[i]) for i, c in enumerate(r)))


def case_row(c):
    if "error" in c:
        return [c["name"], c["params_file"], "ERROR: " + c["error"].split(" in /")[0].split(" in ")[0][:50],
                "", "", "", "", ""]
    gr, wv = c["grid"], c["waves"]
    return [c["name"] + ("" if c["enabled"] else " (off)"), c["params_file"],
            f"{gr['nx']}x{gr['ny']}", g(to_float(param_ci(c, "tstop") or 0) / 3600, 3) + " h",
            wv["type"], c["summary"].get("Morphology", "").split(",")[0],
            f"{g(gr.get('zb_min'))}..{g(gr.get('zb_max'))}" if "zb_min" in gr else "-",
            "yes" if c.get("results") else "-"]


def show_case(c):
    h1(f"XBEACH CASE: {c['name']}" + ("" if c["enabled"] else "   [disabled in .config]"))
    row("Folder", c["folder"])
    row("Params file", c["params_file"])
    if "error" in c:
        print(f"  ERROR: {c['error']}")
        return
    row("Written", " / ".join(f"{k}: {v}" for k, v in c["header"].items()))
    row("Other params files", c["other_params_files"])
    row("OET run settings", c.get("config_run"))
    h2("Model summary")
    for k, v in c["summary"].items():
        row(k, v)
    h2("Grid & bathymetry")
    for k, v in c["grid"].items():
        if not k.startswith("_") and k not in ("zb_min", "zb_max", "nx", "ny"):
            row(k, v)
    w = c["waves"]
    h2(f"Wave boundary ({w['type_key'] or 'type'} = {w['type']})")
    for k in ("bcfile", "format", "random", "rt", "dtbc", "Hm0_range", "Tp_range", "direction", "total_duration"):
        row(k, w.get(k))
    row("Spectrum", w.get("spectrum"))
    row("Wave params in params.txt", w.get("params_in_params_txt"))
    if w.get("rows") and w.get("columns"):
        print("    " + "  ".join(f"{x:>9}" for x in w["columns"]))
        for r in w["rows"][:6]:
            print("    " + "  ".join(f"{g(v):>9}" for v in r))
        if len(w["rows"]) > 6:
            print(f"    ... {len(w['rows']) - 6} more row(s)")
    for f in w.get("filelist", [])[:10]:
        print(f"    {f['file']}  duration {f['duration']}  dtbc {f['dtbc']}{'' if f['exists'] else '  MISSING'}")
    wl = c["water_level"]
    h2("Water level")
    for k in ("mode", "zs0", "zs0file", "tidetype", "tideloc", "tidelen", "records", "time_span"):
        row(k, wl.get(k))
    for L in wl.get("locations", []):
        print(f"    col {L['column']} {L['location']:<18} min {L['min']:>8}  max {L['max']:>8}")
    v = c["vegetation"]
    if v:
        h2("Vegetation")
        row("veggiefile / veggiemapfile", f"{v['veggiefile']} / {v['veggiemapfile']}")
        for s in v["species"]:
            row(s["file"], {k: x for k, x in s.items() if k != "file"}, indent=4)
        row("Map cells", v.get("map_cells_per_species"))
    h2(f"Parameters ({c['n_parameters']} set in {c['params_file']})")
    for cat, items in c["parameters_by_category"].items():
        print(f"  [{cat}]")
        wrap_pairs([(e["key"], e["value"]) for e in items])
    o = c["output"]
    h2("Output request")
    for k, val in o.items():
        row(k, val)
    h2("Referenced files")
    for f in c["files"]:
        print(f"  {f['key']:<16} {f['file']:<24} {'ok' if f['exists'] else 'MISSING'}")
    h2("Results")
    if c["results"]:
        for k, val in c["results"].items():
            row(k, val)
    else:
        print("  No xboutput.nc / XBlog.txt found (case not run here yet).")
    if c["warnings"]:
        h2("Warnings")
        for x in c["warnings"]:
            print(f"  ! {x}")


def show(cfg: dict, details: bool):
    cases = cfg["cases"]
    oc = cfg["oet_config"]
    if len(cases) > 1:
        h1(f"XBEACH CASE SET: {cfg['name']}")
        row("Folder", cfg["root"])
        row("Cases", f"{len(cases)} shown of {cfg['n_cases_total']}"
            + (f" ({sum(c['enabled'] for c in cases)} enabled)" if oc else ""))
        if oc:
            row("OET .config", oc["file"])
            row("Responsible", oc["general"].get("responsible"))
            row("Analysis", oc["general"].get("analysis"))
            row("Categories", oc["categories"])
        if cfg["readme"]:
            first = " ".join(x for x in cfg["readme"].split("\n\n")[1:2] if x).replace("\n", " ")
            row("README", first[:300] + ("..." if len(first) > 300 else ""))
        h2("Cases")
        text_table(["Case", "Params", "nx x ny", "tstop", "Waves", "Morph", "zb range", "Results"],
                   [case_row(c) for c in cases])
        diff = cfg["differing_parameters"]
        if diff:
            h2(f"Parameters that differ between cases ({len(diff)})")
            shown = diff[:8]
            text_table(["Case"] + shown, [[c["name"]] + [param_ci(c, k) or "-" for k in shown]
                                          for c in cases if "params" in c])
            if len(diff) > 8:
                print(f"  ... and {', '.join(diff[8:])}")
        elif len(cases) > 1:
            print("\n  All cases use identical parameter values.")
        if not details and len(cases) > 1:
            print("\n  (Use --case NAME or --details for full per-case configuration, or --html.)")
    if details or len(cases) == 1:
        for c in cases:
            show_case(c)
    print()


# --------------------------------------------------------------------------- #
# SVG charts for the HTML report
# --------------------------------------------------------------------------- #


def _ticks(lo, hi, n=5):
    if hi == lo:
        hi = lo + 1
    span = hi - lo
    step = 10 ** math.floor(math.log10(span / n))
    for m in (1, 2, 2.5, 5, 10):
        if span / (step * m) <= n:
            step *= m
            break
    t = math.ceil(lo / step) * step
    out = []
    while t <= hi + 1e-9 * span:
        out.append(round(t, 10))
        t += step
    return out


def _thin(xs, ys, nmax=1200):
    s = max(1, len(xs) // nmax)
    return xs[::s], ys[::s]


def svg_plot(series, xlabel="", ylabel="", w=820, h=250, hlines=(), legend=True) -> str:
    """series: list of dict(x, y, label, color, fill)."""
    allx = [v for s in series for v in s["x"]]
    ally = [v for s in series for v in s["y"]] + [hl[0] for hl in hlines]
    if not allx:
        return ""
    x0, x1, y0, y1 = min(allx), max(allx), min(ally), max(ally)
    pad = (y1 - y0) * 0.06 or 0.5
    y0, y1 = y0 - pad, y1 + pad
    if x1 == x0:
        x1 = x0 + 1
    L, R, T, B = 58, 14, 12, 38
    X = lambda v: L + (v - x0) / (x1 - x0) * (w - L - R)
    Y = lambda v: T + (y1 - v) / (y1 - y0) * (h - T - B)
    o = [f'<svg class="plot" viewBox="0 0 {w} {h}" preserveAspectRatio="xMidYMid meet">']
    for t in _ticks(y0, y1):
        o.append(f'<line x1="{L}" x2="{w - R}" y1="{Y(t):.1f}" y2="{Y(t):.1f}" class="gl"/>'
                 f'<text x="{L - 6}" y="{Y(t) + 4:.1f}" text-anchor="end">{g(t)}</text>')
    for t in _ticks(x0, x1, 8):
        o.append(f'<line x1="{X(t):.1f}" x2="{X(t):.1f}" y1="{h - B}" y2="{h - B + 4}" class="ax"/>'
                 f'<text x="{X(t):.1f}" y="{h - B + 16}" text-anchor="middle">{g(t)}</text>')
    for yv, lab, col in hlines:
        o.append(f'<line x1="{L}" x2="{w - R}" y1="{Y(yv):.1f}" y2="{Y(yv):.1f}" stroke="{col}" '
                 f'stroke-dasharray="5 4" stroke-width="1.3"/><text x="{w - R - 4}" y="{Y(yv) - 4:.1f}" '
                 f'text-anchor="end" fill="{col}">{_h.escape(lab)}</text>')
    for i, s in enumerate(series):
        xs, ys = _thin(list(s["x"]), list(s["y"]))
        col = s.get("color") or PALETTE[i % len(PALETTE)]
        pts = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in zip(xs, ys))
        tip = f'<title>{_h.escape(s.get("label", ""))}</title>'
        if s.get("fill"):
            poly = f"{X(xs[0]):.1f},{h - B} {pts} {X(xs[-1]):.1f},{h - B}"
            o.append(f'<polygon points="{poly}" fill="{col}" fill-opacity=".35" stroke="none"/>')
        o.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="{s.get("width", 1.6)}">{tip}</polyline>')
    o.append(f'<rect x="{L}" y="{T}" width="{w - L - R}" height="{h - T - B}" class="frame"/>')
    o.append(f'<text x="{(L + w - R) / 2}" y="{h - 4}" text-anchor="middle" class="lab">{_h.escape(xlabel)}</text>')
    o.append(f'<text transform="translate(13,{(T + h - B) / 2}) rotate(-90)" text-anchor="middle" class="lab">{_h.escape(ylabel)}</text>')
    if legend and len(series) > 1 and len(series) <= 8:
        lx = L + 10
        for i, s in enumerate(series):
            col = s.get("color") or PALETTE[i % len(PALETTE)]
            o.append(f'<rect x="{lx}" y="{T + 6}" width="12" height="3" fill="{col}"/>'
                     f'<text x="{lx + 16}" y="{T + 11}">{_h.escape(s.get("label", ""))}</text>')
            lx += 24 + 7 * len(s.get("label", ""))
    o.append("</svg>")
    return "".join(o)


def profile_svg(c) -> str:
    gr = c["grid"]
    pr = gr.get("_plot", {}).get("profile")
    if not pr:
        return ""
    x, zb = pr
    hl = []
    wl = c["water_level"]
    if wl.get("locations"):
        L0 = wl["locations"][0]
        hl += [(to_float(L0["max"]), f"max zs {L0['max']}", "#2f80c0"),
               (to_float(L0["min"]), f"min zs {L0['min']}", "#7fb0e0")]
    else:
        z = to_float(str(wl["zs0"]).split()[0], 0)
        hl.append((z, f"zs0 = {g(z)}", "#2f80c0"))
    series = [{"x": x, "y": zb, "label": "bed level zb", "color": "#b08d57", "fill": True}]
    veg = c.get("vegetation") or {}
    if veg.get("_map") and len(veg["_map"]) == len(x):
        vx = [a for a, m in zip(x, veg["_map"]) if m]
        vz = [b for b, m in zip(zb, veg["_map"]) if m]
        if vx:
            series.append({"x": vx, "y": vz, "label": "vegetated cells", "color": "#3a9d5d", "width": 4})
    return svg_plot(series, "x (m)", "zb (m)", hlines=hl)


def map_canvas(c) -> str:
    grid = c["grid"].get("_plot", {}).get("map")
    if not grid:
        return ""
    flat = [v for r in grid for v in r]
    lo, hi = min(flat), max(flat)
    data = _h.escape(json.dumps([[round(v, 2) for v in r] for r in grid]))
    return (f'<canvas class="zb" data-z="{data}" data-lo="{lo}" data-hi="{hi}"></canvas>'
            f'<div class="muted" style="font-size:12px">bed level zb, grid index space '
            f'(row 0 at top); blue below 0, sand above 0; range {g(lo)} .. {g(hi)} m'
            f'{" (subsampled)" if len(grid[0]) < c["grid"]["nx"] + 1 else ""}</div>')


def waves_svg(w) -> str:
    pl = w.get("_plot")
    if not pl:
        return ""
    return svg_plot([{"x": pl["t_h"], "y": pl["Hm0"], "label": "Hm0 (m)", "color": "#2f80c0"},
                     {"x": pl["t_h"], "y": pl["Tp"], "label": "Tp (s)", "color": "#d9622b"}],
                    "time (h)", "Hm0 (m) / Tp (s)", h=210)


def tide_svg(wl) -> str:
    pl = wl.get("_plot")
    if not pl:
        return ""
    return svg_plot([{"x": pl["t_h"], "y": v, "label": nm} for nm, v in pl["series"]],
                    "time (h)", "zs (m)", h=210)


# --------------------------------------------------------------------------- #
# HTML report  (layout matches ras_config.py / hms_config.py)
# --------------------------------------------------------------------------- #
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<script>(function(){var t;try{t=localStorage.getItem('modelvue-theme')}catch(e){}if(t!=='dark')t='light';document.documentElement.dataset.theme=t;})();</script>
<title>XBeach · @@TITLE@@</title>
<style>
:root{--bg:#f6f7f9;--panel:#fff;--ink:#1d2433;--muted:#667085;--line:#e4e7ec;--acc:#0e5a6b;--accs:#e3f1f3;--warn:#b42318}
:root[data-theme=dark]{color-scheme:dark;--bg:#0f1318;--panel:#171c23;--ink:#e6e9ef;--muted:#98a2b3;--line:#2a313b;--acc:#6cc4d4;--accs:#16313a;--warn:#f97066}
*{box-sizing:border-box} body{margin:0;font:14px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,"Noto Sans TC",sans-serif;background:var(--bg);color:var(--ink)}
header{position:sticky;top:0;z-index:5;background:var(--acc);color:#fff;padding:14px 24px;display:flex;gap:16px;align-items:baseline;flex-wrap:wrap}
header h1{margin:0;font-size:20px} header span{opacity:.85}
[data-theme=dark] header{background:#123c47}
header .theme{margin-left:auto;font:inherit;font-size:12px;padding:3px 12px;border:1px solid rgba(255,255,255,.55);border-radius:14px;background:transparent;color:#fff;cursor:pointer} header .theme:hover{background:rgba(255,255,255,.15)}

.wrap{display:grid;grid-template-columns:240px 1fr;gap:0}
nav{position:sticky;top:56px;align-self:start;height:calc(100vh - 56px);overflow:auto;padding:16px 10px;border-right:1px solid var(--line)}
nav a{display:block;padding:6px 10px;border-radius:6px;color:var(--ink);text-decoration:none;font-size:13px}
nav a:hover,nav a.on{background:var(--accs);color:var(--acc)}
nav .grp{margin:14px 10px 4px;font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}
nav small{color:var(--muted)} nav .off{opacity:.55}
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
td{white-space:nowrap} td.wrap{white-space:normal}
thead th{background:var(--accs);color:var(--acc);cursor:pointer;user-select:none;white-space:nowrap}
thead th:hover{text-decoration:underline} tbody tr:hover{background:var(--accs)}
table.kv{width:auto;min-width:50%} table.kv th{color:var(--muted);font-weight:500;background:none;cursor:auto;padding-right:24px;white-space:nowrap}
table.kv td{white-space:normal}
td.diff{background:rgba(217,98,43,.12)}
.stats{display:flex;flex-wrap:wrap;gap:10px;margin-bottom:12px} .stat{border:1px solid var(--line);border-radius:8px;padding:8px 14px;min-width:100px}
.stat b{display:block;font-size:22px;color:var(--acc)} .stat span{color:var(--muted);font-size:12px}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:0 24px;align-items:start}
details.blk{border:1px solid var(--line);border-radius:8px;margin:12px 0;background:var(--panel)}
details.blk>summary{cursor:pointer;padding:8px 12px;font-weight:600;color:var(--acc);list-style:none}
details.blk>summary::-webkit-details-marker{display:none}
details.blk>summary:before{content:"\\25B8";display:inline-block;width:16px;color:var(--muted)}
details.blk[open]>summary:before{content:"\\25BE"}
details.blk>summary small{color:var(--muted);font-weight:400;margin-left:6px}
details.blk>.cb{padding:0 12px 12px}
.desc{white-space:pre-wrap;color:var(--muted);font-size:13px;max-height:220px;overflow:auto;border:1px solid var(--line);border-radius:6px;padding:8px 12px}
svg.plot{width:100%;max-width:980px;height:auto;display:block;margin:8px 0;font-size:11px}
svg.plot text{fill:var(--muted)} svg.plot .lab{fill:var(--ink)} svg.plot .gl{stroke:var(--line)} svg.plot .ax{stroke:var(--muted)}
svg.plot .frame{fill:none;stroke:var(--line)}
canvas.zb{width:100%;max-width:820px;image-rendering:pixelated;border:1px solid var(--line);border-radius:4px;margin-top:8px}
ul.warns{margin:6px 0;padding-left:20px;color:var(--warn)}
@media (max-width:1000px){.wrap{grid-template-columns:1fr} nav{display:none}}
footer{color:var(--muted);font-size:12px;text-align:center;padding:10px}
@media print{nav,header{display:none} .wrap{display:block} details>*{display:block}}
</style></head>
<body>
<header><h1>XBeach · @@TITLE@@</h1><span>@@SUB@@</span><button class="theme" type="button" title="Switch between the light and dark theme">☾ Dark</button></header>
<div class="wrap"><nav><div class="tools"><button data-x="1">Expand all</button><button data-x="0">Collapse all</button></div>@@NAV@@</nav>
<main>@@BODY@@
<footer>Generated by modelvue.xbvue on @@GEN@@ from <code>@@PATH@@</code></footer></main></div>
<script>
document.querySelectorAll('.tools button').forEach(b=>b.addEventListener('click',()=>{
  const o=b.dataset.x==='1'; document.querySelectorAll('main details').forEach(d=>d.open=o);}));
document.querySelectorAll('table.sortable thead th').forEach((th,i)=>th.addEventListener('click',()=>{
  const tb=th.closest('table').tBodies[0], rows=[...tb.rows], asc=th.dataset.asc!=='1';
  rows.sort((a,b)=>{const x=a.cells[i].innerText,y=b.cells[i].innerText,nx=parseFloat(x),ny=parseFloat(y);
    const r=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y); return asc?r:-r;});
  th.dataset.asc=asc?'1':'0'; rows.forEach(r=>tb.appendChild(r));}));
// bathymetry heat maps (2D cases)
document.querySelectorAll('canvas.zb').forEach(cv=>{
  const z=JSON.parse(cv.dataset.z), lo=+cv.dataset.lo, hi=+cv.dataset.hi, ny=z.length, nx=z[0].length;
  cv.width=nx; cv.height=ny; const ctx=cv.getContext('2d'), im=ctx.createImageData(nx,ny);
  const mix=(a,b,t)=>a.map((v,k)=>Math.round(v+(b[k]-v)*t));
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const v=z[j][i];let c;
    if(v<0){c=mix([198,226,245],[12,52,110],Math.min(1,v/(lo||-1)));}
    else{c=mix([238,222,176],[120,84,40],Math.min(1,v/(hi||1)));}
    const p=4*(j*nx+i); im.data[p]=c[0];im.data[p+1]=c[1];im.data[p+2]=c[2];im.data[p+3]=255;}
  ctx.putImageData(im,0,0);});
const links=[...document.querySelectorAll('nav a')];
const obs=new IntersectionObserver(es=>es.forEach(e=>{if(e.isIntersecting){links.forEach(l=>l.classList.toggle('on',l.getAttribute('href')==='#'+e.target.id));}}),{rootMargin:'-40% 0px -55% 0px'});
document.querySelectorAll('main>section').forEach(s=>obs.observe(s));
const themeBtn=document.querySelector('header .theme');
const themeLabel=()=>{themeBtn.textContent=document.documentElement.dataset.theme==='dark'?'☀ Light':'☾ Dark';};
themeLabel();
themeBtn.addEventListener('click',()=>{const t=document.documentElement.dataset.theme==='dark'?'light':'dark';
  document.documentElement.dataset.theme=t; themeLabel(); try{localStorage.setItem('modelvue-theme',t)}catch(e){}});
</script></body></html>"""


def _kv_table(pairs) -> str:
    rows = "".join(f"<tr><th>{_e(k)}</th><td>{_e(v)}</td></tr>"
                   for k, v in pairs if v not in (None, "", [], {}))
    return f'<div class="tw"><table class="kv">{rows}</table></div>' if rows else ""


def _table(headers, rows, sortable=True) -> str:
    """rows contain pre-escaped HTML cells (or (html, css_class) tuples)."""
    if not rows:
        return ""
    th = "".join(f"<th>{_e(h)}</th>" for h in headers)
    tr = "".join("<tr>" + "".join(f'<td class="{c[1]}">{c[0]}</td>' if isinstance(c, tuple) else f"<td>{c}</td>"
                                  for c in r) + "</tr>" for r in rows)
    return (f'<div class="tw"><table class="{"sortable" if sortable else ""}">'
            f"<thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></div>")


def _blk(title, body, note="", open_=True) -> str:
    if not body:
        return ""
    n = f"<small>{_e(note)}</small>" if note else ""
    return f'<details class="blk"{" open" if open_ else ""}><summary>{title}{n}</summary><div class="cb">{body}</div></details>'


def _cid(c) -> str:
    return "case-" + re.sub(r"[^A-Za-z0-9_-]+", "-", c["name"])


def _case_section(c, many: bool) -> str:
    out = [f'<section id="{_cid(c)}"><h2>{_e(c["name"])}']
    if "error" in c:
        out.append(f'</h2><p class="warn">{_e(c["error"])}</p></section>')
        return "".join(out)
    gr, s = c["grid"], c["summary"]
    out.append(f' <span class="badge o">{"1D" if gr["ny"] == 0 else "2D"}</span>'
               f' <span class="badge o">{_e(s["Wave model"].split()[0])}</span>'
               f' <span class="badge o">morphology {"on" if s["Morphology"].startswith("ON") else "off"}</span>')
    if c["vegetation"]:
        out.append(' <span class="badge o">vegetation</span>')
    if not c["enabled"]:
        out.append(' <span class="badge w">disabled</span>')
    if c["warnings"]:
        out.append(f' <span class="badge w">{len(c["warnings"])} warning(s)</span>')
    out.append("</h2>")
    out.append('<p class="muted">' + _e(c["folder"]) + f' · <code>{_e(c["params_file"])}</code>'
               + (f' · written {_e(c["header"].get("date"))}' if c["header"].get("date") else "")
               + (f' by {_e(c["header"].get("function"))}' if c["header"].get("function") else "")
               + (f' · other params files: {_e(c["other_params_files"])}' if c["other_params_files"] else "")
               + "</p>")
    if c["warnings"]:
        out.append('<ul class="warns">' + "".join(f"<li>{_e(w)}</li>" for w in c["warnings"]) + "</ul>")

    items = list(s.items())
    half = (len(items) + 1) // 2
    out.append(_blk("Model summary", f'<div class="grid2"><div>{_kv_table(items[:half])}</div>'
                                     f'<div>{_kv_table(items[half:])}</div></div>'
                    + (_kv_table([("OET run settings", c["config_run"])]) if c.get("config_run") else "")))

    gk = [(k.replace("_", " "), v) for k, v in gr.items()
          if not k.startswith("_") and k not in ("zb_min", "zb_max", "nx", "ny")]
    hg = (len(gk) + 1) // 2
    out.append(_blk("Grid &amp; bathymetry", f'<div class="grid2"><div>{_kv_table(gk[:hg])}</div>'
                                              f'<div>{_kv_table(gk[hg:])}</div></div>'
                    + profile_svg(c) + map_canvas(c), open_=True))

    w = c["waves"]
    body = '<div class="grid2"><div>' + _kv_table([
        (w["type_key"] or "wave type", w["type"]), ("bcfile", w.get("bcfile")), ("format", w.get("format")),
        ("random", w.get("random")), ("rt", w.get("rt")), ("dtbc", w.get("dtbc"))]) + "</div><div>" + _kv_table([
        ("Hm0 range", w.get("Hm0_range")), ("Tp range", w.get("Tp_range")),
        ("direction(s)", w.get("direction")), ("total duration", w.get("total_duration")),
        ("spectrum", w.get("spectrum")), ("in params.txt", w.get("params_in_params_txt"))]) + "</div></div>"
    body += waves_svg(w)
    if w.get("rows"):
        body += _blk(f'Table rows ({len(w["rows"])})', _table(["#"] + w["columns"],
                     [[str(i + 1)] + [_e(g(v)) for v in r] for i, r in enumerate(w["rows"])]), open_=False)
    if w.get("filelist"):
        body += _table(["File", "Duration", "dtbc", "Exists"],
                       [[_e(f["file"]), _e(f["duration"]), _e(f["dtbc"]),
                         "yes" if f["exists"] else '<span class="warn">MISSING</span>'] for f in w["filelist"]])
    out.append(_blk("Wave boundary", body, w["type"]))

    wl = c["water_level"]
    body = _kv_table([(k, wl.get(k)) for k in ("mode", "zs0", "zs0file", "tidetype", "tideloc",
                                                "tidelen", "records", "time_span")])
    if wl.get("locations"):
        body += _table(["Column", "Location", "min zs", "max zs"],
                       [[str(L["column"]), _e(L["location"]), _e(L["min"]), _e(L["max"])] for L in wl["locations"]])
    body += tide_svg(wl)
    out.append(_blk("Water level", body, wl["mode"]))

    v = c["vegetation"]
    if v:
        keys = []
        for sp in v["species"]:
            keys += [k for k in sp if k not in keys and k != "file"]
        body = _kv_table([("veggiefile", v["veggiefile"]), ("veggiemapfile", v["veggiemapfile"]),
                          ("map cells", v.get("map_cells_per_species"))])
        body += _table(["Species file"] + keys, [[_e(sp["file"])] + [_e(sp.get(k, "")) for k in keys]
                                                 for sp in v["species"]])
        out.append(_blk("Vegetation", body))

    rows = [[_e(cat), f"<b>{_e(e['key'])}</b>", _e(e["value"]), f'<span class="muted">{_e(e["file_section"] or "")}</span>']
            for cat, items in c["parameters_by_category"].items() for e in items]
    out.append(_blk("Parameters", _table(["Category", "Key", "Value", "Section in file"], rows),
                    f'{c["n_parameters"]} set', open_=not many))

    o = c["output"]
    chips = lambda L: "".join(f"<span class=chip>{_e(x)}</span>" for x in L)
    body = _kv_table([(k, val) for k, val in o.items() if not isinstance(val, list)])
    vrows = [[lab, str(len(o[k])), chips(o[k])] for lab, k in
             (("global", "global_vars"), ("mean", "mean_vars"), ("point", "point_vars"),
              ("points (x y)", "points"), ("runup gauges", "runup_gauges")) if o[k]]
    body += _table(["Output", "n", "Variables / locations"], vrows, sortable=False)
    out.append(_blk("Output request", body, o["outputformat"], open_=not many))

    body = _table(["Key", "File", "Status"], [[_e(f["key"]), f'<code>{_e(f["file"])}</code>',
                  "ok" if f["exists"] else '<span class="warn">MISSING</span>'] for f in c["files"]])
    out.append(_blk("Referenced files", body, f'{len(c["files"])} file(s)', open_=not many))

    res = c["results"]
    out.append(_blk("Results", _kv_table(res.items()) if res else
                    '<p class="muted">No xboutput.nc / XBlog.txt found (case not run here yet).</p>', open_=not many))
    out.append("</section>")
    return "".join(out)


def render_html(cfg: dict) -> str:
    import datetime
    cases = cfg["cases"]
    many = len(cases) > 1
    oc = cfg["oet_config"]
    ok = [c for c in cases if "error" not in c]

    ov = []
    if many:
        n1d = sum(1 for c in ok if c["grid"]["ny"] == 0)
        st = [("Cases", len(cases)), ("Enabled", sum(c["enabled"] for c in cases)),
              ("1D / 2D", f"{n1d} / {len(ok) - n1d}"),
              ("Differing params", len(cfg["differing_parameters"])),
              ("With results", sum(1 for c in ok if c["results"])),
              ("Warnings", sum(len(c["warnings"]) for c in ok))]
        ov += ['<section id="overview"><h2>Overview</h2><div class="stats">',
               "".join(f'<div class="stat"><b>{_e(v)}</b><span>{_e(k)}</span></div>' for k, v in st), "</div>"]
        ov.append(_kv_table([("Folder", cfg["root"]),
                             ("OET .config", oc and oc["file"]),
                             ("Responsible", oc and oc["general"].get("responsible")),
                             ("Analysis", oc and oc["general"].get("analysis")),
                             ("Categories", oc and oc["categories"])]))
        if cfg["readme"]:
            ov.append(f'<h4>README</h4><div class="desc">{_e(cfg["readme"])}</div>')
        diff = cfg["differing_parameters"]
        rows = []
        for c in cases:
            r = case_row(c)
            r[0] = f'<a href="#{_cid(c)}">{_e(r[0])}</a>'
            rows.append([r[0]] + [_e(x) for x in r[1:]])
        ov.append("<h4>Cases</h4>" + _table(["Case", "Params", "nx x ny", "tstop", "Waves", "Morphology",
                                              "zb range (m)", "Results"], rows))
        if diff:
            ov.append(f"<h4>Parameters that differ between cases ({len(diff)})</h4>" + _table(
                ["Case"] + diff, [[f'<a href="#{_cid(c)}">{_e(c["name"])}</a>']
                                  + [(_e(param_ci(c, k) or "—"), "diff") for k in diff] for c in ok]))
        else:
            ov.append('<p class="muted">All cases use identical parameter values.</p>')
        profs = [c for c in ok if c["grid"].get("_plot", {}).get("profile")]
        if len(profs) > 1:
            ov.append(f"<h4>Bed profiles ({len(profs)} cases)</h4>" + svg_plot(
                [{"x": c["grid"]["_plot"]["profile"][0], "y": c["grid"]["_plot"]["profile"][1],
                  "label": c["name"], "width": 1.2} for c in profs], "x (m)", "zb (m)", h=300))
        ov.append("</section>")

    nav = (['<a href="#overview">Overview</a>'] if many else []) + ['<div class="grp">Cases</div>']
    for c in cases:
        sub = "error" if "error" in c else f'{"1D" if c["grid"]["ny"] == 0 else "2D"} · {c["waves"]["type"]}'
        nav.append(f'<a href="#{_cid(c)}"{"" if c["enabled"] else " class=off"}>{_e(c["name"])}'
                   f'<br><small>{_e(sub)}</small></a>')

    body = "".join(ov) + "".join(_case_section(c, many) for c in cases)
    if many:
        sub = f'{len(cases)} case(s)' + (' · OET .config' if oc else "")
    else:
        c = cases[0]
        sub = c.get("summary", {}).get("Dimension", "") + " · " + c.get("summary", {}).get("Simulation time (tstop)", "")
    rep = {"@@TITLE@@": _e(cfg["name"]), "@@SUB@@": _e(sub), "@@NAV@@": "".join(nav), "@@BODY@@": body,
           "@@GEN@@": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), "@@PATH@@": _e(cfg["root"])}
    out = HTML_TEMPLATE
    for k, v in rep.items():
        out = out.replace(k, v)
    return out


def write_html(cfg: dict, out: Path, open_browser: bool = True) -> Path:
    out.write_text(render_html(cfg), encoding="utf-8")
    print(f"Webpage written to {out.resolve()}")
    if open_browser:
        import webbrowser
        webbrowser.open(out.resolve().as_uri())
    return out


# --------------------------------------------------------------------------- #
def main(argv=None, prog=None):
    ap = argparse.ArgumentParser(prog=prog, description="Show the configuration of an XBeach model (or case set).")
    ap.add_argument("model", help="params.txt, a case folder, or a folder of cases (.config / sub-folders)")
    ap.add_argument("--case", help="show only this case of a case set, e.g. Raai1085")
    ap.add_argument("--params", help="params file name to use instead of params.txt / the .config entry")
    ap.add_argument("--details", action="store_true", help="terminal: full detail for every case in a set")
    ap.add_argument("--json", metavar="FILE", help="also write the configuration to a JSON file")
    ap.add_argument("--html", nargs="?", const="", metavar="FILE",
                    help="write the configuration as a webpage and open it in the browser, "
                         "instead of printing it in the terminal "
                         "(default file: <folder>_config.html in the model folder)")
    ap.add_argument("--no-browser", action="store_true", help="with --html, write the page but don't open it")
    a = ap.parse_args(argv)

    cfg = build_config(a.model, a.params, a.case)
    if a.html is None:
        show(cfg, a.details)
    else:
        name = cfg["name"] + (f"_{a.case}" if a.case else "")
        out = Path(a.html) if a.html else Path(cfg["root"]) / f"{name}_config.html"
        write_html(cfg, out, open_browser=not a.no_browser)
    if a.json:
        Path(a.json).write_text(json.dumps(strip_private(cfg), indent=2, default=str), encoding="utf-8")
        print(f"Configuration written to {a.json}")


if __name__ == "__main__":
    main()
