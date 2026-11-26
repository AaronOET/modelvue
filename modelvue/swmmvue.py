"""
modelvue.swmmvue - Read an EPA SWMM 5 model and show its configuration.

Parses the input file (.inp) section by section following the SWMM 5.2 input
format (Users Manual, Appendix D / src/solver/keywords.c), checks the network,
and - if the model has been run - reads the report file (.rpt: continuity
errors, warnings, flooding, instability) and the binary results file (.out,
layout as written by src/solver/output.c) for peak values and time series.

No extra packages needed (standard library only).

Usage
-----
    swmmvue  model.inp                    # report in the terminal
    swmmvue  folder                       # every .inp in the folder (and below)
    swmmvue  model.inp --all              # every row of every section
    swmmvue  model.inp --rpt run.rpt --out run.out   # results with other names
    swmmvue  results.out                  # only a binary results file
    swmmvue  model.inp --json config.json
    swmmvue  model.inp --html [report.html] [--no-browser]
"""
from __future__ import annotations

import argparse
import array
import datetime as _dt
import html as _h
import json
import math
import re
import struct
import sys
from pathlib import Path

from ._common import _blk, _e, h1, h2, _kv, strip_private, _table, _ticks, _two
from ._inp import num, table

# --------------------------------------------------------------------------- #
# Section definitions (SWMM 5.2)
# --------------------------------------------------------------------------- #
SECTION_ORDER = [
    "TITLE", "OPTIONS", "FILES", "RAINGAGES", "TEMPERATURE", "EVAPORATION", "ADJUSTMENTS", "EVENTS",
    "SUBCATCHMENTS", "SUBAREAS", "INFILTRATION", "LID_CONTROLS", "LID_USAGE", "AQUIFERS", "GROUNDWATER", "GWF",
    "SNOWPACKS", "JUNCTIONS", "OUTFALLS", "DIVIDERS", "STORAGE", "CONDUITS", "PUMPS", "ORIFICES", "WEIRS",
    "OUTLETS", "XSECTIONS", "TRANSECTS", "STREETS", "INLETS", "INLET_USAGE", "LOSSES", "CONTROLS", "POLLUTANTS",
    "LANDUSES", "COVERAGES", "LOADINGS", "BUILDUP", "WASHOFF", "TREATMENT", "INFLOWS", "DWF", "RDII",
    "HYDROGRAPHS", "CURVES", "TIMESERIES", "PATTERNS", "REPORT", "TAGS", "MAP", "COORDINATES", "VERTICES",
    "POLYGONS", "SYMBOLS", "LABELS", "BACKDROP", "PROFILES"]

# accepted abbreviations (SWMM matches on the leading characters, keywords.c ws_*)
SECTION_ALIASES = {"JUNC": "JUNCTIONS", "XSECT": "XSECTIONS", "INFIL": "INFILTRATION", "EVAP": "EVAPORATION",
                   "POLLUT": "POLLUTANTS", "SNOWMELT": "SNOWPACKS"}

RAW_SECTIONS = {"TITLE", "CONTROLS", "TRANSECTS", "LABELS", "BACKDROP", "FILES", "TEMPERATURE", "TAGS"}

COLS = {
    "RAINGAGES": ["Name", "Format", "Interval", "SCF", "Source", "Data", "Station", "Units"],
    "SUBCATCHMENTS": ["Name", "RainGage", "Outlet", "Area", "%Imperv", "Width", "%Slope", "CurbLen", "SnowPack"],
    "SUBAREAS": ["Subcatch", "N-Imperv", "N-Perv", "S-Imperv", "S-Perv", "%Zero", "RouteTo", "%Routed"],
    "JUNCTIONS": ["Name", "Elevation", "MaxDepth", "InitDepth", "SurDepth", "Aponded"],
    "OUTFALLS": ["Name", "Elevation", "Type", "Stage/Data", "Gated", "RouteTo"],
    "STORAGE": ["Name", "Elevation", "MaxDepth", "InitDepth", "Shape", "Curve / params", "SurDepth", "Fevap", "Psi", "Ksat", "IMD"],
    "DIVIDERS": ["Name", "Elevation", "DivLink", "Type", "Params"],
    "CONDUITS": ["Name", "From", "To", "Length", "Roughness", "InOffset", "OutOffset", "InitFlow", "MaxFlow"],
    "PUMPS": ["Name", "From", "To", "Curve", "Status", "Startup", "Shutoff"],
    "ORIFICES": ["Name", "From", "To", "Type", "Offset", "Qcoeff", "Gated", "CloseTime"],
    "WEIRS": ["Name", "From", "To", "Type", "CrestHt", "Qcoeff", "Gated", "EndCon", "EndCoeff", "Surcharge",
              "RoadWidth", "RoadSurf", "CoeffCurve"],
    "OUTLETS": ["Name", "From", "To", "Offset", "Type", "Qtable/Qcoeff", "Qexpon", "Gated"],
    "XSECTIONS": ["Link", "Shape", "Geom1", "Geom2", "Geom3", "Geom4", "Barrels", "Culvert"],
    "LOSSES": ["Link", "Kentry", "Kexit", "Kavg", "FlapGate", "Seepage"],
    "POLLUTANTS": ["Name", "Units", "Crain", "Cgw", "Crdii", "Kdecay", "SnowOnly", "CoPollutant", "CoFrac", "Cdwf", "Cinit"],
    "LANDUSES": ["Name", "SweepInterval", "Availability", "LastSweep"],
    "COVERAGES": ["Subcatch", "LandUse", "Percent"],
    "LOADINGS": ["Subcatch", "Pollutant", "Buildup"],
    "BUILDUP": ["LandUse", "Pollutant", "Function", "Coeff1", "Coeff2", "Coeff3", "PerUnit"],
    "WASHOFF": ["LandUse", "Pollutant", "Function", "Coeff1", "Coeff2", "SweepRmvl", "BmpRmvl"],
    "TREATMENT": ["Node", "Pollutant", "Function"],
    "INFLOWS": ["Node", "Constituent", "TimeSeries", "Type", "Mfactor", "Sfactor", "Baseline", "Pattern"],
    "DWF": ["Node", "Constituent", "Baseline", "Pattern1", "Pattern2", "Pattern3", "Pattern4"],
    "RDII": ["Node", "UnitHydrograph", "SewerArea"],
    "HYDROGRAPHS": ["Name", "Month/RainGage", "Response", "R", "T", "K", "Dmax", "Drecov", "Dinit"],
    "LID_USAGE": ["Subcatch", "LID", "Number", "Area", "Width", "InitSat", "FromImp", "ToPerv", "RptFile", "DrainTo", "FromPerv"],
    "AQUIFERS": ["Name", "Por", "WP", "FC", "Ksat", "Kslope", "Tslope", "ETu", "ETs", "Seep", "Ebot", "Egw", "Umc", "ETupat"],
    "GROUNDWATER": ["Subcatch", "Aquifer", "Node", "Esurf", "A1", "B1", "A2", "B2", "A3", "Dsw", "Egwt", "Ebot", "Wgr", "Umc"],
    "INLET_USAGE": ["Conduit", "Inlet", "Node", "Number", "%Clogged", "Qmax", "aLocal", "wLocal", "Placement"],
    "STREETS": ["Name", "Tcrown", "Hcurb", "Sx", "nRoad", "a", "W", "Sides", "Tback", "Sback", "nBack"],
    "COORDINATES": ["Node", "X", "Y"], "VERTICES": ["Link", "X", "Y"], "POLYGONS": ["Subcatch", "X", "Y"],
    "SYMBOLS": ["Gage", "X", "Y"],
}
INFIL_COLS = {"HORTON": ["Subcatch", "MaxRate", "MinRate", "Decay", "DryTime", "MaxInfil"],
              "MODIFIED_HORTON": ["Subcatch", "MaxRate", "MinRate", "Decay", "DryTime", "MaxInfil"],
              "GREEN_AMPT": ["Subcatch", "Suction", "Ksat", "IMD"],
              "MODIFIED_GREEN_AMPT": ["Subcatch", "Suction", "Ksat", "IMD"],
              "CURVE_NUMBER": ["Subcatch", "CurveNum", "(unused)", "DryTime"]}

NODE_SECTIONS = ("JUNCTIONS", "OUTFALLS", "DIVIDERS", "STORAGE")
LINK_SECTIONS = ("CONDUITS", "PUMPS", "ORIFICES", "WEIRS", "OUTLETS")

OPTION_HELP = {
    "FLOW_UNITS": "flow units (CFS, GPM, MGD = US; CMS, LPS, MLD = SI)", "INFILTRATION": "infiltration model",
    "FLOW_ROUTING": "routing method (STEADY, KINWAVE, DYNWAVE)", "LINK_OFFSETS": "offsets given as DEPTH or ELEVATION",
    "FORCE_MAIN_EQUATION": "force-main friction (H-W or D-W)", "IGNORE_RAINFALL": "skip rainfall/runoff",
    "IGNORE_SNOWMELT": "skip snowmelt", "IGNORE_GROUNDWATER": "skip groundwater", "IGNORE_RDII": "skip RDII",
    "IGNORE_ROUTING": "skip flow routing", "IGNORE_QUALITY": "skip water quality", "ALLOW_PONDING": "ponding at flooded nodes",
    "SKIP_STEADY_STATE": "skip steady periods", "SYS_FLOW_TOL": "steady-state system flow tolerance (%)",
    "LAT_FLOW_TOL": "steady-state lateral flow tolerance (%)", "START_DATE": "simulation start", "START_TIME": "",
    "END_DATE": "simulation end", "END_TIME": "", "REPORT_START_DATE": "reporting start", "REPORT_START_TIME": "",
    "SWEEP_START": "street sweeping start", "SWEEP_END": "street sweeping end", "DRY_DAYS": "antecedent dry days",
    "REPORT_STEP": "reporting time step", "WET_STEP": "runoff step, wet periods", "DRY_STEP": "runoff step, dry periods",
    "ROUTING_STEP": "routing time step", "RULE_STEP": "control-rule evaluation step", "INERTIAL_DAMPING": "DYNWAVE inertial damping",
    "NORMAL_FLOW_LIMITED": "DYNWAVE normal-flow limitation", "VARIABLE_STEP": "DYNWAVE variable step safety factor (0 = fixed)",
    "LENGTHENING_STEP": "DYNWAVE conduit lengthening step (s)", "MIN_SURFAREA": "DYNWAVE minimum nodal surface area",
    "MAX_TRIALS": "DYNWAVE trials per step", "HEAD_TOLERANCE": "DYNWAVE head convergence tolerance",
    "MINIMUM_STEP": "DYNWAVE minimum variable step (s)", "THREADS": "parallel threads", "MIN_SLOPE": "minimum conduit slope (%)",
    "SURCHARGE_METHOD": "EXTRAN or SLOT", "TEMPDIR": "folder for scratch files", "COMPATIBILITY": "", "SLOPE_WEIGHTING": "",
}

FLOW_UNITS = ["CFS", "GPM", "MGD", "CMS", "LPS", "MLD"]


def units_for(flow_units: str) -> dict:
    us = flow_units.upper() in ("CFS", "GPM", "MGD")
    return {"system": "US" if us else "SI", "flow": flow_units.upper(), "length": "ft" if us else "m",
            "area": "ac" if us else "ha", "rain": "in" if us else "mm", "rate": "in/hr" if us else "mm/hr",
            "volume": "ft3" if us else "m3"}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def g(v, nd=4):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "" if v is None else str(v)
    if math.isnan(f):
        return "nan"
    if f.is_integer() and abs(f) < 1e12:
        return str(int(f))
    return f"{f:.{nd}g}"


def tokenize(line: str) -> list[str]:
    """Whitespace split; "quoted names" kept whole, "" kept as an empty token, ; starts a comment."""
    toks, cur, inq, quoted = [], "", False, False
    for ch in line:
        if ch == '"':
            inq = not inq
            quoted = True
            continue
        if ch == ";" and not inq:
            break
        if ch.isspace() and not inq:
            if cur or quoted:
                toks.append(cur)
            cur, quoted = "", False
            continue
        cur += ch
    if cur or quoted:
        toks.append(cur)
    return toks


class CI(dict):
    """dict with case-insensitive keys - SWMM object names are not case sensitive (hash.c: samestr)."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self._u = {str(x).upper(): x for x in dict.keys(self)}

    def _k(self, k):
        return self._u.get(str(k).upper(), k)

    def __setitem__(self, k, v):
        super().__setitem__(self._k(k) if str(k).upper() in self._u else k, v)
        self._u.setdefault(str(k).upper(), k)

    def __getitem__(self, k):
        return super().__getitem__(self._k(k))

    def __contains__(self, k):
        return str(k).upper() in self._u

    def get(self, k, d=None):
        return super().get(self._k(k), d)

    def setdefault(self, k, d=None):
        if k not in self:
            self[k] = d
        return self[k]


def rng(vals, unit=""):
    v = [x for x in vals if x is not None]
    if not v:
        return "-"
    lo, hi = min(v), max(v)
    return f"{g(lo)}{unit}" if lo == hi else f"{g(lo)} .. {g(hi)}{unit}"


def parse_hours(t: str) -> float | None:
    """'1:30', '1:30:00', '1.5' -> hours."""
    if t is None:
        return None
    if ":" in t:
        p = [num(x, 0) for x in t.split(":")]
        return p[0] + (p[1] if len(p) > 1 else 0) / 60 + (p[2] if len(p) > 2 else 0) / 3600
    return num(t)


def parse_date(d: str):
    for fmt in ("%m/%d/%Y", "%m/%d/%y", "%Y-%m-%d"):
        try:
            return _dt.datetime.strptime(d, fmt)
        except (ValueError, TypeError):
            pass
    return None


def swmm_datetime(days: float) -> _dt.datetime:
    """SWMM DateTime: decimal days since 12/30/1899 (datetime.c, DateDelta = 693594)."""
    return _dt.datetime(1899, 12, 30) + _dt.timedelta(days=days)


# --------------------------------------------------------------------------- #
# INP file
# --------------------------------------------------------------------------- #
def read_inp(path: Path) -> dict:
    secs, raw, order, cur = {}, {}, [], None
    for ln in path.read_text(errors="replace").splitlines():
        s = ln.strip()
        if s.startswith("[") and "]" in s:
            name = s[1:s.index("]")].strip().upper()
            name = SECTION_ALIASES.get(name, name)
            for full in SECTION_ORDER:          # SWMM accepts leading-substring matches
                if full.startswith(name) or name.startswith(full):
                    name = full
                    break
            cur = name
            secs.setdefault(cur, [])
            raw.setdefault(cur, [])
            if cur not in order:
                order.append(cur)
            continue
        if cur is None:
            continue
        if cur in RAW_SECTIONS:
            if s and not s.startswith(";;"):
                raw[cur].append(ln.rstrip())
            continue
        t = tokenize(ln)
        if t:
            secs[cur].append(t)
    return {"sections": secs, "raw": raw, "order": order}


def options(inp) -> dict:
    out = {}
    for t in table(inp, "OPTIONS"):
        out[t[0].upper()] = " ".join(t[1:])
    return out


def rows_as_dicts(rows, cols):
    out = []
    for t in rows:
        d = {c: (t[i] if i < len(t) else "") for i, c in enumerate(cols)}
        if len(t) > len(cols):
            d["(extra)"] = " ".join(t[len(cols):])
        out.append(d)
    return out


def parse_outfalls(rows):
    out = []
    for t in rows:
        typ = t[2].upper() if len(t) > 2 else ""
        if typ in ("FREE", "NORMAL"):
            d = {"Name": t[0], "Elevation": t[1], "Type": typ, "Stage/Data": "", "Gated": t[3] if len(t) > 3 else "",
                 "RouteTo": t[4] if len(t) > 4 else ""}
        else:
            d = {"Name": t[0], "Elevation": t[1] if len(t) > 1 else "", "Type": typ, "Stage/Data": t[3] if len(t) > 3 else "",
                 "Gated": t[4] if len(t) > 4 else "", "RouteTo": t[5] if len(t) > 5 else ""}
        out.append(d)
    return out


STORAGE_NPAR = {"FUNCTIONAL": 3, "TABULAR": 1, "CYLINDRICAL": 3, "CONICAL": 3, "PARABOLIC": 3, "PYRAMIDAL": 3}


def parse_storage(rows):
    out = []
    for t in rows:
        shape = t[4].upper() if len(t) > 4 else ""
        n = STORAGE_NPAR.get(shape, 1)
        rest = t[5 + n:]
        d = {"Name": t[0], "Elevation": t[1] if len(t) > 1 else "", "MaxDepth": t[2] if len(t) > 2 else "",
             "InitDepth": t[3] if len(t) > 3 else "", "Shape": shape, "Curve / params": " ".join(t[5:5 + n]),
             "SurDepth": rest[0] if rest else "", "Fevap": rest[1] if len(rest) > 1 else "",
             "Psi": rest[2] if len(rest) > 2 else "", "Ksat": rest[3] if len(rest) > 3 else "", "IMD": rest[4] if len(rest) > 4 else ""}
        out.append(d)
    return out


def parse_curves(rows) -> dict:
    curves = CI()
    for t in rows:
        name = t[0]
        vals = t[1:]
        if name not in curves:
            typ = vals[0].upper() if vals and num(vals[0]) is None else ""
            if typ:
                vals = vals[1:]
            curves[name] = {"type": typ, "points": []}
        for i in range(0, len(vals) - 1, 2):
            x, y = num(vals[i]), num(vals[i + 1])
            if x is not None and y is not None:
                curves[name]["points"].append((x, y))

    return curves


def parse_timeseries(rows) -> dict:
    ts = CI()
    for t in rows:
        name = t[0]
        e = ts.setdefault(name, {"points": [], "file": None, "dates": False})
        if len(t) >= 3 and t[1].upper() == "FILE":
            e["file"] = t[2]
            continue
        vals, i, date = t[1:], 0, None
        while i < len(vals):
            tok = vals[i]
            if ("/" in tok or "-" in tok[1:]) and parse_date(tok):
                date = parse_date(tok)
                e["dates"] = True
                i += 1
                continue
            if i + 1 >= len(vals):
                break
            h = parse_hours(tok)
            v = num(vals[i + 1])
            if h is not None and v is not None:
                if date:
                    e["points"].append((date + _dt.timedelta(hours=h), v))
                else:
                    e["points"].append((h, v))
            i += 2
    for e in ts.values():
        p = e["points"]
        if p and isinstance(p[0][0], _dt.datetime):
            t0 = p[0][0]
            e["start"] = t0.strftime("%Y-%m-%d %H:%M")
            e["points"] = [((d - t0).total_seconds() / 3600, v) for d, v in p]
    return ts


def parse_patterns(rows) -> dict:
    pats = CI()
    for t in rows:
        name = t[0]
        vals = t[1:]
        if name not in pats:
            typ = vals[0].upper() if vals and num(vals[0]) is None else ""
            pats[name] = {"type": typ, "factors": []}
            if typ:
                vals = vals[1:]
        pats[name]["factors"] += [num(v) for v in vals if num(v) is not None]
    return pats


def parse_lid_controls(rows) -> dict:
    lids = {}
    for t in rows:
        name = t[0]
        if len(t) == 2:
            lids[name] = {"type": t[1].upper(), "layers": {}}
        elif name in lids and len(t) > 1:
            lids[name]["layers"][t[1].upper()] = " ".join(t[2:])
    return lids


def parse_profiles(rows) -> dict:
    prof = {}
    for t in rows:
        prof.setdefault(t[0], []).extend(t[1:])
    return prof


LID_TYPES = {"BC": "bio-retention cell", "RG": "rain garden", "GR": "green roof", "IT": "infiltration trench",
             "PP": "permeable pavement", "RB": "rain barrel", "RD": "rooftop disconnection", "VS": "vegetative swale"}


def build_model(inp: dict) -> dict:
    opt = options(inp)
    fu = opt.get("FLOW_UNITS", "CFS").upper()
    U = units_for(fu)
    infil_model = opt.get("INFILTRATION", "HORTON").upper()
    m = {"options": opt, "units": U}
    m["title"] = [ln.strip() for ln in inp["raw"].get("TITLE", []) if ln.strip()]
    m["raingages"] = rows_as_dicts(table(inp, "RAINGAGES"), COLS["RAINGAGES"])
    m["subcatchments"] = rows_as_dicts(table(inp, "SUBCATCHMENTS"), COLS["SUBCATCHMENTS"])
    m["subareas"] = CI({r["Subcatch"]: r for r in rows_as_dicts(table(inp, "SUBAREAS"), COLS["SUBAREAS"])})
    inf = CI()
    for t in table(inp, "INFILTRATION"):
        method = infil_model
        if t[-1].upper() in INFIL_COLS:          # SWMM 5.2: optional per-subcatchment method at the end
            method, t = t[-1].upper(), t[:-1]
        cols = INFIL_COLS.get(method, ["Subcatch"] + [f"p{i}" for i in range(1, 6)])
        d = {c: (t[i] if i < len(t) else "") for i, c in enumerate(cols)}
        d["method"] = method
        inf[t[0]] = d
    m["infiltration"] = inf
    m["junctions"] = rows_as_dicts(table(inp, "JUNCTIONS"), COLS["JUNCTIONS"])
    m["outfalls"] = parse_outfalls(table(inp, "OUTFALLS"))
    m["storage"] = parse_storage(table(inp, "STORAGE"))
    m["dividers"] = rows_as_dicts(table(inp, "DIVIDERS"), COLS["DIVIDERS"])
    for sec in LINK_SECTIONS:
        m[sec.lower()] = rows_as_dicts(table(inp, sec), COLS[sec])
    m["xsections"] = CI({t[0]: {"Shape": t[1].upper() if len(t) > 1 else "", "geom": t[2:]} for t in table(inp, "XSECTIONS")})
    m["losses"] = rows_as_dicts(table(inp, "LOSSES"), COLS["LOSSES"])
    m["curves"] = parse_curves(table(inp, "CURVES"))
    m["timeseries"] = parse_timeseries(table(inp, "TIMESERIES"))
    m["patterns"] = parse_patterns(table(inp, "PATTERNS"))
    m["pollutants"] = rows_as_dicts(table(inp, "POLLUTANTS"), COLS["POLLUTANTS"])
    m["landuses"] = rows_as_dicts(table(inp, "LANDUSES"), COLS["LANDUSES"])
    m["inflows"] = rows_as_dicts(table(inp, "INFLOWS"), COLS["INFLOWS"])
    m["dwf"] = rows_as_dicts(table(inp, "DWF"), COLS["DWF"])
    m["rdii"] = rows_as_dicts(table(inp, "RDII"), COLS["RDII"])
    m["lid_controls"] = parse_lid_controls(table(inp, "LID_CONTROLS"))
    m["lid_usage"] = rows_as_dicts(table(inp, "LID_USAGE"), COLS["LID_USAGE"])
    m["controls"] = inp["raw"].get("CONTROLS", [])
    m["rules"] = [ln.split(None, 1)[1] if len(ln.split(None, 1)) > 1 else "" for ln in m["controls"]
                  if ln.strip().upper().startswith("RULE")]
    m["evaporation"] = {t[0].upper(): " ".join(t[1:]) for t in table(inp, "EVAPORATION")}
    m["report"] = {t[0].upper(): " ".join(t[1:]) for t in table(inp, "REPORT")}
    m["profiles"] = parse_profiles(table(inp, "PROFILES"))
    m["coords"] = CI({t[0]: (num(t[1]), num(t[2])) for t in table(inp, "COORDINATES") if len(t) >= 3})
    vert = CI()
    for t in table(inp, "VERTICES"):
        if len(t) >= 3:
            vert.setdefault(t[0], []).append((num(t[1]), num(t[2])))
    m["vertices"] = vert
    poly = CI()
    for t in table(inp, "POLYGONS"):
        if len(t) >= 3:
            poly.setdefault(t[0], []).append((num(t[1]), num(t[2])))
    m["polygons"] = poly
    m["symbols"] = {t[0]: (num(t[1]), num(t[2])) for t in table(inp, "SYMBOLS") if len(t) >= 3}
    m["map"] = {t[0].upper(): " ".join(t[1:]) for t in table(inp, "MAP")}
    return m


# --------------------------------------------------------------------------- #
# Derived information & checks
# --------------------------------------------------------------------------- #
def node_table(m) -> dict:
    nodes = CI()
    for typ, rows in (("JUNCTION", m["junctions"]), ("OUTFALL", m["outfalls"]), ("STORAGE", m["storage"]),
                      ("DIVIDER", m["dividers"])):
        for r in rows:
            nodes[r["Name"]] = {"type": typ, "invert": num(r.get("Elevation")),
                                "maxdepth": num(r.get("MaxDepth"), 0.0) if typ != "OUTFALL" else 0.0}
    return nodes


def link_table(m) -> dict:
    links = CI()
    for sec, typ in (("conduits", "CONDUIT"), ("pumps", "PUMP"), ("orifices", "ORIFICE"), ("weirs", "WEIR"),
                     ("outlets", "OUTLET")):
        for r in m[sec]:
            links[r["Name"]] = {"type": typ, "from": r["From"], "to": r["To"], "row": r}
    return links


def conduit_geometry(m, nodes) -> list[dict]:
    offs_elev = m["options"].get("LINK_OFFSETS", "DEPTH").upper() == "ELEVATION"
    out = []
    for r in m["conduits"]:
        n1, n2 = nodes.get(r["From"]), nodes.get(r["To"])
        L = num(r["Length"])
        o1, o2 = num(r["InOffset"], 0.0), num(r["OutOffset"], 0.0)
        x = m["xsections"].get(r["Name"], {})
        g1 = num(x.get("geom", [None])[0]) if x.get("geom") else None
        e = {"name": r["Name"], "from": r["From"], "to": r["To"], "length": L, "n": num(r["Roughness"]),
             "shape": x.get("Shape", ""), "geom1": g1, "barrels": num(x["geom"][4], 1) if len(x.get("geom", [])) > 4 else 1}
        if n1 and n2 and n1["invert"] is not None and n2["invert"] is not None and L:
            z1 = o1 if offs_elev else n1["invert"] + o1
            z2 = o2 if offs_elev else n2["invert"] + o2
            e["z_up"], e["z_dn"] = z1, z2
            e["slope_pct"] = 100 * (z1 - z2) / L
        out.append(e)
    return out


def rain_summary(m) -> list[dict]:
    out = []
    for r in m["raingages"]:
        e = {"name": r["Name"], "format": r["Format"].upper(), "interval": r["Interval"], "scf": r["SCF"],
             "source": r["Source"].upper(), "data": r["Data"]}
        if e["source"] == "TIMESERIES":
            ts = m["timeseries"].get(r["Data"])
            if ts is None:
                e["problem"] = f"time series '{r['Data']}' not defined"
            elif ts["file"]:
                e["note"] = f"time series read from file {ts['file']}"
            elif ts["points"]:
                p = ts["points"]
                dt = parse_hours(r["Interval"]) or ((p[1][0] - p[0][0]) if len(p) > 1 else 1)
                vals = [v for _, v in p]
                if e["format"] == "INTENSITY":
                    # each value holds until the next time (the last one for one interval)
                    depth = sum(p[i][1] * (p[i + 1][0] - p[i][0]) for i in range(len(p) - 1)) + p[-1][1] * dt
                    e["peak"] = f"{g(max(vals))} {m['units']['rate']}"
                elif e["format"] == "VOLUME":
                    depth = sum(vals)
                    e["peak"] = f"{g(max(vals) / dt)} {m['units']['rate']} (max interval depth / interval)"
                else:  # CUMULATIVE
                    depth = vals[-1]
                e["total_depth"] = f"{g(depth * num(r['SCF'], 1.0))} {m['units']['rain']}"
                e["duration"] = f"{g(p[-1][0] - p[0][0])} h"
                if ts.get("start"):
                    e["starts"] = ts["start"]
                e["_series"] = (p, dt, e["format"])
        elif e["source"] == "FILE":
            e["note"] = f"rainfall file {r['Data']} station {r['Station']} ({r['Units']})"
        out.append(e)
    return out


def checks(m, nodes, links, cond) -> list[str]:
    W = []
    names = CI()
    for sec, rows in (("JUNCTIONS", m["junctions"]), ("OUTFALLS", m["outfalls"]), ("STORAGE", m["storage"]),
                      ("DIVIDERS", m["dividers"])):
        for r in rows:
            names.setdefault(r["Name"], []).append(sec)
    dup = [n for n, v in names.items() if len(v) > 1]
    if dup:
        W.append(f"duplicate node names: {', '.join(dup[:10])}")
    for n, l in links.items():
        for end in ("from", "to"):
            if l[end] not in nodes:
                W.append(f"{l['type'].lower()} {n}: {end} node '{l[end]}' is not defined")
    sub_names = {r["Name"].upper() for r in m["subcatchments"]}
    gages = {r["Name"].upper() for r in m["raingages"]}
    for r in m["subcatchments"]:
        if r["Outlet"] not in nodes and r["Outlet"].upper() not in sub_names:
            W.append(f"subcatchment {r['Name']}: outlet '{r['Outlet']}' is neither a node nor a subcatchment")
        if r["RainGage"].upper() not in gages:
            W.append(f"subcatchment {r['Name']}: rain gage '{r['RainGage']}' is not defined")
        if r["Name"] not in m["subareas"]:
            W.append(f"subcatchment {r['Name']}: no [SUBAREAS] entry")
        if m["infiltration"] and r["Name"] not in m["infiltration"]:
            W.append(f"subcatchment {r['Name']}: no [INFILTRATION] entry")
    for c in m["conduits"]:
        if c["Name"] not in m["xsections"]:
            W.append(f"conduit {c['Name']}: no [XSECTIONS] entry")
    for o in m["orifices"] + m["weirs"]:
        if o["Name"] not in m["xsections"]:
            W.append(f"{'orifice' if o in m['orifices'] else 'weir'} {o['Name']}: no [XSECTIONS] entry")
    used = {l["from"].upper() for l in links.values()} | {l["to"].upper() for l in links.values()}
    orphan = [n for n in nodes if n.upper() not in used]
    if orphan:
        W.append(f"{len(orphan)} node(s) not connected to any link: {', '.join(orphan[:10])}{' ...' if len(orphan) > 10 else ''}")
    if nodes and not m["outfalls"]:
        W.append("no outfall in the network")
    # referenced curves / time series / patterns
    ts_ref, cv_ref, pat_ref = set(), set(), set()
    for r in m["raingages"]:
        if r["Source"].upper() == "TIMESERIES":
            ts_ref.add(r["Data"])
    for r in m["inflows"]:
        if r["TimeSeries"]:
            ts_ref.add(r["TimeSeries"])
        if r["Pattern"]:
            pat_ref.add(r["Pattern"])
    for r in m["dwf"]:
        for k in ("Pattern1", "Pattern2", "Pattern3", "Pattern4"):
            if r[k]:
                pat_ref.add(r[k])
    for r in m["outfalls"]:
        if r["Type"] == "TIDAL":
            cv_ref.add(r["Stage/Data"])
        if r["Type"] == "TIMESERIES":
            ts_ref.add(r["Stage/Data"])
    for r in m["storage"]:
        if r["Shape"] == "TABULAR":
            cv_ref.add(r["Curve / params"])
    for r in m["pumps"]:
        if r["Curve"] != "*":
            cv_ref.add(r["Curve"])
    for n, x in m["xsections"].items():
        if x["Shape"] in ("CUSTOM",) and len(x["geom"]) > 1:
            cv_ref.add(x["geom"][1])
    for r in m["outlets"]:
        if r["Type"].upper().startswith("TABULAR"):
            cv_ref.add(r["Qtable/Qcoeff"])
    for n in sorted(x for x in ts_ref if x and x not in m["timeseries"]):
        W.append(f"time series '{n}' is referenced but not defined")
    for n in sorted(x for x in cv_ref if x and x != "*" and x not in m["curves"]):
        W.append(f"curve '{n}' is referenced but not defined")
    for n in sorted(x for x in pat_ref if x and x not in m["patterns"]):
        W.append(f"pattern '{n}' is referenced but not defined")
    # conduit geometry
    adverse = [c["name"] for c in cond if c.get("slope_pct") is not None and c["slope_pct"] < 0]
    flat = [c["name"] for c in cond if c.get("slope_pct") is not None and c["slope_pct"] == 0]
    if adverse:
        W.append(f"{len(adverse)} conduit(s) with adverse slope: {', '.join(adverse[:8])}{' ...' if len(adverse) > 8 else ''}")
    if flat and m["options"].get("FLOW_ROUTING", "").upper() != "DYNWAVE":
        W.append(f"{len(flat)} flat conduit(s) with {m['options'].get('FLOW_ROUTING', 'KINWAVE')} routing")
    for c in cond:
        if c["geom1"] and c.get("z_up") is not None:
            n1 = nodes[c["from"]]
            if n1["type"] == "JUNCTION" and n1["maxdepth"] and c["z_up"] + c["geom1"] > n1["invert"] + n1["maxdepth"] + 1e-6:
                W.append(f"conduit {c['name']}: crown above the rim of upstream node {c['from']}")
    # Dynamic wave, fixed step: conduits shorter than the gravity-wave travel distance in one step
    if m["options"].get("FLOW_ROUTING", "").upper() == "DYNWAVE" and num(m["options"].get("VARIABLE_STEP"), 0) == 0:
        rs_txt = m["options"].get("ROUTING_STEP", "")
        rs = (parse_hours(rs_txt) * 3600) if ":" in rs_txt else num(rs_txt)
        f = 0.3048 if m["units"]["system"] == "US" else 1.0
        if rs:
            short = [c["name"] for c in cond if c["length"] and c["geom1"]
                     and c["length"] * f < math.sqrt(9.81 * c["geom1"] * f) * rs]
            if short:
                W.append(f"{len(short)} conduit(s) shorter than the gravity-wave travel distance in one fixed "
                         f"{g(rs)} s routing step (e.g. {', '.join(short[:5])}) - consider VARIABLE_STEP or LENGTHENING_STEP")
    return W


def summarize(path: Path, inp: dict) -> dict:
    m = build_model(inp)
    nodes, links = node_table(m), link_table(m)
    cond = conduit_geometry(m, nodes)
    U = m["units"]
    opt = m["options"]
    area = sum(num(r["Area"], 0) for r in m["subcatchments"])
    imp = sum(num(r["Area"], 0) * num(r["%Imperv"], 0) for r in m["subcatchments"]) / area if area else None
    sd = parse_date(opt.get("START_DATE", ""))
    ed = parse_date(opt.get("END_DATE", ""))
    period = None
    if sd and ed:
        s = sd + _dt.timedelta(hours=parse_hours(opt.get("START_TIME", "0:00")) or 0)
        e = ed + _dt.timedelta(hours=parse_hours(opt.get("END_TIME", "0:00")) or 0)
        period = f"{s:%Y-%m-%d %H:%M} -> {e:%Y-%m-%d %H:%M} ({(e - s).total_seconds() / 3600:g} h)"
    process = []
    for flag, lab, has in (("IGNORE_RAINFALL", "rainfall/runoff", bool(m["subcatchments"])),
                           ("IGNORE_RDII", "RDII", bool(m["rdii"])),
                           ("IGNORE_SNOWMELT", "snowmelt", bool(table(inp, "SNOWPACKS"))),
                           ("IGNORE_GROUNDWATER", "groundwater", bool(table(inp, "GROUNDWATER"))),
                           ("IGNORE_ROUTING", "flow routing", bool(links)),
                           ("IGNORE_QUALITY", "water quality", bool(m["pollutants"]))):
        if has and opt.get(flag, "NO").upper() != "YES":
            process.append(lab)
    lens = [c["length"] for c in cond]
    sizes = [c["geom1"] for c in cond]
    slopes = [c.get("slope_pct") for c in cond]
    shapes = {}
    for c in cond:
        shapes[c["shape"]] = shapes.get(c["shape"], 0) + 1
    S = {
        "Title": " / ".join(m["title"][:2]) if m["title"] else None,
        "Units": f"{U['system']} - flow {U['flow']}, length {U['length']}, area {U['area']}",
        "Simulation period": period,
        "Report / wet / dry / routing step": " / ".join(opt.get(k, "-") for k in ("REPORT_STEP", "WET_STEP", "DRY_STEP", "ROUTING_STEP")),
        "Processes": ", ".join(process) or "none",
        "Infiltration": opt.get("INFILTRATION", "HORTON (default)") if m["subcatchments"] else None,
        "Flow routing": opt.get("FLOW_ROUTING", "KINWAVE (default)") + (
            f", {opt.get('SURCHARGE_METHOD', 'EXTRAN')} surcharge, variable step {opt.get('VARIABLE_STEP', '0')}"
            if opt.get("FLOW_ROUTING", "").upper() == "DYNWAVE" else ""),
        "Rain gages": len(m["raingages"]) or None,
        "Subcatchments": f"{len(m['subcatchments'])}, {g(area)} {U['area']}, {g(imp)} % impervious (area-weighted)" if m["subcatchments"] else None,
        "Nodes": ", ".join(f"{len(m[k])} {k}" for k in ("junctions", "outfalls", "storage", "dividers") if m[k]) or None,
        "Links": ", ".join(f"{len(m[k])} {k}" for k in ("conduits", "pumps", "orifices", "weirs", "outlets") if m[k]) or None,
        "Conduits": (f"total {g(sum(x for x in lens if x))} {U['length']}; size {rng(sizes, ' ' + U['length'])}; "
                     f"slope {rng([round(s, 3) if s is not None else None for s in slopes], ' %')}; "
                     + ", ".join(f"{n} {s.lower()}" for s, n in shapes.items())) if cond else None,
        "LID controls / usage": f"{len(m['lid_controls'])} / {len(m['lid_usage'])}" if m["lid_controls"] or m["lid_usage"] else None,
        "Pollutants": ", ".join(p["Name"] for p in m["pollutants"]) or None,
        "Inflows (external / DWF / RDII)": f"{len(m['inflows'])} / {len(m['dwf'])} / {len(m['rdii'])}"
        if m["inflows"] or m["dwf"] or m["rdii"] else None,
        "Control rules": len(m["rules"]) or None,
        "Curves / time series / patterns": f"{len(m['curves'])} / {len(m['timeseries'])} / {len(m['patterns'])}",
    }
    rain = rain_summary(m)
    for r in rain:
        if r.get("total_depth"):
            S.setdefault("Design rainfall", f"{r['name']}: {r['total_depth']} in {r['duration']}, peak {r.get('peak')}")
    W = checks(m, nodes, links, cond)
    for r in rain:
        if r.get("problem"):
            W.append(f"rain gage {r['name']}: {r['problem']}")
    sections = []
    for sname in inp["order"]:
        n = len(inp["raw"][sname]) if sname in RAW_SECTIONS else len(inp["sections"].get(sname, []))
        sections.append({"section": sname, "rows": n})
    return {"name": path.stem, "path": str(path), "folder": str(path.parent), "model": m, "nodes": nodes,
            "links": links, "conduit_geometry": cond, "rain": rain, "summary": {k: v for k, v in S.items() if v not in (None, "")},
            "sections": sections, "warnings": W, "_inp": inp}

# --------------------------------------------------------------------------- #
# Report file (.rpt)
# --------------------------------------------------------------------------- #
def rpt_blocks(txt: str) -> dict:
    """Split the report into its '*****'-framed sections."""
    lines = txt.splitlines()
    blocks, i, cur = {}, 0, "header"
    blocks[cur] = []
    while i < len(lines):
        ln = lines[i]
        if re.match(r"\s*\*{5,}", ln) and i + 2 < len(lines) and re.match(r"\s*\*{5,}", lines[i + 2]):
            title = re.sub(r"\s{2,}.*$", "", lines[i + 1].strip())
            cur = title
            blocks[cur] = [lines[i + 1]]
            i += 3
            continue
        blocks.setdefault(cur, []).append(ln)
        i += 1
    return blocks


def _table_rows(block, min_cols=2):
    rows, seen_dash = [], 0
    for ln in block:
        if re.match(r"\s*-{10,}", ln):
            seen_dash += 1
            continue
        if seen_dash >= 2 and ln.strip():
            t = ln.split()
            if len(t) >= min_cols:
                rows.append(t)
        elif seen_dash >= 2 and not ln.strip() and rows:
            break
    return rows


def parse_rpt(path: Path) -> dict:
    txt = path.read_text(errors="replace")
    out = {"file": path.name}
    m = re.search(r"VERSION\s+([\d.]+\s*\(Build[^)]*\))", txt)
    out["engine"] = m.group(1) if m else None
    errs, warns, lines = [], [], txt.splitlines()
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("ERROR"):
            nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
            errs.append(s + (f"  >> {nxt}" if nxt and not nxt.startswith(("ERROR", "WARNING")) else ""))
        elif s.startswith("WARNING"):
            warns.append(s)
    out["errors"], out["warnings"] = errs, list(dict.fromkeys(warns))
    for key, pat in (("analysis begun", r"Analysis begun on:\s*(.+)"), ("analysis ended", r"Analysis ended on:\s*(.+)"),
                     ("elapsed time", r"Total elapsed time:\s*(.+)")):
        mm = re.search(pat, txt)
        if mm:
            out[key] = mm.group(1).strip()
    B = rpt_blocks(txt)
    out["sections"] = [k for k in B if k != "header"]
    cont = {}
    for name in ("Runoff Quantity Continuity", "Flow Routing Continuity", "Runoff Quality Continuity",
                 "Quality Routing Continuity", "Groundwater Continuity"):
        blk = next((v for k, v in B.items() if k.startswith(name)), None)
        if blk:
            vals = []
            terms = {}
            for ln in blk:
                mm = re.match(r"\s*(.+?)\s*\.{2,}\s+(.+)$", ln)
                if mm:
                    nums_ = mm.group(2).split()
                    terms[mm.group(1).strip()] = nums_
                    if mm.group(1).strip().startswith("Continuity Error"):
                        vals = [num(x) for x in nums_]
            cont[name] = {"continuity_error_pct": vals, "terms": terms}
    out["continuity"] = cont
    rts = B.get("Routing Time Step Summary")
    if rts:
        d = {}
        for ln in rts:
            mm = re.match(r"\s*(Minimum Time Step|Average Time Step|Maximum Time Step|% of Time in Steady State|"
                          r"Average Iterations per Step|% of Steps Not Converging)\s*:\s*(.+)", ln)
            if mm:
                d[mm.group(1)] = mm.group(2).strip()
        out["routing_time_step"] = d
    for name, key in (("Highest Continuity Errors", "highest_continuity_errors"),
                      ("Time-Step Critical Elements", "time_step_critical"),
                      ("Highest Flow Instability Indexes", "flow_instability"),
                      ("Most Frequent Nonconverging Nodes", "nonconverging_nodes")):
        blk = B.get(name)
        if blk:
            items = [ln.strip() for ln in blk[1:] if ln.strip()]
            out[key] = items[:10]
    for name, key in (("Node Flooding Summary", "flooding"), ("Node Surcharge Summary", "node_surcharge"),
                      ("Conduit Surcharge Summary", "conduit_surcharge"), ("Outfall Loading Summary", "outfalls"),
                      ("Storage Volume Summary", "storage"), ("Pumping Summary", "pumping")):
        blk = B.get(name)
        if blk:
            if any(re.search(r"No nodes were|No conduits were|No ", ln) for ln in blk[1:4]):
                out[key] = []
            else:
                out[key] = _table_rows(blk)
    return out


# --------------------------------------------------------------------------- #
# Binary results file (.out) - layout from src/solver/output.c
# --------------------------------------------------------------------------- #
MAGIC = 516114522
SUB_VARS = ["rainfall", "snow depth", "evap loss", "infil loss", "runoff", "gw outflow", "gw elevation", "soil moisture"]
NODE_VARS = ["depth", "head", "volume", "lateral inflow", "total inflow", "flooding"]
LINK_VARS = ["flow", "depth", "velocity", "volume", "capacity"]
SYS_VARS = ["air temperature", "rainfall", "snow depth", "infil loss", "runoff", "dry weather inflow", "groundwater inflow",
            "RDII inflow", "external inflow", "total lateral inflow", "flooding", "outfall outflow", "storage volume",
            "evap rate", "potential ET"]
NODE_TYPES = ["JUNCTION", "OUTFALL", "STORAGE", "DIVIDER"]
LINK_TYPES = ["CONDUIT", "PUMP", "ORIFICE", "WEIR", "OUTLET"]


class SwmmOut:
    def __init__(self, path: Path):
        self.path = path
        f = self.f = open(path, "rb")
        f.seek(0, 2)
        size = f.tell()
        if size < 64:
            raise ValueError("file too small")
        f.seek(size - 24)
        self.id_pos, self.in_pos, self.out_pos, self.nperiods, self.error_code, magic2 = struct.unpack("<6i", f.read(24))
        f.seek(0)
        magic, self.version, self.flow_units, self.nsub, self.nnode, self.nlink, self.npoll = struct.unpack("<7i", f.read(28))
        if magic != MAGIC:
            raise ValueError("not a SWMM binary output file (bad magic number)")
        self.complete = magic2 == MAGIC
        names = [self._id() for _ in range(self.nsub + self.nnode + self.nlink + self.npoll)]
        self.sub_names = names[:self.nsub]
        self.node_names = names[self.nsub:self.nsub + self.nnode]
        self.link_names = names[self.nsub + self.nnode:self.nsub + self.nnode + self.nlink]
        self.poll_names = names[self.nsub + self.nnode + self.nlink:]
        self.poll_units = list(struct.unpack(f"<{self.npoll}i", f.read(4 * self.npoll))) if self.npoll else []
        # input data
        n = self._i()
        f.read(4 * n)
        self.sub_area = list(struct.unpack(f"<{self.nsub}f", f.read(4 * self.nsub))) if self.nsub else []
        n = self._i()
        f.read(4 * n)
        self.node_info = []
        for _ in range(self.nnode):
            t = self._i()
            inv, md = struct.unpack("<2f", f.read(8))
            self.node_info.append((t, inv, md))
        n = self._i()
        f.read(4 * n)
        self.link_info = []
        for _ in range(self.nlink):
            t = self._i()
            self.link_info.append((t,) + struct.unpack("<4f", f.read(16)))
        self.nsv = self._i()
        f.read(4 * self.nsv)
        self.nnv = self._i()
        f.read(4 * self.nnv)
        self.nlv = self._i()
        f.read(4 * self.nlv)
        self.nsys = self._i()
        f.read(4 * self.nsys)
        self.start = struct.unpack("<d", f.read(8))[0]
        self.report_step = self._i()
        self.data_pos = f.tell()
        self.nvals = self.nsub * self.nsv + self.nnode * self.nnv + self.nlink * self.nlv + self.nsys
        self.bytes_per_period = 8 + 4 * self.nvals
        if not self.complete:      # run aborted: count the periods actually written
            self.nperiods = max(0, (size - self.data_pos) // self.bytes_per_period)

    def close(self):
        self.f.close()

    def _i(self):
        return struct.unpack("<i", self.f.read(4))[0]

    def _id(self):
        n = self._i()
        return self.f.read(n).decode("utf-8", "replace")

    def periods(self):
        f = self.f
        f.seek(self.data_pos)
        for _ in range(self.nperiods):
            buf = f.read(self.bytes_per_period)
            if len(buf) < self.bytes_per_period:
                return
            t = struct.unpack("<d", buf[:8])[0]
            a = array.array("f")
            a.frombytes(buf[8:])
            if sys.byteorder != "little":
                a.byteswap()
            yield t, a

    def offsets(self):
        s = 0
        n = s + self.nsub * self.nsv
        l = n + self.nnode * self.nnv
        y = l + self.nlink * self.nlv
        return s, n, l, y


def read_out(path: Path, max_series=4000) -> dict:
    o = SwmmOut(path)
    try:
        res = {"file": path.name, "size_MB": round(path.stat().st_size / 1e6, 2), "version": o.version,
               "flow_units": FLOW_UNITS[o.flow_units] if 0 <= o.flow_units < 6 else o.flow_units,
               "counts": {"subcatchments": o.nsub, "nodes": o.nnode, "links": o.nlink, "pollutants": o.npoll},
               "pollutants": o.poll_names, "report_start": swmm_datetime(o.start).strftime("%Y-%m-%d %H:%M:%S"),
               "report_step_s": o.report_step, "periods": o.nperiods, "error_code": o.error_code,
               "complete": "yes" if o.complete else "NO - file not closed properly (run aborted?)",
               "system_variables": o.nsys}
        s0, n0, l0, y0 = o.offsets()
        sub_max = [[-1e30, 0] for _ in range(o.nsub)]
        node_max = [[-1e30, 0, -1e30, -1e30, 0.0] for _ in range(o.nnode)]   # depth, t, head, inflow, flood volume
        node_tflood = [0] * o.nnode
        link_max = [[0.0, 0, 0.0, 0.0] for _ in range(o.nlink)]               # |flow|, t, velocity, capacity
        sys_series, times = [], []
        stride = max(1, o.nperiods // max_series)
        dt = o.report_step
        for k, (t, a) in enumerate(o.periods()):
            for i in range(o.nsub):
                v = a[s0 + i * o.nsv + 4]
                if v > sub_max[i][0]:
                    sub_max[i] = [v, t]
            for i in range(o.nnode):
                b = n0 + i * o.nnv
                d, h, qin, fl = a[b], a[b + 1], a[b + 4], a[b + 5]
                nm = node_max[i]
                if d > nm[0]:
                    nm[0], nm[1] = d, t
                if h > nm[2]:
                    nm[2] = h
                if qin > nm[3]:
                    nm[3] = qin
                if fl > 0:
                    nm[4] += fl * dt
                    node_tflood[i] += dt
            for i in range(o.nlink):
                b = l0 + i * o.nlv
                q, vel, cap = abs(a[b]), abs(a[b + 2]), a[b + 4]
                lm = link_max[i]
                if q > lm[0]:
                    lm[0], lm[1] = q, t
                if vel > lm[2]:
                    lm[2] = vel
                if cap > lm[3]:
                    lm[3] = cap
            if k % stride == 0:
                times.append(t)
                sys_series.append(list(a[y0:y0 + o.nsys]))
        res["_sub_max"] = {o.sub_names[i]: sub_max[i] for i in range(o.nsub)}
        res["_node_max"] = {o.node_names[i]: node_max[i] + [node_tflood[i]] for i in range(o.nnode)}
        res["_link_max"] = {o.link_names[i]: link_max[i] for i in range(o.nlink)}
        res["_node_info"] = {o.node_names[i]: o.node_info[i] for i in range(o.nnode)}
        res["_link_info"] = {o.link_names[i]: o.link_info[i] for i in range(o.nlink)}
        res["_sys"] = {"hours": [(t - o.start) * 24 for t in times], "values": sys_series}
        if times:
            res["first_result"] = swmm_datetime(times[0]).strftime("%Y-%m-%d %H:%M")
            res["last_result"] = swmm_datetime(times[-1]).strftime("%Y-%m-%d %H:%M")
        # headline numbers
        if sys_series:
            col = lambda j: [r[j] for r in sys_series] if o.nsys > j else []
            res["peak system runoff"] = g(max(col(4) or [0]))
            res["peak outfall outflow"] = g(max(col(11) or [0]))
            res["peak flooding rate"] = g(max(col(10) or [0]))
        fl = [(n, v[4]) for n, v in res["_node_max"].items() if v[4] > 0]
        res["flooded nodes"] = len(fl)
        sur = [n for n, v in res["_link_max"].items() if v[3] >= 0.999]
        res["links at full capacity"] = len(sur)
        return res
    finally:
        o.close()


def _safe_read_out(path):
    try:
        return read_out(path)
    except Exception as e:
        return {"file": path.name, "error": f"{type(e).__name__}: {e}"}


# --------------------------------------------------------------------------- #
# Profiles
# --------------------------------------------------------------------------- #
def profile_data(cfg, links_seq) -> dict | None:
    nodes, links, cond = cfg["nodes"], cfg["links"], {c["name"]: c for c in cfg["conduit_geometry"]}
    out = cfg.get("out") or {}
    nmax = out.get("_node_max", {})
    pts, x = [], 0.0
    prev = None
    for ln in links_seq:
        L = links.get(ln)
        if not L:
            return None
        a, b = L["from"], L["to"]
        if prev is not None and a != prev and b == prev:
            a, b = b, a
        c = cond.get(ln, {})
        length = c.get("length") or 0.0
        na, nb = nodes.get(a), nodes.get(b)
        if not na or not nb:
            return None
        if not pts:
            pts.append({"node": a, "x": 0.0, "invert": na["invert"], "rim": na["invert"] + na["maxdepth"] if na["invert"] is not None else None,
                        "hgl": nmax.get(a, [None, None, None])[2]})
        x += length
        pts.append({"node": b, "x": x, "invert": nb["invert"], "rim": nb["invert"] + nb["maxdepth"] if nb["invert"] is not None else None,
                    "hgl": nmax.get(b, [None, None, None])[2]})
        pts[-1]["link"] = {"name": ln, "z1": c.get("z_up"), "z2": c.get("z_dn"), "h": c.get("geom1"), "x1": x - length, "x2": x}
        prev = b
    return {"points": pts}


# --------------------------------------------------------------------------- #
# Model discovery & build
# --------------------------------------------------------------------------- #
def find_inputs(arg: str) -> list[Path]:
    p = Path(arg).expanduser().resolve()
    if not p.exists():
        sys.exit(f"Path not found: {p}")
    if p.is_file():
        return [p]
    fs = sorted(p.rglob("*.inp"))
    if not fs:
        sys.exit(f"No .inp file found under {p}")
    return fs


def build_one(path: Path, rpt=None, out=None) -> dict:
    if path.suffix.lower() == ".out":
        cfg = {"name": path.stem, "path": str(path), "folder": str(path.parent), "summary": {}, "warnings": [],
               "sections": [], "only_out": True}
        cfg["out"] = _safe_read_out(path)
        return cfg
    inp = read_inp(path)
    cfg = summarize(path, inp)
    rp = Path(rpt) if rpt else path.with_suffix(".rpt")
    op = Path(out) if out else path.with_suffix(".out")
    if rp.exists():
        try:
            cfg["rpt"] = parse_rpt(rp)
        except OSError as e:
            cfg["rpt"] = {"file": rp.name, "error": str(e)}
    if op.exists():
        cfg["out"] = _safe_read_out(op)
    W = cfg["warnings"]
    r = cfg.get("rpt")
    if r and not r.get("error"):
        if r["errors"]:
            W.append(f"{r['file']}: {len(r['errors'])} ERROR(s) - {r['errors'][0][:120]}")
        for k, v in r.get("continuity", {}).items():
            errs = [x for x in v["continuity_error_pct"] if x is not None]
            if errs and abs(errs[0]) > 1:
                W.append(f"{k} error {g(errs[0])} % (> 1 %)")
        if r.get("flooding"):
            W.append(f"{len(r['flooding'])} node(s) flooded (see Node Flooding Summary)")
        nc = (r.get("routing_time_step") or {}).get("% of Steps Not Converging")
        if nc and num(nc.split()[0], 0) > 1:
            W.append(f"{nc} of routing steps did not converge")
    o = cfg.get("out")
    if o and o.get("error"):
        W.append(f"{o['file']}: {o['error']}")
    elif o:
        if o.get("error_code"):
            W.append(f"{o['file']}: run ended with error code {o['error_code']}")
        if not str(o.get("complete", "")).startswith("yes"):
            W.append(f"{o['file']}: {o['complete']}")
        inp_n = (len(cfg["model"]["subcatchments"]), len(cfg["nodes"]), len(cfg["links"]))
        if (o["counts"]["subcatchments"], o["counts"]["nodes"], o["counts"]["links"]) != inp_n:
            W.append(f"{o['file']}: reports {o['counts']['subcatchments']}/{o['counts']['nodes']}/{o['counts']['links']} "
                     f"subcatchments/nodes/links, the .inp has {inp_n[0]}/{inp_n[1]}/{inp_n[2]} - output from another version of the model, "
                     f"or [REPORT] limits the elements saved")
        try:
            if path.stat().st_mtime > op.stat().st_mtime + 1:
                W.append(f"{path.name} is newer than {op.name} - results may be out of date")
        except OSError:
            pass
    prof = {}
    for name, seq in cfg["model"]["profiles"].items():
        d = profile_data(cfg, seq)
        if d:
            prof[name] = d
        else:
            W.append(f"profile '{name}' references unknown links or nodes")
    cfg["profiles"] = prof
    return cfg


def build_config(arg, rpt=None, out=None) -> dict:
    files = find_inputs(arg)
    root = Path(arg).expanduser().resolve()
    return {"root": str(root if root.is_dir() else root.parent), "name": root.stem if root.is_file() else root.name,
            "models": [build_one(f, rpt, out) for f in files]}


# --------------------------------------------------------------------------- #
# HTML template & helpers (same layout as ras_config.py / hms_config.py)
# --------------------------------------------------------------------------- #
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<script>(function(){var t;try{t=localStorage.getItem('modelvue-theme')}catch(e){}if(t!=='dark')t='light';document.documentElement.dataset.theme=t;})();</script>
<title>SWMM · @@TITLE@@</title>
<style>
:root{--bg:#f6f7f9;--panel:#fff;--ink:#1d2433;--muted:#667085;--line:#e4e7ec;--acc:#1e5aa8;--accs:#e6eef9;--warn:#b42318}
:root[data-theme=dark]{color-scheme:dark;--bg:#0f1318;--panel:#171c23;--ink:#e6e9ef;--muted:#98a2b3;--line:#2a313b;--acc:#8fb4ec;--accs:#18243a;--warn:#f97066}
*{box-sizing:border-box} body{margin:0;font:14px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,"Noto Sans TC",sans-serif;background:var(--bg);color:var(--ink)}
header{position:sticky;top:0;z-index:5;background:var(--acc);color:#fff;padding:14px 24px;display:flex;gap:16px;align-items:baseline;flex-wrap:wrap}
header h1{margin:0;font-size:20px} header span{opacity:.85}
[data-theme=dark] header{background:#173a6b}
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
<header><h1>SWMM · @@TITLE@@</h1><span>@@SUB@@</span><button class="theme" type="button" title="Switch between the light and dark theme">☾ Dark</button></header>
<div class="wrap"><nav><div class="tools"><button data-x="1">Expand all</button><button data-x="0">Collapse all</button></div>@@NAV@@</nav>
<main>@@BODY@@
<footer>Generated by modelvue.swmmvue on @@GEN@@ from <code>@@PATH@@</code></footer></main></div>
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


def row(label, value, indent=2):
    if value in (None, "", [], {}):
        return
    if isinstance(value, dict):
        value = ", ".join(f"{k}: {v}" for k, v in value.items())
    elif isinstance(value, list):
        value = ", ".join(map(str, value))
    print(f"{' ' * indent}{label:<32}: {value}")


def text_table(headers, rows, indent=2, maxw=24, limit=None):
    if not rows:
        return
    shown = rows if limit is None else rows[:limit]
    w = [min(maxw, max(len(str(h)), *(len(str(r[i] if i < len(r) else "")) for r in shown))) for i, h in enumerate(headers)]
    pad = " " * indent
    print(pad + "  ".join(str(h)[:w[i]].ljust(w[i]) for i, h in enumerate(headers)))
    print(pad + "  ".join("-" * x for x in w))
    for r in shown:
        print(pad + "  ".join(str(r[i] if i < len(r) else "")[:w[i]].ljust(w[i]) for i in range(len(headers))))
    if limit is not None and len(rows) > limit:
        print(f"{pad}... {len(rows) - limit} more (use --all)")


def section_rows(cfg, sec):
    inp = cfg["_inp"]
    m = cfg["model"]
    if sec == "OUTFALLS":
        return COLS["OUTFALLS"], [[r[c] for c in COLS["OUTFALLS"]] for r in m["outfalls"]]
    if sec == "STORAGE":
        return COLS["STORAGE"], [[r[c] for c in COLS["STORAGE"]] for r in m["storage"]]
    if sec == "INFILTRATION":
        rows = list(m["infiltration"].values())
        if not rows:
            return [], []
        cols = [c for c in rows[0] if c != "method"] + (["method"] if len({r["method"] for r in rows}) > 1 else [])
        return cols, [[r.get(c, "") for c in cols] for r in rows]
    rows = table(inp, sec)
    cols = COLS.get(sec) or [f"c{i + 1}" for i in range(max((len(r) for r in rows), default=0))]
    n = max((len(r) for r in rows), default=0)
    if n > len(cols):
        cols = cols + [f"c{i + 1}" for i in range(len(cols), n)]
    return cols, [r + [""] * (len(cols) - len(r)) for r in rows]


def show_model(cfg, full):
    h1(f"SWMM MODEL: {cfg['name']}")
    row("File", cfg["path"])
    if cfg.get("only_out"):
        show_out(cfg)
        return
    m = cfg["model"]
    for t in m["title"]:
        print(f"  | {t}")
    h2("Summary")
    for k, v in cfg["summary"].items():
        row(k, v)
    h2("Options")
    for k, v in m["options"].items():
        row(k, f"{v}" + (f"   ({OPTION_HELP[k]})" if OPTION_HELP.get(k) else ""))
    h2("Sections in the file")
    print("  " + ", ".join(f"{s['section']} ({s['rows']})" for s in cfg["sections"]))
    if cfg["rain"]:
        h2("Rainfall")
        for r in cfg["rain"]:
            print(f"  {r['name']}: {r['format']} every {r['interval']}, SCF {r['scf']}, {r['source']} {r['data']}")
            for k in ("total_depth", "peak", "duration", "starts", "note", "problem"):
                row(k, r.get(k), indent=6)
    lim = None if full else 12
    for sec in ("SUBCATCHMENTS", "SUBAREAS", "INFILTRATION", "JUNCTIONS", "OUTFALLS", "STORAGE", "DIVIDERS", "CONDUITS",
                "PUMPS", "ORIFICES", "WEIRS", "OUTLETS", "XSECTIONS", "LOSSES", "LID_USAGE", "INFLOWS", "DWF", "RDII",
                "POLLUTANTS", "LANDUSES"):
        cols, rows = section_rows(cfg, sec)
        if rows:
            h2(f"[{sec}] ({len(rows)})")
            text_table(cols, rows, limit=lim, maxw=16)
    cg = cfg["conduit_geometry"]
    if cg:
        h2("Conduit geometry (computed)")
        text_table(["Conduit", "From", "To", "Length", "Shape", "Geom1", "Inlet z", "Outlet z", "Slope %"],
                   [[c["name"], c["from"], c["to"], g(c["length"]), c["shape"], g(c["geom1"]), g(c.get("z_up")),
                     g(c.get("z_dn")), g(c.get("slope_pct"))] for c in cg], limit=lim, maxw=14)
    if m["lid_controls"]:
        h2("LID controls")
        for n, d in m["lid_controls"].items():
            print(f"  {n}: {d['type']} ({LID_TYPES.get(d['type'], '?')}) - layers {', '.join(d['layers'])}")
    if m["curves"]:
        h2("Curves")
        for n, c in m["curves"].items():
            print(f"  {n:<20} {c['type']:<10} {len(c['points'])} points")
    if m["timeseries"]:
        h2("Time series")
        for n, t in m["timeseries"].items():
            p = t["points"]
            desc = f"file {t['file']}" if t["file"] else (f"{len(p)} points, {g(p[0][0])} .. {g(p[-1][0])} h, value {rng([v for _, v in p])}" if p else "empty")
            print(f"  {n:<20} {desc}" + (f", starts {t['start']}" if t.get("start") else ""))
    if m["patterns"]:
        h2("Patterns")
        for n, p in m["patterns"].items():
            print(f"  {n:<20} {p['type']:<8} {len(p['factors'])} factors, {rng(p['factors'])}")
    if m["rules"]:
        h2(f"Control rules ({len(m['rules'])})")
        if full:
            for ln in m["controls"]:
                print("  " + ln)
        else:
            print("  " + ", ".join(m["rules"]))
    if m["profiles"]:
        h2("Profiles")
        for n, seq in m["profiles"].items():
            print(f"  {n}: {' -> '.join(seq)}")
    if m["report"]:
        h2("[REPORT]")
        row("settings", m["report"])
    show_rpt(cfg)
    show_out(cfg)
    if cfg["warnings"]:
        h2("Warnings")
        for w in cfg["warnings"]:
            print(f"  ! {w}")


def show_rpt(cfg):
    r = cfg.get("rpt")
    h2("Report file (.rpt)")
    if not r:
        print("  none found (model not run, or report saved elsewhere - use --rpt)")
        return
    row("file", r["file"])
    for k in ("engine", "analysis begun", "analysis ended", "elapsed time"):
        row(k, r.get(k))
    for x in r["errors"]:
        print(f"    {x}")
    for x in r["warnings"][:15]:
        print(f"    {x}")
    for k, v in r.get("continuity", {}).items():
        row(f"{k} error %", ", ".join(g(x) for x in v["continuity_error_pct"]))
    row("routing time step", r.get("routing_time_step"))
    row("flooded nodes", len(r["flooding"]) if "flooding" in r else None)
    if r.get("flooding"):
        text_table(["Node", "Hours", "Max rate", "Day", "hh:mm", "Total vol", "Max ponded"], r["flooding"][:12], indent=4)
    row("surcharged conduits", len(r.get("conduit_surcharge", [])) if "conduit_surcharge" in r else None)
    for k in ("highest_continuity_errors", "time_step_critical", "flow_instability", "nonconverging_nodes"):
        if r.get(k):
            row(k.replace("_", " "), "; ".join(r[k][:5]))


def show_out(cfg):
    o = cfg.get("out")
    h2("Binary results (.out)")
    if not o:
        print("  none found")
        return
    if o.get("error"):
        print(f"  {o['file']}: {o['error']}")
        return
    for k, v in o.items():
        if not k.startswith("_"):
            row(k, v)
    nm = o["_node_max"]
    top = sorted(nm.items(), key=lambda kv: -kv[1][4])
    fl = [(n, v) for n, v in top if v[4] > 0][:10]
    if fl:
        print("    Most flooded nodes (volume from .out, flow units x s)")
        for n, v in fl:
            print(f"      {n:<20} {g(v[4]):>12}   flooded {g(v[5] / 3600)} h")
    lk = sorted(o["_link_max"].items(), key=lambda kv: -kv[1][3])[:10]
    if lk:
        print("    Links by max capacity used (depth/full depth)")
        for n, v in lk:
            print(f"      {n:<20} capacity {g(v[3]):>6}   max |flow| {g(v[0]):>10}   max velocity {g(v[2])}")


def show(cfg, full=False):
    for m in cfg["models"]:
        show_model(m, full)
    print()


# --------------------------------------------------------------------------- #
# SVG graphics
# --------------------------------------------------------------------------- #
LINK_COLORS = {"CONDUIT": "#1f5fa8", "PUMP": "#c2185b", "ORIFICE": "#00897b", "WEIR": "#ef6c00", "OUTLET": "#6d4c41"}
NODE_STYLE = {"JUNCTION": ("circle", "#1f5fa8"), "OUTFALL": ("tri", "#2e7d32"), "STORAGE": ("rect", "#6a1b9a"),
              "DIVIDER": ("diamond", "#ef6c00")}


def _cap_color(c):
    c = max(0.0, min(1.0, c))
    stops = [(0, (44, 123, 182)), (0.5, (255, 214, 0)), (0.8, (244, 109, 67)), (1, (165, 0, 38))]
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if c <= t1:
            f = (c - t0) / (t1 - t0)
            return "#%02x%02x%02x" % tuple(round(a + (b - a) * f) for a, b in zip(c0, c1))
    return "#a50026"


def svg_map(cfg, results=False, w=900) -> str:
    m, nodes, links = cfg["model"], cfg["nodes"], cfg["links"]
    C = m["coords"]
    pts = list(C.values()) + [p for v in m["vertices"].values() for p in v] + [p for v in m["polygons"].values() for p in v] \
        + list(m["symbols"].values())
    pts = [p for p in pts if p[0] is not None and p[1] is not None]
    if not pts:
        return '<p class="muted">No [COORDINATES] in the input file - no map.</p>'
    x0, x1 = min(p[0] for p in pts), max(p[0] for p in pts)
    y0, y1 = min(p[1] for p in pts), max(p[1] for p in pts)
    span = max(x1 - x0, y1 - y0) or 1
    pad = span * 0.06
    x0, x1, y0, y1 = x0 - pad, x1 + pad, y0 - pad, y1 + pad
    sc = w / (x1 - x0)
    h = (y1 - y0) * sc
    if h > 760:
        sc = 760 / (y1 - y0)
        h, w = 760, (x1 - x0) * sc
    X = lambda v: (v - x0) * sc
    Y = lambda v: (y1 - v) * sc
    o = [f'<svg class="map" viewBox="0 0 {w:.0f} {h:.0f}">']
    out = cfg.get("out") if results else None
    lmax = (out or {}).get("_link_max", {})
    nmax = (out or {}).get("_node_max", {})
    # subcatchments
    sub_out = CI({r["Name"]: r["Outlet"] for r in m["subcatchments"]})
    for name, poly in m["polygons"].items():
        P = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in poly)
        o.append(f'<polygon points="{P}" class="sub"><title>subcatchment {_h.escape(name)}</title></polygon>')
        cx, cy = sum(a for a, _ in poly) / len(poly), sum(b for _, b in poly) / len(poly)
        tgt = sub_out.get(name)
        tp = C.get(tgt) or (lambda pp: (sum(a for a, _ in pp) / len(pp), sum(b for _, b in pp) / len(pp)))(m["polygons"][tgt]) \
            if tgt in C or tgt in m["polygons"] else None
        if tp:
            o.append(f'<line x1="{X(cx):.1f}" y1="{Y(cy):.1f}" x2="{X(tp[0]):.1f}" y2="{Y(tp[1]):.1f}" class="subl"/>')
        o.append(f'<text x="{X(cx):.1f}" y="{Y(cy):.1f}" class="slab">{_h.escape(name)}</text>')
    # links
    for name, L in links.items():
        a, b = C.get(L["from"]), C.get(L["to"])
        if not a or not b:
            continue
        P = [a] + m["vertices"].get(name, []) + [b]
        col = LINK_COLORS.get(L["type"], "#555")
        tip = f"{L['type'].lower()} {name}: {L['from']} -> {L['to']}"
        width = 2.4
        if results and name in lmax:
            cap = lmax[name][3]
            col = _cap_color(cap) if L["type"] == "CONDUIT" else col
            width = 2 + 3 * min(1, cap)
            tip += f"; max capacity {g(cap)}, max |flow| {g(lmax[name][0])}"
        dash = ' stroke-dasharray="6 3"' if L["type"] in ("PUMP", "ORIFICE", "WEIR", "OUTLET") else ""
        o.append(f'<polyline points="{" ".join(f"{X(p[0]):.1f},{Y(p[1]):.1f}" for p in P)}" fill="none" stroke="{col}" '
                 f'stroke-width="{width:.1f}"{dash} stroke-linecap="round"><title>{_h.escape(tip)}</title></polyline>')
    # nodes
    many = len(nodes) > 80
    for name, n in nodes.items():
        p = C.get(name)
        if not p:
            continue
        shape, col = NODE_STYLE.get(n["type"], ("circle", "#555"))
        x, y = X(p[0]), Y(p[1])
        r = 3 if many else 5
        tip = f"{n['type'].lower()} {name}: invert {g(n['invert'])}"
        extra = ""
        if results and name in nmax:
            v = nmax[name]
            tip += f"; max depth {g(v[0])}, max head {g(v[2])}" + (f"; FLOODED {g(v[5] / 3600)} h" if v[4] > 0 else "")
            if v[4] > 0:
                extra = f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r + 5}" class="flood"/>'
        o.append(extra)
        t = f"<title>{_h.escape(tip)}</title>"
        if shape == "circle":
            o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{col}" class="nd">{t}</circle>')
        elif shape == "tri":
            o.append(f'<polygon points="{x:.1f},{y - r - 2:.1f} {x - r - 2:.1f},{y + r:.1f} {x + r + 2:.1f},{y + r:.1f}" fill="{col}" class="nd">{t}</polygon>')
        elif shape == "rect":
            o.append(f'<rect x="{x - r - 1:.1f}" y="{y - r - 1:.1f}" width="{2 * r + 2}" height="{2 * r + 2}" fill="{col}" class="nd">{t}</rect>')
        else:
            o.append(f'<polygon points="{x:.1f},{y - r - 2:.1f} {x + r + 2:.1f},{y:.1f} {x:.1f},{y + r + 2:.1f} {x - r - 2:.1f},{y:.1f}" fill="{col}" class="nd">{t}</polygon>')
        if not many:
            o.append(f'<text x="{x + 7:.1f}" y="{y - 6:.1f}" class="nlab">{_h.escape(name)}</text>')
    for name, p in m["symbols"].items():
        x, y = X(p[0]), Y(p[1])
        o.append(f'<g class="gage"><rect x="{x - 7:.1f}" y="{y - 9:.1f}" width="14" height="10" rx="2"/>'
                 f'<line x1="{x - 4:.1f}" y1="{y + 3:.1f}" x2="{x - 6:.1f}" y2="{y + 7:.1f}"/>'
                 f'<line x1="{x:.1f}" y1="{y + 3:.1f}" x2="{x - 2:.1f}" y2="{y + 7:.1f}"/>'
                 f'<line x1="{x + 4:.1f}" y1="{y + 3:.1f}" x2="{x + 2:.1f}" y2="{y + 7:.1f}"/>'
                 f'<title>rain gage {_h.escape(name)}</title></g>')
    o.append("</svg>")
    leg = ['<div class="legend">']
    glyph = {"circle": "\u25cf", "tri": "\u25b2", "rect": "\u25a0", "diamond": "\u25c6"}
    for t, (shape, col) in NODE_STYLE.items():
        if any(n["type"] == t for n in nodes.values()):
            leg.append(f'<span><b style="color:{col};margin-right:5px">{glyph[shape]}</b>{t.lower()}</span>')
    if results and lmax:
        grad = ",".join(_cap_color(t) for t in (0, .25, .5, .75, 1))
        leg.append(f'<span><i class="sw" style="background:linear-gradient(90deg,{grad})"></i>conduit max capacity 0 .. 1 (width too)</span>')
        leg.append('<span><i class="mk flood"></i>flooded node</span>')
    else:
        for t, col in LINK_COLORS.items():
            if any(L["type"] == t for L in links.values()):
                leg.append(f'<span><i class="ln" style="border-color:{col}"></i>{t.lower()}</span>')
    if m["polygons"]:
        leg.append('<span><i class="mk" style="background:#cfe3c8;border:1px solid #7a9f6f"></i>subcatchment (dashed: to outlet)</span>')
    leg.append("</div>")
    return "".join(o) + "".join(leg)


def _frame(w, h, L, R, T, B, x0, x1, y0, y1, xlabel, ylabel, flip=False):
    X = lambda v: L + (v - x0) / ((x1 - x0) or 1) * (w - L - R)
    Y = (lambda v: T + (v - y0) / ((y1 - y0) or 1) * (h - T - B)) if flip else (lambda v: T + (y1 - v) / ((y1 - y0) or 1) * (h - T - B))
    o = []
    for t in _ticks(y0, y1):
        o.append(f'<line x1="{L}" x2="{w - R}" y1="{Y(t):.1f}" y2="{Y(t):.1f}" class="gl"/>'
                 f'<text x="{L - 6}" y="{Y(t) + 4:.1f}" text-anchor="end">{g(t)}</text>')
    for t in _ticks(x0, x1, 8):
        o.append(f'<text x="{X(t):.1f}" y="{h - B + 15}" text-anchor="middle">{g(t)}</text>')
    o.append(f'<rect x="{L}" y="{T}" width="{w - L - R}" height="{h - T - B}" class="frame"/>')
    o.append(f'<text x="{(L + w - R) / 2}" y="{h - 3}" text-anchor="middle" class="lab">{_h.escape(xlabel)}</text>')
    o.append(f'<text transform="translate(13,{(T + h - B) / 2}) rotate(-90)" text-anchor="middle" class="lab">{_h.escape(ylabel)}</text>')
    return X, Y, o


def svg_hyeto(r, U, w=820, h=200) -> str:
    s = r.get("_series")
    if not s:
        return ""
    p, dt, fmt = s
    vals = [v for _, v in p]
    if fmt == "VOLUME":
        vals = [v / dt for v in vals]
    elif fmt == "CUMULATIVE":
        vals = [max(0, (p[i][1] - (p[i - 1][1] if i else 0))) / dt for i in range(len(p))]
    x0, x1 = p[0][0], p[-1][0] + dt
    X, Y, o = _frame(w, h, 52, 12, 10, 34, x0, x1, 0, max(vals) * 1.1 or 1, "hours", f"intensity ({U['rate']})", flip=True)
    bars = []
    for (t, _), v in zip(p, vals):
        bars.append(f'<rect x="{X(t):.1f}" y="{Y(0):.1f}" width="{max(1, X(t + dt) - X(t) - 0.5):.1f}" height="{Y(v) - Y(0):.1f}" '
                    f'class="bar"><title>{g(t)} h: {g(v)}</title></rect>')
    return f'<svg class="plot" viewBox="0 0 {w} {h}">' + "".join(bars + o) + "</svg>"


def svg_system(out, U, w=820, h=240) -> str:
    S = out.get("_sys")
    if not S or not S["hours"]:
        return ""
    hrs, V = S["hours"], S["values"]
    want = [(4, "runoff", "#2e7d32"), (9, "total lateral inflow", "#1f5fa8"), (11, "outfall outflow", "#6a1b9a"),
            (10, "flooding", "#c62828")]
    series = [(lab, col, [r[j] for r in V]) for j, lab, col in want if len(V[0]) > j and any(r[j] for r in V)]
    if not series:
        return ""
    ymax = max(max(s[2]) for s in series) * 1.08 or 1
    X, Y, o = _frame(w, h, 58, 12, 14, 34, hrs[0], hrs[-1], 0, ymax, "hours since report start", f"flow ({out['flow_units']})")
    for lab, col, ys in series:
        o.append(f'<polyline points="{" ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in zip(hrs, ys))}" fill="none" stroke="{col}" '
                 f'stroke-width="1.8"><title>{lab}</title></polyline>')
    lx = 68
    for lab, col, _ in series:
        o.append(f'<rect x="{lx}" y="20" width="12" height="3" fill="{col}"/><text x="{lx + 16}" y="25">{lab}</text>')
        lx += 30 + 6.5 * len(lab)
    rain = [r[1] for r in V] if len(V[0]) > 1 else []
    return f'<svg class="plot" viewBox="0 0 {w} {h}">' + "".join(o) + "</svg>" + (
        f'<p class="muted" style="font-size:12px">area-averaged rainfall peak {g(max(rain))} {U["rate"]}</p>' if rain and max(rain) else "")


def svg_profile(name, d, U, w=880, h=300) -> str:
    P = d["points"]
    zs = [v for p in P for v in (p["invert"], p["rim"], p["hgl"]) if v is not None]
    for p in P:
        lk = p.get("link")
        if lk and lk["z1"] is not None and lk["h"]:
            zs += [lk["z1"] + lk["h"], lk["z2"] + lk["h"]]
    if not zs:
        return ""
    z0, z1 = min(zs), max(zs)
    pad = (z1 - z0) * 0.08 or 1
    X, Y, o = _frame(w, h, 60, 14, 14, 34, P[0]["x"], P[-1]["x"] or 1, z0 - pad, z1 + pad,
                     f"distance ({U['length']})", f"elevation ({U['length']})")
    ground = " ".join(f"{X(p['x']):.1f},{Y(p['rim']):.1f}" for p in P if p["rim"] is not None)
    o.append(f'<polyline points="{ground}" fill="none" stroke="#8d6e63" stroke-width="1.5" stroke-dasharray="5 3"><title>node rim</title></polyline>')
    for p in P:
        lk = p.get("link")
        if lk and lk["z1"] is not None:
            a, b = lk["x1"], lk["x2"]
            hh = lk["h"] or 0
            o.append(f'<polygon points="{X(a):.1f},{Y(lk["z1"]):.1f} {X(b):.1f},{Y(lk["z2"]):.1f} {X(b):.1f},{Y(lk["z2"] + hh):.1f} '
                     f'{X(a):.1f},{Y(lk["z1"] + hh):.1f}" class="pipe"><title>{_h.escape(lk["name"])}</title></polygon>')
    for p in P:
        if p["invert"] is not None and p["rim"] is not None:
            o.append(f'<line x1="{X(p["x"]):.1f}" x2="{X(p["x"]):.1f}" y1="{Y(p["invert"]):.1f}" y2="{Y(p["rim"]):.1f}" class="mh"/>'
                     f'<text x="{X(p["x"]):.1f}" y="{Y(p["rim"]) - 5:.1f}" text-anchor="middle" class="nlab">{_h.escape(p["node"])}</text>')
    hg = [(p["x"], p["hgl"]) for p in P if p["hgl"] is not None]
    if len(hg) > 1:
        o.append(f'<polyline points="{" ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in hg)}" fill="none" stroke="#1565c0" '
                 f'stroke-width="2.2"><title>maximum HGL (from .out)</title></polyline>')
    leg = ('<div class="legend"><span><i class="ln" style="border-color:#8d6e63;border-top-style:dashed"></i>rim / ground</span>'
           '<span><i class="mk" style="background:#d7ccc8;border:1px solid #6d4c41"></i>conduit</span>'
           + ('<span><i class="ln" style="border-color:#1565c0"></i>maximum HGL</span>' if len(hg) > 1 else "") + "</div>")
    return f'<h4>Profile "{_h.escape(name)}"</h4><svg class="plot" viewBox="0 0 {w} {h}">' + "".join(o) + "</svg>" + leg


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #
EXTRA_CSS = """
svg.map{width:100%;max-width:1000px;height:auto;display:block;margin:8px 0;background:var(--bg);border:1px solid var(--line);border-radius:6px}
svg.map .sub{fill:#cfe3c8;fill-opacity:.55;stroke:#7a9f6f;stroke-width:1} svg.map .subl{stroke:#7a9f6f;stroke-dasharray:4 3}
svg.map .slab{font-size:11px;fill:#33691e;text-anchor:middle} svg.map .nlab{font-size:10px;fill:var(--ink)}
svg.map .nd{stroke:#fff;stroke-width:1} svg.map .flood{fill:none;stroke:#c62828;stroke-width:2.5}
svg.map .gage rect{fill:#90caf9;stroke:#1565c0} svg.map .gage line{stroke:#1565c0;stroke-width:1.5}
svg.plot{width:100%;max-width:980px;height:auto;display:block;margin:8px 0;font-size:11px}
svg.plot text{fill:var(--muted)} svg.plot .lab{fill:var(--ink)} svg.plot .gl{stroke:var(--line)} svg.plot .frame{fill:none;stroke:var(--line)}
svg.plot .bar{fill:#5c9bd5} svg.plot .pipe{fill:#d7ccc8;stroke:#6d4c41;stroke-width:1} svg.plot .mh{stroke:#6d4c41;stroke-width:3}
svg.plot .nlab{fill:var(--ink);font-size:10px}
.legend{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:12px;color:var(--muted);margin-bottom:8px}
.legend i{display:inline-block;vertical-align:middle;margin-right:5px}
.legend .mk{width:10px;height:10px} .legend .mk.flood{border:2px solid #c62828;border-radius:50%;background:none}
.legend .ln{width:22px;height:0;border-top:3px solid} .legend .sw{width:110px;height:10px;border-radius:2px}
td.msg{white-space:normal}
.desc{white-space:pre-wrap;color:var(--muted);font-size:13px;max-height:260px;overflow:auto;border:1px solid var(--line);border-radius:6px;padding:8px 12px;margin:6px 0}
"""


def _rows_html(cols, rows, limit=400):
    body = _table(cols, [[_e(c) for c in r] for r in rows[:limit]])
    if len(rows) > limit:
        body += f'<p class="muted">{len(rows) - limit} more rows not shown.</p>'
    return body


def _model_html(cfg, mi):
    pid = lambda k: f"m{mi}-{k}"
    present, out = [], []

    def add(key, title, body, note="", open_=True):
        if body:
            present.append((pid(key), title))
            out.append(_blk(_e(title), body, note, open_, pid(key)))

    if cfg.get("only_out"):
        add("out", "Binary results", _out_html(cfg))
        return "".join(out), present
    m, U = cfg["model"], cfg["model"]["units"]
    if m["title"]:
        out.append(f'<div class="desc">{_e(chr(10).join(m["title"]))}</div>')
    add("sum", "Summary", _two(list(cfg["summary"].items())))
    add("map", "Network map", svg_map(cfg), f"{len(m['coords'])} nodes with coordinates")
    add("opt", "Options", _table(["Option", "Value", "Meaning"],
                                 [[f"<b>{_e(k)}</b>", _e(v), (_e(OPTION_HELP.get(k, "")), "msg")] for k, v in m["options"].items()])
        + (_kv(list(m["evaporation"].items())) and "<h4>Evaporation</h4>" + _kv(list(m["evaporation"].items()))))
    if cfg["rain"]:
        body = ""
        for r in cfg["rain"]:
            body += f"<h4>Rain gage {_e(r['name'])}</h4>" + _kv([(k, r.get(k)) for k in (
                "format", "interval", "scf", "source", "data", "total_depth", "peak", "duration", "starts", "note", "problem")])
            body += svg_hyeto(r, U)
        add("rain", "Rainfall", body)
    groups = [("Subcatchments", ("SUBCATCHMENTS", "SUBAREAS", "INFILTRATION", "LID_USAGE", "COVERAGES", "GROUNDWATER")),
              ("Nodes", ("JUNCTIONS", "OUTFALLS", "STORAGE", "DIVIDERS", "INFLOWS", "DWF", "RDII", "TREATMENT")),
              ("Links", ("CONDUITS", "PUMPS", "ORIFICES", "WEIRS", "OUTLETS", "XSECTIONS", "LOSSES")),
              ("Water quality", ("POLLUTANTS", "LANDUSES", "BUILDUP", "WASHOFF", "LOADINGS"))]
    for title, secs in groups:
        body, n = "", 0
        for sec in secs:
            cols, rows = section_rows(cfg, sec)
            if rows:
                n += len(rows)
                body += _blk(f"[{sec}]", _rows_html(cols, rows), f"{len(rows)} rows", open_=sec in secs[:1])
        if title == "Links" and cfg["conduit_geometry"]:
            cg = cfg["conduit_geometry"]
            body += _blk("Conduit geometry (computed from node inverts and offsets)", _rows_html(
                ["Conduit", "From", "To", f"Length ({U['length']})", "Shape", "Geom1", "Inlet z", "Outlet z", "Slope %"],
                [[c["name"], c["from"], c["to"], g(c["length"]), c["shape"], g(c["geom1"]), g(c.get("z_up")),
                  g(c.get("z_dn")), g(c.get("slope_pct"))] for c in cg]), open_=False)
        add(title.lower().replace(" ", ""), title, body, f"{n} rows")
    body = ""
    if m["lid_controls"]:
        body += "<h4>LID controls</h4>" + _table(["Name", "Type", "Layers"], [[_e(n), _e(f"{d['type']} ({LID_TYPES.get(d['type'], '?')})"),
                                                                            (_e("; ".join(f"{k}: {v}" for k, v in d["layers"].items())), "msg")]
                                                                           for n, d in m["lid_controls"].items()])
    if m["curves"]:
        body += "<h4>Curves</h4>" + _table(["Name", "Type", "Points", "X range", "Y range"],
                                           [[_e(n), _e(c["type"]), _e(len(c["points"])), _e(rng([p[0] for p in c["points"]])),
                                             _e(rng([p[1] for p in c["points"]]))] for n, c in m["curves"].items()])
    if m["timeseries"]:
        body += "<h4>Time series</h4>" + _table(["Name", "Points / file", "Time range (h)", "Values", "Starts"],
                                                [[_e(n), _e(t["file"] or len(t["points"])),
                                                  _e(rng([p[0] for p in t["points"]])), _e(rng([p[1] for p in t["points"]])),
                                                  _e(t.get("start", ""))] for n, t in m["timeseries"].items()])
    if m["patterns"]:
        body += "<h4>Patterns</h4>" + _table(["Name", "Type", "Factors", "Range"],
                                             [[_e(n), _e(p["type"]), _e(len(p["factors"])), _e(rng(p["factors"]))] for n, p in m["patterns"].items()])
    if m["controls"]:
        body += f'<h4>Control rules ({len(m["rules"])})</h4><pre class="desc">{_e(chr(10).join(m["controls"]))}</pre>'
    add("data", "LIDs, curves, time series, patterns & controls", body, open_=False)
    if cfg["profiles"]:
        add("prof", "Profiles", "".join(svg_profile(n, d, U) for n, d in cfg["profiles"].items()),
            "with maximum HGL" if cfg.get("out") and not cfg["out"].get("error") else "")
    add("secs", "Sections in the file", _table(["Section", "Rows"], [[_e(s["section"]), _e(s["rows"])] for s in cfg["sections"]]),
        f"{len(cfg['sections'])} sections", open_=False)
    add("rpt", "Report file (.rpt)", _rpt_html(cfg))
    add("out", "Binary results (.out)", _out_html(cfg))
    return "".join(out), present


def _rpt_html(cfg):
    r = cfg.get("rpt")
    if not r:
        return '<p class="muted">No report file found (model not run, or saved elsewhere - use --rpt).</p>'
    body = _two([(k, r.get(k)) for k in ("file", "engine", "analysis begun", "analysis ended", "elapsed time")]
                + [(f"{k} error %", ", ".join(g(x) for x in v["continuity_error_pct"])) for k, v in r.get("continuity", {}).items()]
                + list((r.get("routing_time_step") or {}).items()))
    if r["errors"] or r["warnings"]:
        body += _table(["Level", "Message"], [['<span class="warn">ERROR</span>', (_e(x), "msg")] for x in r["errors"]]
                       + [["warning", (_e(x), "msg")] for x in r["warnings"]], sortable=False)
    for k, v in r.get("continuity", {}).items():
        body += _blk(_e(k), _kv([(a, " / ".join(b)) for a, b in v["terms"].items()]), open_=False)
    for k, lab, cols in (("flooding", "Node flooding", ["Node", "Hours flooded", "Max rate", "Day", "hh:mm", "Total volume", "Max ponded"]),
                         ("node_surcharge", "Node surcharge", ["Node", "Type", "Hours", "Max height above crown", "Min depth below rim"]),
                         ("conduit_surcharge", "Conduit surcharge", ["Conduit", "Both ends full", "Upstream full", "Dnstream full",
                                                                     "Above normal flow", "Capacity limited"]),
                         ("outfalls", "Outfall loading", ["Outfall", "Flow freq %", "Avg flow", "Max flow", "Total volume"])):
        if k in r:
            body += _blk(lab, _rows_html(cols, r[k]) if r[k] else '<p class="muted">none</p>', f"{len(r[k])} rows", open_=bool(r[k]) and k == "flooding")
    for k in ("highest_continuity_errors", "time_step_critical", "flow_instability", "nonconverging_nodes"):
        if r.get(k):
            body += f"<h4>{_e(k.replace('_', ' '))}</h4><pre class=\"desc\">{_e(chr(10).join(r[k]))}</pre>"
    return body


def _out_html(cfg):
    o = cfg.get("out")
    if not o:
        return '<p class="muted">No binary results file found.</p>'
    if o.get("error"):
        return f'<p class="warn">{_e(o["file"])}: {_e(o["error"])}</p>'
    U = cfg["model"]["units"] if cfg.get("model") else units_for(str(o.get("flow_units", "CFS")))
    body = _two([(k, v) for k, v in o.items() if not k.startswith("_")])
    body += svg_system(o, U)
    if cfg.get("model"):
        body += "<h4>Map - maximum conduit capacity used and flooded nodes</h4>" + svg_map(cfg, results=True)
    nm = sorted(o["_node_max"].items(), key=lambda kv: (-kv[1][4], -kv[1][0]))
    body += _blk("Node maxima", _rows_html(["Node", "Max depth", "Time", "Max head", "Max total inflow", "Flood volume", "Hours flooded"],
                 [[n, g(v[0]), swmm_datetime(v[1]).strftime("%m-%d %H:%M") if v[1] else "", g(v[2]), g(v[3]), g(v[4]), g(v[5] / 3600)]
                  for n, v in nm]), f"{len(nm)} nodes", open_=False)
    lk = sorted(o["_link_max"].items(), key=lambda kv: -kv[1][3])
    body += _blk("Link maxima", _rows_html(["Link", "Max |flow|", "Time", "Max velocity", "Max capacity"],
                 [[n, g(v[0]), swmm_datetime(v[1]).strftime("%m-%d %H:%M") if v[1] else "", g(v[2]), g(v[3])] for n, v in lk]),
                 f"{len(lk)} links", open_=False)
    sm = sorted(o["_sub_max"].items(), key=lambda kv: -kv[1][0])
    if sm:
        body += _blk("Subcatchment peak runoff", _rows_html(["Subcatchment", "Peak runoff", "Time"],
                     [[n, g(v[0]), swmm_datetime(v[1]).strftime("%m-%d %H:%M") if v[1] else ""] for n, v in sm]), open_=False)
    return body


def render_html(cfg) -> str:
    nav, body = [], []
    for i, m in enumerate(cfg["models"]):
        sec = [f'<section id="m{i}"><h2>{_e(m["name"])}']
        if m.get("only_out"):
            sec.append(' <span class="badge o">results file only</span>')
        elif m.get("rpt") or m.get("out"):
            ok = not (m.get("rpt") or {}).get("errors") and not (m.get("out") or {}).get("error_code")
            sec.append(f' <span class="badge{"" if ok else " w"}">{"run finished" if ok else "run errors"}</span>')
        else:
            sec.append(' <span class="badge o">not run</span>')
        if m["warnings"]:
            sec.append(f' <span class="badge w">{len(m["warnings"])} warning(s)</span>')
        sec.append(f'</h2><p class="muted"><code>{_e(m["path"])}</code></p>')
        if m["warnings"]:
            sec.append('<ul class="warns">' + "".join(f"<li>{_e(w)}</li>" for w in m["warnings"]) + "</ul>")
        html_, present = _model_html(m, i)
        sec.append(html_ + "</section>")
        body.append("".join(sec))
        nav.append(f'<a href="#m{i}"><b>{_e(m["name"])}</b></a>')
        nav += [f'<a class="sub" href="#{b}">{_e(t)}</a>' for b, t in present]
    m0 = cfg["models"][0]
    sub = (f'{len(cfg["models"])} models' if len(cfg["models"]) > 1 else
           " · ".join(x for x in (m0["summary"].get("Units", "").split(" - ")[0], m0["summary"].get("Flow routing", "").split(",")[0],
                                  m0["summary"].get("Simulation period")) if x))
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


# --------------------------------------------------------------------------- #
def main(argv=None, prog=None):
    ap = argparse.ArgumentParser(prog=prog, description="Show the configuration of an EPA SWMM 5 model.")
    ap.add_argument("model", help=".inp file, a folder of .inp files, or a binary .out file")
    ap.add_argument("--rpt", help="report file (default: <model>.rpt)")
    ap.add_argument("--out", help="binary results file (default: <model>.out)")
    ap.add_argument("--all", action="store_true", help="terminal: every row of every section, and the control rules")
    ap.add_argument("--json", metavar="FILE", help="also write the configuration to a JSON file")
    ap.add_argument("--html", nargs="?", const="", metavar="FILE",
                    help="write the configuration as a webpage and open it in the browser, instead of printing it "
                         "(default file: <model>_config.html next to the .inp)")
    ap.add_argument("--no-browser", action="store_true", help="with --html, write the page but don't open it")
    a = ap.parse_args(argv)

    cfg = build_config(a.model, a.rpt, a.out)
    if a.html is None:
        show(cfg, a.all)
    else:
        m0 = cfg["models"][0]
        out = Path(a.html) if a.html else (Path(m0["folder"]) / f"{m0['name']}_config.html" if len(cfg["models"]) == 1
                                           else Path(cfg["root"]) / f"{cfg['name']}_config.html")
        write_html(cfg, out, open_browser=not a.no_browser)
    if a.json:
        Path(a.json).write_text(json.dumps(strip_private(cfg), indent=2, default=str), encoding="utf-8")
        print(f"Configuration written to {a.json}")


if __name__ == "__main__":
    main()
