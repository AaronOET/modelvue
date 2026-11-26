"""
modelvue.epanetvue - Read an EPANET 2.x water distribution model and show its configuration.

Parses the input file (.inp) following the EPANET 2.2 input format (Users Manual,
Appendix C), checks the network (references, connectivity, supply sources,
valve placement, curves and patterns), and - if the model has been run - reads
the report file (.rpt: warnings, errors, hydraulic status) and the binary
output file (.out, layout per the OWA-EPANET "Output File" reference: prolog,
energy usage, dynamic results, epilog) for pressures, velocities, tank levels,
pump energy and water quality.

No extra packages needed (standard library only).

Usage
-----
    epanetvue  network.inp                    # report in the terminal
    epanetvue  folder                         # every .inp in the folder (and below)
    epanetvue  network.inp --all              # every row of every section
    epanetvue  network.inp --rpt x.rpt --out x.out
    epanetvue  results.out                    # only a binary output file
    epanetvue  network.inp --json config.json
    epanetvue  network.inp --html [report.html] [--no-browser]

EPANET-GUI project files (.net) are a private binary format: open them in EPANET
and use File > Export > Network to write an .inp.
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

from ._common import _blk, _e, h1, h2, row, strip_private, _table, _ticks, _two
from ._inp import num, table

# --------------------------------------------------------------------------- #
# Format definitions (EPANET 2.2 / OWA-EPANET 2.3)
# --------------------------------------------------------------------------- #
SECTIONS = ["TITLE", "JUNCTIONS", "RESERVOIRS", "TANKS", "PIPES", "PUMPS", "VALVES", "TAGS", "DEMANDS", "STATUS",
            "PATTERNS", "CURVES", "CONTROLS", "RULES", "ENERGY", "EMITTERS", "LEAKAGE", "QUALITY", "SOURCES",
            "REACTIONS", "MIXING", "TIMES", "REPORT", "OPTIONS", "COORDINATES", "VERTICES", "LABELS", "BACKDROP", "END"]
RAW_SECTIONS = {"TITLE", "CONTROLS", "RULES", "LABELS", "BACKDROP", "TAGS"}

US_FLOW = ["CFS", "GPM", "MGD", "IMGD", "AFD"]
SI_FLOW = ["LPS", "LPM", "MLD", "CMH", "CMD", "CMS"]
FLOW_CODES = ["CFS", "GPM", "MGD", "IMGD", "AFD", "LPS", "LPM", "MLD", "CMH", "CMD", "CMS"]          # EN_FlowUnits
LINK_CODES = ["CV pipe", "pipe", "pump", "PRV", "PSV", "PBV", "FCV", "TCV", "GPV", "PCV"]          # EN_LinkType
QUAL_CODES = ["none", "chemical", "age", "trace"]                                               # EN_QualityType
STAT_CODES = ["none (time series)", "average", "minimum", "maximum", "range"]                    # EN_StatisticType
PRESS_CODES = ["psi", "kPa", "m", "bar", "ft"]
STATUS_CODES = {0: "closed (max head exceeded)", 1: "temporarily closed", 2: "closed", 3: "open", 4: "active",
                5: "open (max flow exceeded)", 6: "open (FCV can't supply flow)", 7: "open (can't meet pressure)"}

COLS = {
    "JUNCTIONS": ["ID", "Elevation", "Demand", "Pattern"],
    "RESERVOIRS": ["ID", "Head", "Pattern"],
    "TANKS": ["ID", "Elevation", "InitLevel", "MinLevel", "MaxLevel", "Diameter", "MinVol", "VolCurve", "Overflow"],
    "PIPES": ["ID", "Node1", "Node2", "Length", "Diameter", "Roughness", "MinorLoss", "Status"],
    "PUMPS": ["ID", "Node1", "Node2", "Parameters"],
    "VALVES": ["ID", "Node1", "Node2", "Diameter", "Type", "Setting", "MinorLoss"],
    "DEMANDS": ["Junction", "Demand", "Pattern", "Category"],
    "STATUS": ["ID", "Status/Setting"],
    "EMITTERS": ["Junction", "Coefficient"],
    "LEAKAGE": ["Pipe", "C1", "C2"],
    "QUALITY": ["Node", "InitQual"],
    "SOURCES": ["Node", "Type", "Strength", "Pattern"],
    "MIXING": ["Tank", "Model", "Fraction"],
    "COORDINATES": ["Node", "X", "Y"], "VERTICES": ["Link", "X", "Y"],
}

OPTION_KEYS = ["UNITS", "PRESSURE", "HEADLOSS", "HYDRAULICS", "QUALITY", "VISCOSITY", "DIFFUSIVITY", "SPECIFIC GRAVITY",
               "TRIALS", "ACCURACY", "HEADERROR", "FLOWCHANGE", "CHECKFREQ", "MAXCHECK", "DAMPLIMIT", "UNBALANCED",
               "PATTERN", "DEMAND MODEL", "MINIMUM PRESSURE", "REQUIRED PRESSURE", "PRESSURE EXPONENT",
               "DEMAND MULTIPLIER", "EMITTER EXPONENT", "EMITTER BACKFLOW", "TOLERANCE", "MAP", "QUALITY TOLERANCE"]
OPTION_DEFAULTS = {"UNITS": "GPM", "HEADLOSS": "H-W", "QUALITY": "NONE", "VISCOSITY": "1.0", "DIFFUSIVITY": "1.0",
                   "TRIALS": "200", "ACCURACY": "0.001", "HEADERROR": "0", "FLOWCHANGE": "0", "UNBALANCED": "STOP",
                   "PATTERN": "1", "DEMAND MULTIPLIER": "1.0", "DEMAND MODEL": "DDA", "MINIMUM PRESSURE": "0",
                   "REQUIRED PRESSURE": "0.1", "PRESSURE EXPONENT": "0.5", "EMITTER EXPONENT": "0.5", "TOLERANCE": "0.01"}
OPTION_HELP = {"UNITS": "flow units (CFS, GPM, MGD, IMGD, AFD = US; LPS, LPM, MLD, CMH, CMD, CMS = SI)",
               "HEADLOSS": "H-W Hazen-Williams, D-W Darcy-Weisbach, C-M Chezy-Manning",
               "QUALITY": "NONE, CHEMICAL name units, AGE, or TRACE node", "UNBALANCED": "what to do if not converged",
               "PATTERN": "default demand pattern", "DEMAND MODEL": "DDA demand-driven / PDA pressure-driven",
               "MINIMUM PRESSURE": "PDA: no demand below this", "REQUIRED PRESSURE": "PDA: full demand above this",
               "PRESSURE EXPONENT": "PDA exponent", "TRIALS": "max hydraulic trials", "ACCURACY": "convergence (flow change ratio)",
               "HEADERROR": "max head-loss error (0 = off)", "FLOWCHANGE": "max flow change (0 = off)",
               "TOLERANCE": "water quality tolerance", "SPECIFIC GRAVITY": "relative to water at 4 C",
               "VISCOSITY": "relative kinematic viscosity", "DIFFUSIVITY": "relative chemical diffusivity",
               "DEMAND MULTIPLIER": "scales all base demands", "EMITTER EXPONENT": "emitter pressure exponent",
               "HYDRAULICS": "USE/SAVE a hydraulics file", "CHECKFREQ": "status check frequency", "MAXCHECK": "",
               "DAMPLIMIT": "", "PRESSURE": "pressure units (2.3)"}
TIME_KEYS = ["DURATION", "HYDRAULIC TIMESTEP", "QUALITY TIMESTEP", "RULE TIMESTEP", "PATTERN TIMESTEP", "PATTERN START",
             "REPORT TIMESTEP", "REPORT START", "START CLOCKTIME", "STATISTIC"]
TIME_DEFAULTS = {"DURATION": "0", "HYDRAULIC TIMESTEP": "1:00", "PATTERN TIMESTEP": "1:00", "PATTERN START": "0",
                 "REPORT TIMESTEP": "1:00", "REPORT START": "0", "START CLOCKTIME": "12 AM", "STATISTIC": "NONE"}


def units_for(flow: str, pressure: str | None = None) -> dict:
    us = flow.upper() in US_FLOW
    return {"system": "US" if us else "SI", "flow": flow.upper(), "length": "ft" if us else "m",
            "diameter": "in" if us else "mm", "pressure": (pressure or ("psi" if us else "m")).lower(),
            "head": "ft" if us else "m", "velocity": "ft/s" if us else "m/s", "volume": "ft3" if us else "m3"}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def g(v, nd=4):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "" if v is None else str(v)
    if math.isnan(f) or math.isinf(f):
        return str(f)
    if f.is_integer() and abs(f) < 1e12:
        return str(int(f))
    return f"{f:.{nd}g}"


def rng(vals, unit=""):
    v = [x for x in vals if x is not None]
    if not v:
        return "-"
    lo, hi = min(v), max(v)
    return f"{g(lo)}{unit}" if lo == hi else f"{g(lo)} .. {g(hi)}{unit}"


def parse_seconds(txt: str, default_unit="HOURS") -> float | None:
    """EPANET time: '1:30', '1:30:00', '24', '24 HOURS', '30 MIN', '12 AM', '1:30 PM'."""
    if txt is None:
        return None
    t = txt.strip().upper().split()
    if not t:
        return None
    ampm = t[-1] if t[-1] in ("AM", "PM") else None
    unit = t[1] if len(t) > 1 and t[1] not in ("AM", "PM") else default_unit
    v = t[0]
    if ":" in v:
        p = [num(x, 0) for x in v.split(":")]
        s = p[0] * 3600 + (p[1] if len(p) > 1 else 0) * 60 + (p[2] if len(p) > 2 else 0)
    else:
        x = num(v)
        if x is None:
            return None
        mult = {"SEC": 1, "SECONDS": 1, "MIN": 60, "MINUTES": 60, "HOUR": 3600, "HOURS": 3600, "DAY": 86400, "DAYS": 86400}
        s = x * mult.get(unit, 3600)
    if ampm:
        h = s % 43200
        s = h + (43200 if ampm == "PM" else 0)
    return s


def hhmm(sec):
    if sec is None:
        return "-"
    sec = int(round(sec))
    return f"{sec // 3600}:{(sec % 3600) // 60:02d}" + (f":{sec % 60:02d}" if sec % 60 else "")


def tokenize(line: str) -> list[str]:
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


def keyed(tokens: list[str], keys: list[str]):
    """Match a multi-word keyword at the start of a line: returns (KEY, rest)."""
    up = [t.upper() for t in tokens]
    for k in sorted(keys, key=lambda k: -len(k.split())):
        kw = k.split()
        if up[:len(kw)] == kw:
            return k, tokens[len(kw):]
    return up[0], tokens[1:]


# --------------------------------------------------------------------------- #
# INP file
# --------------------------------------------------------------------------- #
def read_inp(path: Path) -> dict:
    head = path.read_bytes()[:16]
    if b"<EPANET2>" in head or b"<EPANET" in head:
        raise ValueError("this is an EPANET-GUI project file (.net, binary); open it in EPANET and use "
                         "File > Export > Network to write an .inp")
    secs, raw, order, cur, comments = {}, {}, [], None, {}
    last_comment = None
    for ln in path.read_text(errors="replace").splitlines():
        s = ln.strip()
        if s.startswith("[") and "]" in s:
            cur = s[1:s.index("]")].strip().upper()
            secs.setdefault(cur, [])
            raw.setdefault(cur, [])
            if cur not in order and cur != "END":
                order.append(cur)
            last_comment = None
            continue
        if cur is None or cur == "END":
            continue
        if cur in RAW_SECTIONS:
            if s and not s.startswith(";"):
                raw[cur].append(ln.rstrip())
            continue
        if s.startswith(";"):
            last_comment = s[1:].strip()
            continue
        t = tokenize(ln)
        if t:
            secs[cur].append(t)
            if cur == "CURVES" and last_comment is not None:
                comments.setdefault(t[0], last_comment)
            last_comment = None
    return {"sections": secs, "raw": raw, "order": order, "curve_comments": comments}


def rows_as_dicts(rows, cols):
    out = []
    for t in rows:
        d = {c: (t[i] if i < len(t) else "") for i, c in enumerate(cols)}
        if len(t) > len(cols):
            d[cols[-1]] = " ".join(t[len(cols) - 1:])
        out.append(d)
    return out


def build_model(inp: dict) -> dict:
    opt = {}
    for t in table(inp, "OPTIONS"):
        k, rest = keyed(t, OPTION_KEYS)
        opt[k] = " ".join(rest)
    times = {}
    for t in table(inp, "TIMES"):
        k, rest = keyed(t, TIME_KEYS)
        times[k] = " ".join(rest)
    U = units_for(opt.get("UNITS", "GPM"), opt.get("PRESSURE"))
    m = {"options": opt, "times": times, "units": U}
    m["title"] = [ln.strip() for ln in inp["raw"].get("TITLE", []) if ln.strip()]
    for sec in ("JUNCTIONS", "RESERVOIRS", "TANKS", "PIPES", "VALVES", "EMITTERS", "LEAKAGE", "QUALITY", "SOURCES", "MIXING",
                "STATUS"):
        m[sec.lower()] = rows_as_dicts(table(inp, sec), COLS[sec])
    pumps = []
    for t in table(inp, "PUMPS"):
        p = {"ID": t[0], "Node1": t[1] if len(t) > 1 else "", "Node2": t[2] if len(t) > 2 else "", "params": {}}
        kv = t[3:]
        for i in range(0, len(kv) - 1, 2):
            p["params"][kv[i].upper()] = kv[i + 1]
        if len(kv) == 1 and num(kv[0]) is not None:      # very old format: power only
            p["params"]["POWER"] = kv[0]
        pumps.append(p)
    m["pumps"] = pumps
    dem = {}
    for t in table(inp, "DEMANDS"):
        dem.setdefault(t[0], []).append({"demand": num(t[1], 0.0) if len(t) > 1 else 0.0,
                                         "pattern": t[2] if len(t) > 2 else "", "category": " ".join(t[3:])})
    m["demands"] = dem
    pats = {}
    for t in table(inp, "PATTERNS"):
        pats.setdefault(t[0], []).extend(v for v in (num(x) for x in t[1:]) if v is not None)
    m["patterns"] = pats
    curves = {}
    for t in table(inp, "CURVES"):
        c = curves.setdefault(t[0], {"type": "", "points": []})
        vals = [num(x) for x in t[1:]]
        for i in range(0, len(vals) - 1, 2):
            if vals[i] is not None and vals[i + 1] is not None:
                c["points"].append((vals[i], vals[i + 1]))
    for cid, com in inp["curve_comments"].items():
        mm = re.match(r"\s*(PUMP|EFFICIENCY|VOLUME|HEADLOSS|GENERAL|VALVE)\s*:", com, re.I)
        if cid in curves and mm:
            curves[cid]["type"] = mm.group(1).upper()
        if cid in curves:
            curves[cid]["comment"] = com
    m["curves"] = curves
    m["energy"] = [" ".join(t) for t in table(inp, "ENERGY")]
    m["reactions"] = [" ".join(t) for t in table(inp, "REACTIONS")]
    m["report"] = {(t[0].upper()): " ".join(t[1:]) for t in table(inp, "REPORT")}
    m["controls"] = [ln.strip() for ln in inp["raw"].get("CONTROLS", []) if ln.strip()]
    m["rules_text"] = [ln.rstrip() for ln in inp["raw"].get("RULES", []) if ln.strip()]
    m["rules"] = [ln.split(None, 1)[1].strip() if len(ln.split(None, 1)) > 1 else "" for ln in m["rules_text"]
                  if ln.strip().upper().startswith("RULE")]
    m["coords"] = {t[0]: (num(t[1]), num(t[2])) for t in table(inp, "COORDINATES") if len(t) >= 3}
    vert = {}
    for t in table(inp, "VERTICES"):
        if len(t) >= 3:
            vert.setdefault(t[0], []).append((num(t[1]), num(t[2])))
    m["vertices"] = vert
    labels = []
    for ln in inp["raw"].get("LABELS", []):
        t = tokenize(ln)
        if len(t) >= 3 and num(t[0]) is not None:
            labels.append((num(t[0]), num(t[1]), t[2]))
    m["labels"] = labels
    return m


# --------------------------------------------------------------------------- #
# Derived data & checks
# --------------------------------------------------------------------------- #
def node_table(m) -> dict:
    nodes = {}
    for r in m["junctions"]:
        nodes[r["ID"]] = {"type": "junction", "elev": num(r["Elevation"])}
    for r in m["reservoirs"]:
        nodes[r["ID"]] = {"type": "reservoir", "elev": num(r["Head"])}
    for r in m["tanks"]:
        nodes[r["ID"]] = {"type": "tank", "elev": num(r["Elevation"])}
    return nodes


def link_table(m) -> dict:
    links = {}
    for r in m["pipes"]:
        links[r["ID"]] = {"type": "CV pipe" if r["Status"].upper() == "CV" else "pipe", "n1": r["Node1"], "n2": r["Node2"],
                          "length": num(r["Length"]), "diameter": num(r["Diameter"])}
    for p in m["pumps"]:
        links[p["ID"]] = {"type": "pump", "n1": p["Node1"], "n2": p["Node2"]}
    for r in m["valves"]:
        links[r["ID"]] = {"type": r["Type"].upper(), "n1": r["Node1"], "n2": r["Node2"], "diameter": num(r["Diameter"])}
    return links


def junction_demands(m) -> dict:
    """Base demands per junction: [DEMANDS] replaces the [JUNCTIONS] demand (Users Manual C.3)."""
    default_pat = m["options"].get("PATTERN", "1")
    out = {}
    for r in m["junctions"]:
        if r["ID"] in m["demands"]:
            out[r["ID"]] = [(d["demand"], d["pattern"] or default_pat, d["category"]) for d in m["demands"][r["ID"]]]
        else:
            out[r["ID"]] = [(num(r["Demand"], 0.0), r["Pattern"] or default_pat, "")]
    return out


def pump_curve(m, p) -> dict:
    e = {"id": p["ID"], "from": p["Node1"], "to": p["Node2"]}
    prm = p["params"]
    if "POWER" in prm:
        e["kind"] = f"constant power {prm['POWER']} {'hp' if m['units']['system'] == 'US' else 'kW'}"
    if "HEAD" in prm:
        c = m["curves"].get(prm["HEAD"])
        e["curve"] = prm["HEAD"]
        if c is None:
            e["problem"] = f"head curve '{prm['HEAD']}' not defined"
        else:
            pts = c["points"]
            if len(pts) == 1:
                q, h = pts[0]
                e["kind"] = "single-point curve"
                e["fitted"] = "EPANET fits h = A - B q^2: shutoff head 1.33 x design head, max flow 2 x design flow"
                e["_curve"] = [(x * 2 * q / 40, 1.33333 * h - (1.33333 * h - h) * (x * 2 * q / 40 / q) ** 2) for x in range(41)]
                e["design"] = f"{g(q)} {m['units']['flow']} at {g(h)} {m['units']['head']}"
            elif len(pts) == 3:
                e["kind"] = "three-point curve (fitted to h = A - B q^C)"
                e["_curve"] = pts
            else:
                e["kind"] = f"multi-point curve ({len(pts)} points, linear interpolation)"
                e["_curve"] = pts
            if len(pts) > 1 and any(pts[i + 1][1] >= pts[i][1] for i in range(len(pts) - 1)):
                e["problem"] = "head does not decrease with flow (EPANET error 227)"
            e["_points"] = pts
    if "SPEED" in prm:
        e["speed"] = prm["SPEED"]
    if "PATTERN" in prm:
        e["speed pattern"] = prm["PATTERN"]
    return e


def components(nodes, links) -> list[set]:
    adj = {n: set() for n in nodes}
    for l in links.values():
        if l["n1"] in adj and l["n2"] in adj:
            adj[l["n1"]].add(l["n2"])
            adj[l["n2"]].add(l["n1"])
    seen, comps = set(), []
    for n in nodes:
        if n in seen:
            continue
        stack, comp = [n], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(adj[x] - comp)
        seen |= comp
        comps.append(comp)
    return comps


def checks(m, nodes, links, dem) -> list[str]:
    W = []
    ids = {}
    for sec in ("junctions", "reservoirs", "tanks"):
        for r in m[sec]:
            ids.setdefault(r["ID"], []).append(sec)
    dup = [k for k, v in ids.items() if len(v) > 1]
    if dup:
        W.append(f"duplicate node IDs (EPANET error 215): {', '.join(dup[:10])}")
    lid = {}
    for sec in ("pipes", "valves"):
        for r in m[sec]:
            lid.setdefault(r["ID"], []).append(sec)
    for p in m["pumps"]:
        lid.setdefault(p["ID"], []).append("pumps")
    dup = [k for k, v in lid.items() if len(v) > 1]
    if dup:
        W.append(f"duplicate link IDs (EPANET error 215): {', '.join(dup[:10])}")
    upper_nodes = {n.upper(): n for n in nodes}
    for lk, l in links.items():
        for end in ("n1", "n2"):
            if l[end] not in nodes:
                hint = f" (a node '{upper_nodes[l[end].upper()]}' exists - IDs are case sensitive?)" if l[end].upper() in upper_nodes else ""
                W.append(f"{l['type']} {lk}: node '{l[end]}' is not defined{hint}")
        if l["n1"] == l["n2"]:
            W.append(f"{l['type']} {lk}: both ends at the same node")
    used = {l["n1"] for l in links.values()} | {l["n2"] for l in links.values()}
    orphan = [n for n in nodes if n not in used]
    if orphan:
        W.append(f"{len(orphan)} node(s) not connected to any link (EPANET error 233): {', '.join(orphan[:10])}")
    sources = {n for n, v in nodes.items() if v["type"] in ("reservoir", "tank")}
    if nodes and not sources:
        W.append("no reservoir or tank - the network has no fixed-grade node (EPANET error 223)")
    for comp in components(nodes, links):
        if len(comp) > 1 and not comp & sources:
            W.append(f"{len(comp)} connected node(s) with no reservoir or tank, e.g. {', '.join(sorted(comp)[:5])}")
    # valves next to tanks/reservoirs (error 219) and valves in series of the same kind (error 220)
    for r in m["valves"]:
        t = r["Type"].upper()
        if t in ("PRV", "PSV", "FCV"):
            for end in ("Node1", "Node2"):
                n = nodes.get(r[end])
                if n and n["type"] in ("reservoir", "tank"):
                    W.append(f"{t} {r['ID']} is connected directly to {n['type']} {r[end]} (EPANET error 219)")
    # references
    pats = set(m["patterns"])
    for j, ds in dem.items():
        for d, p, _ in ds:
            if d and p and p not in pats and p != "1":
                W.append(f"junction {j}: demand pattern '{p}' is not defined")
    dp = m["options"].get("PATTERN")
    if dp and dp not in pats:
        W.append(f"[OPTIONS] PATTERN '{dp}' is not defined - demands without a pattern stay constant")
    for r in m["reservoirs"]:
        if r["Pattern"] and r["Pattern"] not in pats:
            W.append(f"reservoir {r['ID']}: head pattern '{r['Pattern']}' is not defined")
    for r in m["tanks"]:
        lv = [num(r[k]) for k in ("MinLevel", "InitLevel", "MaxLevel")]
        if None not in lv and not (lv[0] <= lv[1] <= lv[2]):
            W.append(f"tank {r['ID']}: initial level {g(lv[1])} not between min {g(lv[0])} and max {g(lv[2])}")
        if r["VolCurve"] and r["VolCurve"] not in ("*",) and r["VolCurve"] not in m["curves"]:
            W.append(f"tank {r['ID']}: volume curve '{r['VolCurve']}' is not defined")
    for p in m["pumps"]:
        e = pump_curve(m, p)
        if e.get("problem"):
            W.append(f"pump {p['ID']}: {e['problem']}")
        if "PATTERN" in p["params"] and p["params"]["PATTERN"] not in pats:
            W.append(f"pump {p['ID']}: speed pattern '{p['params']['PATTERN']}' is not defined")
        if "POWER" not in p["params"] and "HEAD" not in p["params"]:
            W.append(f"pump {p['ID']}: neither POWER nor HEAD given")
    for r in m["valves"]:
        if r["Type"].upper() == "GPV" and r["Setting"] not in m["curves"]:
            W.append(f"GPV {r['ID']}: head-loss curve '{r['Setting']}' is not defined")
    odd_c, bad_ld = [], []
    hl = m["options"].get("HEADLOSS", "H-W").upper()
    for r in m["pipes"]:
        L, D, C = num(r["Length"]), num(r["Diameter"]), num(r["Roughness"])
        if (L is not None and L <= 0) or (D is not None and D <= 0):
            bad_ld.append(r["ID"])
        if C is not None and hl == "H-W" and not (40 <= C <= 160):
            odd_c.append(f"{r['ID']} ({g(C)})")
    if bad_ld:
        W.append(f"{len(bad_ld)} pipe(s) with non-positive length or diameter: {', '.join(bad_ld[:8])}")
    if odd_c:
        W.append(f"{len(odd_c)} pipe(s) with Hazen-Williams C outside 40-160 (often deliberate for short connector pipes): "
                 f"{', '.join(odd_c[:6])}{' ...' if len(odd_c) > 6 else ''}")
    # controls referencing unknown objects
    for ln in m["controls"]:
        t = ln.split()
        up = [x.upper() for x in t]
        if len(t) >= 2 and up[0] == "LINK" and t[1] not in links:
            W.append(f"control refers to undefined link '{t[1]}': {ln}")
        if "NODE" in up:
            i = up.index("NODE")
            if i + 1 < len(t) and t[i + 1] not in nodes:
                W.append(f"control refers to undefined node '{t[i + 1]}': {ln}")
    for ln in m["rules_text"]:
        t = ln.split()
        up = [x.upper() for x in t]
        for i, w in enumerate(up[:-1]):
            if w in ("NODE", "JUNCTION", "TANK", "RESERVOIR") and up[0] in ("IF", "AND", "OR") and t[i + 1] not in nodes:
                W.append(f"rule refers to undefined node '{t[i + 1]}': {ln.strip()}")
            if w in ("LINK", "PIPE", "PUMP", "VALVE") and up[0] in ("IF", "AND", "OR", "THEN", "ELSE") and t[i + 1] not in links:
                W.append(f"rule refers to undefined link '{t[i + 1]}': {ln.strip()}")
    for src in m["sources"]:
        if src["Node"] not in nodes:
            W.append(f"source at undefined node '{src['Node']}'")
    return W


def summarize(path: Path, inp: dict) -> dict:
    m = build_model(inp)
    U = m["units"]
    nodes, links = node_table(m), link_table(m)
    dem = junction_demands(m)
    mult = num(m["options"].get("DEMAND MULTIPLIER"), 1.0)
    total = sum(d for ds in dem.values() for d, _, _ in ds) * mult
    consume = sum(d for ds in dem.values() for d, _, _ in ds if d > 0) * mult
    inflow = -sum(d for ds in dem.values() for d, _, _ in ds if d < 0) * mult
    pats = m["patterns"]
    # peak demand factor (all junctions on the same pattern only)
    used_pats = {p for ds in dem.values() for d, p, _ in ds if d}
    peak = None
    if len(used_pats) == 1 and next(iter(used_pats)) in pats and pats[next(iter(used_pats))]:
        pf = pats[next(iter(used_pats))]
        peak = f"{g(total * min(pf))} .. {g(total * max(pf))} {U['flow']} (pattern {next(iter(used_pats))}: {g(min(pf))} .. {g(max(pf))})"
    tim = m["times"]
    dur = parse_seconds(tim.get("DURATION", "0"))
    hyd = parse_seconds(tim.get("HYDRAULIC TIMESTEP", "1:00"))
    rep = parse_seconds(tim.get("REPORT TIMESTEP", "1:00"))
    pipes = m["pipes"]
    length = sum(num(r["Length"], 0) for r in pipes)
    by_d = {}
    for r in pipes:
        by_d[num(r["Diameter"])] = by_d.get(num(r["Diameter"]), 0) + num(r["Length"], 0)
    qual = m["options"].get("QUALITY", "NONE")
    tanks = []
    for r in m["tanks"]:
        D = num(r["Diameter"])
        lo, hi = num(r["MinLevel"]), num(r["MaxLevel"])
        vol = math.pi / 4 * D * D * (hi - lo) if None not in (D, lo, hi) and not r["VolCurve"] else None
        tanks.append({"id": r["ID"], "base": num(r["Elevation"]), "levels": f"{r['MinLevel']} / {r['InitLevel']} / {r['MaxLevel']}",
                      "diameter": D, "working volume": vol, "curve": r["VolCurve"],
                      "overflow": r.get("Overflow", "")})
    S = {
        "Title": " / ".join(m["title"][:2]) or None,
        "Units": f"{U['system']} - flow {U['flow']}, pressure {U['pressure']}, length {U['length']}, diameter {U['diameter']}",
        "Head loss": m["options"].get("HEADLOSS", "H-W (default)"),
        "Demand model": m["options"].get("DEMAND MODEL", "DDA (default)") + (
            f" (min {m['options'].get('MINIMUM PRESSURE', '0')}, required {m['options'].get('REQUIRED PRESSURE', '0.1')} {U['pressure']})"
            if m["options"].get("DEMAND MODEL", "").upper() == "PDA" else ""),
        "Simulation": ("steady state (single period)" if not dur else
                       f"extended period {hhmm(dur)} h, hydraulic step {hhmm(hyd)}, report step {hhmm(rep)}"),
        "Water quality": qual,
        "Nodes": ", ".join(f"{len(m[k])} {k}" for k in ("junctions", "reservoirs", "tanks") if m[k]) or None,
        "Links": ", ".join(f"{len(m[k])} {k}" for k in ("pipes", "pumps", "valves") if m[k]) or None,
        "Pipes": (f"total {g(length)} {U['length']}; diameter {rng([num(r['Diameter']) for r in pipes], ' ' + U['diameter'])}; "
                  f"roughness {rng([num(r['Roughness']) for r in pipes])}") if pipes else None,
        "Base demand": f"{g(consume)} {U['flow']}" + (f" x multiplier {g(mult)}" if mult != 1 else "")
                       + (f"; plus {g(inflow)} {U['flow']} injected at nodes with negative demand" if inflow else ""),
        "Demand range over the day": peak,
        "Elevations (junctions)": rng([num(r["Elevation"]) for r in m["junctions"]], " " + U["length"]),
        "Patterns / curves": f"{len(m['patterns'])} / {len(m['curves'])}",
        "Controls / rules": f"{len(m['controls'])} / {len(m['rules'])}" if m["controls"] or m["rules"] else None,
        "Emitters / sources": f"{len(m['emitters'])} / {len(m['sources'])}" if m["emitters"] or m["sources"] else None,
    }
    W = checks(m, nodes, links, dem)
    sections = [{"section": s, "rows": len(inp["raw"][s]) if s in RAW_SECTIONS else len(inp["sections"].get(s, []))}
                for s in inp["order"]]
    return {"name": path.stem, "path": str(path), "folder": str(path.parent), "model": m, "nodes": nodes, "links": links,
            "demands": dem, "pumps": [pump_curve(m, p) for p in m["pumps"]], "tanks": tanks,
            "pipe_diameters": sorted(((d, l) for d, l in by_d.items() if d is not None)),
            "summary": {k: v for k, v in S.items() if v not in (None, "")}, "warnings": W, "sections": sections, "_inp": inp}

# --------------------------------------------------------------------------- #
# Binary output file (.out) - OWA-EPANET "Output File" layout
# --------------------------------------------------------------------------- #
MAGIC = 516114521


class EpaOut:
    def __init__(self, path: Path):
        self.path = path
        f = self.f = open(path, "rb")
        f.seek(0, 2)
        size = f.tell()
        if size < 912:
            raise ValueError("file too small for an EPANET output file")
        f.seek(0)
        h = struct.unpack("<15i", f.read(60))
        (magic, self.version, self.nnodes, self.ntanks, self.nlinks, self.npumps, self.nvalves, self.qual, self.trace,
         self.flow_units, self.press_units, self.statistic, self.rep_start, self.rep_step, self.duration) = h
        if magic != MAGIC:
            raise ValueError("not an EPANET binary output file (bad magic number)")
        txt = lambda n: f.read(n).split(b"\0", 1)[0].decode("latin-1").strip()
        self.title = [txt(80) for _ in range(3)]
        self.inp_file, self.rpt_file = txt(260), txt(260)
        self.chem, self.chem_units = txt(32), txt(32)
        self.node_ids = [txt(32) for _ in range(self.nnodes)]
        self.link_ids = [txt(32) for _ in range(self.nlinks)]
        rd_i = lambda n: list(struct.unpack(f"<{n}i", f.read(4 * n))) if n else []
        rd_f = lambda n: list(struct.unpack(f"<{n}f", f.read(4 * n))) if n else []
        self.link_n1, self.link_n2, self.link_type = rd_i(self.nlinks), rd_i(self.nlinks), rd_i(self.nlinks)
        self.tank_index, self.tank_area = rd_i(self.ntanks), rd_f(self.ntanks)
        self.node_elev = rd_f(self.nnodes)
        self.link_len, self.link_diam = rd_f(self.nlinks), rd_f(self.nlinks)
        self.energy = []
        for _ in range(self.npumps):
            idx = struct.unpack("<i", f.read(4))[0]
            self.energy.append((idx,) + struct.unpack("<6f", f.read(24)))
        self.demand_charge = struct.unpack("<f", f.read(4))[0]
        self.data_pos = f.tell()
        self.nvals = 4 * self.nnodes + 8 * self.nlinks
        f.seek(size - 28)
        e = struct.unpack("<4f3i", f.read(28))
        self.react_bulk, self.react_wall, self.react_tank, self.source_in = e[:4]
        self.nperiods, self.warning_flag, magic2 = e[4:]
        self.complete = magic2 == MAGIC
        if not self.complete:
            self.nperiods = max(0, (size - self.data_pos) // (4 * self.nvals))

    def close(self):
        self.f.close()

    def periods(self):
        self.f.seek(self.data_pos)
        for _ in range(self.nperiods):
            buf = self.f.read(4 * self.nvals)
            if len(buf) < 4 * self.nvals:
                return
            a = array.array("f")
            a.frombytes(buf)
            if sys.byteorder != "little":
                a.byteswap()
            yield a


def read_out(path: Path, max_series=3000) -> dict:
    o = EpaOut(path)
    try:
        N, L = o.nnodes, o.nlinks
        tank_res = {o.tank_index[i] - 1: o.tank_area[i] for i in range(o.ntanks)}     # area 0 = reservoir
        res = {"file": path.name, "size_MB": round(path.stat().st_size / 1e6, 2), "format version": o.version,
               "flow units": FLOW_CODES[o.flow_units] if 0 <= o.flow_units < len(FLOW_CODES) else o.flow_units,
               "pressure units": PRESS_CODES[o.press_units] if 0 <= o.press_units < len(PRESS_CODES) else o.press_units,
               "nodes / links": f"{N} / {L} ({o.ntanks - sum(1 for a in o.tank_area if a > 0)} reservoirs, "
                                f"{sum(1 for a in o.tank_area if a > 0)} tanks, {o.npumps} pumps, {o.nvalves} valves)",
               "water quality": QUAL_CODES[o.qual] if 0 <= o.qual < 4 else o.qual,
               "report start / step": f"{hhmm(o.rep_start)} / {hhmm(o.rep_step)}", "duration": hhmm(o.duration),
               "periods": o.nperiods, "statistic": STAT_CODES[o.statistic] if 0 <= o.statistic < 5 else o.statistic,
               "warnings during run": "yes" if o.warning_flag else "no",
               "complete": "yes" if o.complete else "NO - epilog missing (run aborted?)"}
        if o.qual == 1:
            res["chemical"] = f"{o.chem} ({o.chem_units})"
        if o.title and o.title[0]:
            res["title"] = o.title[0]
        junc = [i for i in range(N) if i not in tank_res]
        pmin = [(1e30, 0)] * N
        pmax = [(-1e30, 0)] * N
        dmax = [0.0] * N
        qmin, qmax = [1e30] * N, [-1e30] * N
        lq = [(0.0, 0)] * L
        lv = [(0.0, 0)] * L
        lhl = [0.0] * L
        lclosed = [0] * L
        pump_pts = {i: [] for i in range(L) if o.link_type[i] == 2}
        times, sysdem, minp, tanklev, resflow = [], [], [], {i: [] for i in tank_res}, {}
        stride = max(1, o.nperiods // max_series)
        tank_elev = {i: o.node_elev[i] for i in tank_res}
        for k, a in enumerate(o.periods()):
            t = o.rep_start + k * o.rep_step if o.statistic == 0 else 0
            D, H, P, Q = 0, N, 2 * N, 3 * N
            F, V, HL, ST = 4 * N, 4 * N + L, 4 * N + 2 * L, 4 * N + 4 * L
            sd, mp = 0.0, 1e30
            for i in junc:
                p = a[P + i]
                if p < pmin[i][0]:
                    pmin[i] = (p, t)
                if p > pmax[i][0]:
                    pmax[i] = (p, t)
                d = a[D + i]
                if d > dmax[i]:
                    dmax[i] = d
                sd += d
                mp = min(mp, p)
            for i in range(N):
                qq = a[Q + i]
                qmin[i] = min(qmin[i], qq)
                qmax[i] = max(qmax[i], qq)
            for j in range(L):
                q, v = abs(a[F + j]), abs(a[V + j])
                if q > lq[j][0]:
                    lq[j] = (q, t)
                if v > lv[j][0]:
                    lv[j] = (v, t)
                hl = abs(a[HL + j])
                if hl > lhl[j]:
                    lhl[j] = hl
                if a[ST + j] <= 2:
                    lclosed[j] += 1
                if j in pump_pts and a[ST + j] > 2 and k % stride == 0:
                    pump_pts[j].append((a[F + j], abs(a[HL + j])))
            if k % stride == 0:
                times.append(t / 3600)
                sysdem.append(sd)
                minp.append(mp if junc else None)
                for i in tank_res:
                    if tank_res[i] > 0:
                        tanklev[i].append(a[H + i] - tank_elev[i])
                    else:
                        resflow.setdefault(i, []).append(-a[D + i])
        nid, lid = o.node_ids, o.link_ids
        res["_times_h"] = times
        res["_system_demand"] = sysdem
        res["_min_pressure"] = minp
        res["_tank_level"] = {nid[i]: v for i, v in tanklev.items() if tank_res[i] > 0}
        res["_reservoir_out"] = {nid[i]: v for i, v in resflow.items()}
        res["_node"] = {nid[i]: {"pmin": pmin[i][0], "pmin_t": pmin[i][1], "pmax": pmax[i][0], "dmax": dmax[i],
                                 "qmin": qmin[i], "qmax": qmax[i]} for i in junc}
        res["_link"] = {lid[j]: {"type": LINK_CODES[o.link_type[j]] if 0 <= o.link_type[j] < len(LINK_CODES) else o.link_type[j],
                                 "qmax": lq[j][0], "qmax_t": lq[j][1], "vmax": lv[j][0], "vmax_t": lv[j][1], "hlmax": lhl[j],
                                 "closed_periods": lclosed[j]} for j in range(L)}
        res["_pump_points"] = {lid[j]: v for j, v in pump_pts.items()}
        res["_energy"] = [{"pump": lid[e[0] - 1] if 0 < e[0] <= L else e[0], "utilization %": e[1], "efficiency %": e[2],
                           "kWh per unit volume": e[3], "average kW": e[4], "peak kW": e[5], "cost per day": e[6]}
                          for e in o.energy]
        if res["_energy"]:
            res["demand charge"] = g(o.demand_charge)
        if o.qual:
            res["reaction rates bulk / wall / tank / source (mass/h)"] = " / ".join(g(x) for x in
                                                                                  (o.react_bulk, o.react_wall, o.react_tank, o.source_in))
        neg = sorted(((n, v["pmin"], v["pmin_t"]) for n, v in res["_node"].items() if v["pmin"] < 0), key=lambda x: x[1])
        res["junctions with negative pressure"] = len(neg)
        res["_negative"] = neg
        if res["_node"]:
            lo = min(res["_node"].items(), key=lambda kv: kv[1]["pmin"])
            hi = max(res["_node"].items(), key=lambda kv: kv[1]["pmax"])
            res["junction pressure range"] = f"{g(lo[1]['pmin'])} ({lo[0]} at {hhmm(lo[1]['pmin_t'])}) .. {g(hi[1]['pmax'])} ({hi[0]})"
        pv = [(n, v) for n, v in res["_link"].items() if v["type"] in ("pipe", "CV pipe")]
        if pv:
            vm = max(pv, key=lambda kv: kv[1]["vmax"])
            res["max pipe velocity"] = f"{g(vm[1]['vmax'])} in {vm[0]} at {hhmm(vm[1]['vmax_t'])}"
        if o.qual and res["_node"]:
            res["junction quality range"] = f"{g(min(v['qmin'] for v in res['_node'].values()))} .. " \
                                            f"{g(max(v['qmax'] for v in res['_node'].values()))}" + (" h" if o.qual == 2 else "")
        return res
    finally:
        o.close()


def _safe_out(p):
    try:
        return read_out(p)
    except Exception as e:
        return {"file": p.name, "error": f"{type(e).__name__}: {e}"}


# --------------------------------------------------------------------------- #
# Report file (.rpt)
# --------------------------------------------------------------------------- #
def parse_rpt(path: Path) -> dict:
    txt = path.read_text(errors="replace")
    lines = txt.splitlines()
    out = {"file": path.name}
    m = re.search(r"Version\s+([\d.]+)", txt)
    out["engine"] = f"EPANET {m.group(1)}" if m else None
    for k, pat in (("analysis begun", r"Analysis begun\s+(.+)"), ("analysis ended", r"Analysis ended\s+(.+)")):
        mm = re.search(pat, txt)
        if mm:
            out[k] = mm.group(1).strip()
    out["errors"] = list(dict.fromkeys(ln.strip() for ln in lines if re.match(r"\s*Error\s+\d+", ln)))
    out["warnings"] = [ln.strip() for ln in lines if re.match(r"\s*WARNING", ln)]
    status = [ln.strip() for ln in lines if re.match(r"\s*\d+:\d\d:\d\d:", ln)]
    out["status messages"] = len(status)
    out["pump/valve changes"] = [s for s in status if re.search(r"changed", s)]
    out["unbalanced"] = [s for s in status if re.search(r"unbalanced|not converge", s, re.I)]
    trials = [int(x) for x in re.findall(r"Balanced after (\d+) trials", txt)]
    if trials:
        out["trials per step"] = f"max {max(trials)}, mean {sum(trials) / len(trials):.1f} over {len(trials)} solutions"
    # Energy usage table
    i = next((k for k, ln in enumerate(lines) if "Energy Usage" in ln), None)
    if i is not None:
        rows = []
        for ln in lines[i + 1:i + 40]:
            t = ln.split()
            if len(t) >= 7 and all(num(x) is not None for x in t[1:7]):
                rows.append(t)
        out["energy"] = rows
    return out


# --------------------------------------------------------------------------- #
# Build
# --------------------------------------------------------------------------- #
def find_inputs(arg: str) -> list[Path]:
    p = Path(arg).expanduser().resolve()
    if not p.exists():
        sys.exit(f"Path not found: {p}")
    if p.is_file():
        return [p]
    fs = sorted(x for x in p.rglob("*") if x.suffix.lower() in (".inp", ".net"))
    inp_stems = {f.with_suffix("").as_posix().lower() for f in fs if f.suffix.lower() == ".inp"}
    # a .net project with an exported .inp of the same name adds nothing - show the .inp only
    fs = [f for f in fs if f.suffix.lower() == ".inp" or f.with_suffix("").as_posix().lower() not in inp_stems]
    if not fs:
        sys.exit(f"No .inp file found under {p}")
    return sorted(set(fs))


def build_one(path: Path, rpt=None, out=None) -> dict:
    base = {"name": path.stem, "path": str(path), "folder": str(path.parent), "summary": {}, "warnings": [], "sections": []}
    if path.suffix.lower() == ".out":
        base["only_out"] = True
        base["out"] = _safe_out(path)
        return base
    try:
        inp = read_inp(path)
    except ValueError as e:
        base["error"] = str(e)
        if path.with_suffix(".inp").exists():
            base["error"] += f" - an .inp with the same name exists ({path.with_suffix('.inp').name})"
        return base
    cfg = summarize(path, inp)
    rp = Path(rpt) if rpt else path.with_suffix(".rpt")
    op = Path(out) if out else path.with_suffix(".out")
    if not op.exists() and path.with_suffix(".bin").exists():
        op = path.with_suffix(".bin")
    if rp.exists():
        cfg["rpt"] = parse_rpt(rp)
    if op.exists():
        cfg["out"] = _safe_out(op)
    W = cfg["warnings"]
    r = cfg.get("rpt")
    if r:
        for e in r["errors"]:
            W.append(f"{r['file']}: {e}")
        if r["warnings"]:
            W.append(f"{r['file']}: {len(r['warnings'])} warning(s) - {r['warnings'][0]}")
    o = cfg.get("out")
    if o and o.get("error"):
        W.append(f"{o['file']}: {o['error']}")
    elif o:
        if o.get("_negative"):
            n = o["_negative"]
            W.append(f"{len(n)} junction(s) reach negative pressure, lowest {g(n[0][1])} {cfg['model']['units']['pressure']} "
                     f"at {n[0][0]} ({hhmm(n[0][2])}): {', '.join(x[0] for x in n[:8])}")
        if not str(o.get("complete")).startswith("yes"):
            W.append(f"{o['file']}: {o['complete']}")
        U = cfg["model"]["units"]
        vlim = 3.0 if U["system"] == "SI" else 10.0
        fast = [(n, v["vmax"]) for n, v in o["_link"].items() if v["type"] in ("pipe", "CV pipe") and v["vmax"] > vlim]
        if fast:
            W.append(f"{len(fast)} pipe(s) above {g(vlim)} {U['velocity']}: " + ", ".join(f"{n} {g(v)}" for n, v in fast[:6]))
        try:
            if path.stat().st_mtime > op.stat().st_mtime + 1:
                W.append(f"{path.name} is newer than {op.name} - results may be out of date")
        except OSError:
            pass
    return cfg


def build_config(arg, rpt=None, out=None) -> dict:
    files = find_inputs(arg)
    root = Path(arg).expanduser().resolve()
    return {"root": str(root if root.is_dir() else root.parent), "name": root.stem if root.is_file() else root.name,
            "models": [build_one(f, rpt, out) for f in files]}


# --------------------------------------------------------------------------- #
# Terminal output
# --------------------------------------------------------------------------- #


def text_table(headers, rows, indent=2, maxw=18, limit=None):
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
    m = cfg["model"]
    if sec == "PUMPS":
        return ["ID", "Node1", "Node2", "Parameters"], [[p["ID"], p["Node1"], p["Node2"],
                                                          " ".join(f"{k} {v}" for k, v in p["params"].items())] for p in m["pumps"]]
    if sec == "DEMANDS":
        return COLS["DEMANDS"], [[j, g(d["demand"]), d["pattern"], d["category"]] for j, ds in m["demands"].items() for d in ds]
    rows = m.get(sec.lower(), [])
    cols = COLS.get(sec, [])
    return cols, [[r.get(c, "") for c in cols] for r in rows]


def show_model(cfg, full):
    h1(f"EPANET MODEL: {cfg['name']}")
    row("File", cfg["path"])
    if cfg.get("error"):
        print(f"  {cfg['error']}")
        return
    if cfg.get("only_out"):
        show_out(cfg)
        return
    m, U = cfg["model"], cfg["model"]["units"]
    for t in m["title"]:
        print(f"  | {t}")
    h2("Summary")
    for k, v in cfg["summary"].items():
        row(k, v)
    h2("Options")
    for k in OPTION_KEYS:
        if k in m["options"] or k in ("UNITS", "HEADLOSS", "DEMAND MODEL", "TRIALS", "ACCURACY", "UNBALANCED", "PATTERN"):
            v = m["options"].get(k)
            row(k, (v if v is not None else f"{OPTION_DEFAULTS.get(k, '-')} (default)") + (f"   ({OPTION_HELP[k]})" if OPTION_HELP.get(k) else ""))
    for k, v in m["options"].items():
        if k not in OPTION_KEYS:
            row(k, v + "   (unrecognised keyword)")
    h2("Times")
    for k in TIME_KEYS:
        v = m["times"].get(k)
        if v is not None or k in TIME_DEFAULTS:
            row(k, v if v is not None else f"{TIME_DEFAULTS.get(k, '1/10 of hydraulic step')} (default)")
    h2("Sections in the file")
    print("  " + ", ".join(f"{s['section']} ({s['rows']})" for s in cfg["sections"]))
    lim = None if full else 15
    for sec in ("JUNCTIONS", "RESERVOIRS", "TANKS", "PIPES", "PUMPS", "VALVES", "DEMANDS", "EMITTERS", "LEAKAGE", "STATUS",
                "QUALITY", "SOURCES", "MIXING"):
        cols, rows = section_rows(cfg, sec)
        if rows:
            h2(f"[{sec}] ({len(rows)})")
            text_table(cols, rows, limit=lim)
    if cfg["tanks"]:
        h2("Tanks (computed)")
        for t in cfg["tanks"]:
            row(t["id"], f"base {g(t['base'])} {U['length']}, levels min/init/max {t['levels']}, diameter {g(t['diameter'])}"
                + (f", working volume {g(t['working volume'])} {U['volume']}" if t["working volume"] else f", volume curve {t['curve']}"))
    if cfg["pumps"]:
        h2("Pumps")
        for p in cfg["pumps"]:
            row(p["id"], {k: v for k, v in p.items() if not k.startswith("_") and k != "id"})
    if cfg["pipe_diameters"]:
        h2("Pipe length by diameter")
        text_table([f"Diameter ({U['diameter']})", f"Length ({U['length']})"], [[g(d), g(l)] for d, l in cfg["pipe_diameters"]])
    if m["patterns"]:
        h2("Patterns")
        for n, f in m["patterns"].items():
            print(f"  {n:<16} {len(f)} factors, {rng(f)}, mean {g(sum(f) / len(f)) if f else '-'}")
    if m["curves"]:
        h2("Curves")
        for n, c in m["curves"].items():
            print(f"  {n:<16} {c['type'] or '?':<10} {len(c['points'])} point(s): " + ", ".join(f"({g(x)}, {g(y)})" for x, y in c["points"][:6]))
    if m["controls"]:
        h2(f"Simple controls ({len(m['controls'])})")
        for c in m["controls"]:
            print(f"  {c}")
    if m["rules"]:
        h2(f"Rule-based controls ({len(m['rules'])})")
        if full:
            for ln in m["rules_text"]:
                print("  " + ln)
        else:
            print("  " + ", ".join(m["rules"]))
    for lab, key in (("Energy", "energy"), ("Reactions", "reactions")):
        if m[key]:
            h2(lab)
            for ln in m[key]:
                print("  " + ln)
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
    for k in ("file", "engine", "analysis begun", "analysis ended", "trials per step", "status messages"):
        row(k, r.get(k))
    for x in r["errors"]:
        print(f"    {x}")
    for x in r["warnings"][:15]:
        print(f"    {x}")
    for x in r["pump/valve changes"][:12]:
        print(f"    {x}")


def show_out(cfg):
    o = cfg.get("out")
    h2("Binary output (.out)")
    if not o:
        print("  none found")
        return
    if o.get("error"):
        print(f"  {o['file']}: {o['error']}")
        return
    for k, v in o.items():
        if not k.startswith("_"):
            row(k, v)
    if o.get("_negative"):
        print("    Junctions with negative pressure (lowest first)")
        for n, p, t in o["_negative"][:12]:
            print(f"      {n:<16} {g(p):>10} at {hhmm(t)}")
    for n, v in o["_tank_level"].items():
        if v:
            print(f"    tank {n}: level {g(min(v))} .. {g(max(v))}")
    for e in o["_energy"]:
        print(f"    pump {e['pump']}: utilization {g(e['utilization %'])} %, efficiency {g(e['efficiency %'])} %, "
              f"avg {g(e['average kW'])} kW, peak {g(e['peak kW'])} kW, cost/day {g(e['cost per day'])}")


def show(cfg, full=False):
    for m in cfg["models"]:
        show_model(m, full)
    print()


# --------------------------------------------------------------------------- #
# HTML template & helpers (same layout as ras_config.py / hms_config.py)
# --------------------------------------------------------------------------- #
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<script>(function(){var t;try{t=localStorage.getItem('modelvue-theme')}catch(e){}if(t!=='dark')t='light';document.documentElement.dataset.theme=t;})();</script>
<title>EPANET · @@TITLE@@</title>
<style>
:root{--bg:#f6f7f9;--panel:#fff;--ink:#1d2433;--muted:#667085;--line:#e4e7ec;--acc:#0d6e8a;--accs:#e1f1f6;--warn:#b42318}
:root[data-theme=dark]{color-scheme:dark;--bg:#0f1318;--panel:#171c23;--ink:#e6e9ef;--muted:#98a2b3;--line:#2a313b;--acc:#7fcbe0;--accs:#14303a;--warn:#f97066}
*{box-sizing:border-box} body{margin:0;font:14px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,"Noto Sans TC",sans-serif;background:var(--bg);color:var(--ink)}
header{position:sticky;top:0;z-index:5;background:var(--acc);color:#fff;padding:14px 24px;display:flex;gap:16px;align-items:baseline;flex-wrap:wrap}
header h1{margin:0;font-size:20px} header span{opacity:.85}
[data-theme=dark] header{background:#0b4152}
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
<header><h1>EPANET · @@TITLE@@</h1><span>@@SUB@@</span><button class="theme" type="button" title="Switch between the light and dark theme">☾ Dark</button></header>
<div class="wrap"><nav><div class="tools"><button data-x="1">Expand all</button><button data-x="0">Collapse all</button></div>@@NAV@@</nav>
<main>@@BODY@@
<footer>Generated by modelvue.epanetvue on @@GEN@@ from <code>@@PATH@@</code></footer></main></div>
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
# SVG graphics
# --------------------------------------------------------------------------- #
PALETTE = ["#1f5fa8", "#d9622b", "#2e7d32", "#8e5bb5", "#c62828", "#b8860b", "#16a2a2", "#7f8c8d"]
LINK_STYLE = {"pump": "#c2185b", "PRV": "#00897b", "PSV": "#00897b", "PBV": "#00897b", "FCV": "#00897b", "TCV": "#00897b",
              "GPV": "#00897b", "PCV": "#00897b"}


def _frame(w, h, L, R, T, B, x0, x1, y0, y1, xlabel, ylabel):
    X = lambda v: L + (v - x0) / ((x1 - x0) or 1) * (w - L - R)
    Y = lambda v: T + (y1 - v) / ((y1 - y0) or 1) * (h - T - B)
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


def svg_lines(series, xlabel, ylabel, w=820, h=230, hlines=(), ymin=None) -> str:
    xs = [x for _, _, pts in series for x, _ in pts]
    ys = [y for _, _, pts in series for _, y in pts if y is not None] + [v for v, _, _ in hlines]
    if not xs or not ys:
        return ""
    y0, y1 = min(ys), max(ys)
    if ymin is not None:
        y0 = min(y0, ymin)
    pad = (y1 - y0) * 0.08 or 1
    X, Y, o = _frame(w, h, 58, 12, 14, 34, min(xs), max(xs), y0 - pad, y1 + pad, xlabel, ylabel)
    for v, lab, col in hlines:
        o.append(f'<line x1="58" x2="{w - 12}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="{col}" stroke-dasharray="5 4"/>'
                 f'<text x="{w - 16}" y="{Y(v) - 4:.1f}" text-anchor="end" fill="{col}">{_h.escape(lab)}</text>')
    lx = 68
    for i, (lab, col, pts) in enumerate(series):
        col = col or PALETTE[i % len(PALETTE)]
        P = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in pts if b is not None)
        o.append(f'<polyline points="{P}" fill="none" stroke="{col}" stroke-width="1.8"><title>{_h.escape(lab)}</title></polyline>')
        if len(series) > 1:
            o.append(f'<rect x="{lx}" y="20" width="12" height="3" fill="{col}"/><text x="{lx + 16}" y="25">{_h.escape(lab)}</text>')
            lx += 30 + 6.5 * len(lab)
    return f'<svg class="plot" viewBox="0 0 {w} {h}">' + "".join(o) + "</svg>"


def svg_pattern(name, f, step_h, w=820, h=190) -> str:
    if not f:
        return ""
    n = len(f)
    X, Y, o = _frame(w, h, 52, 12, 10, 34, 0, n * step_h, 0, max(f) * 1.1, "hours", "multiplier")
    mean = sum(f) / n
    for i, v in enumerate(f):
        o.append(f'<rect x="{X(i * step_h) + 0.5:.1f}" y="{Y(v):.1f}" width="{max(1, X((i + 1) * step_h) - X(i * step_h) - 1):.1f}" '
                 f'height="{Y(0) - Y(v):.1f}" class="bar"><title>{g(i * step_h)} h: {g(v)}</title></rect>')
    o.append(f'<line x1="52" x2="{w - 12}" y1="{Y(mean):.1f}" y2="{Y(mean):.1f}" stroke="#c62828" stroke-dasharray="5 4"/>'
             f'<text x="{w - 16}" y="{Y(mean) - 4:.1f}" text-anchor="end" fill="#c62828">mean {g(mean)}</text>')
    return f'<h4>Pattern {_h.escape(name)}</h4><svg class="plot" viewBox="0 0 {w} {h}">' + "".join(o) + "</svg>"


def svg_pump(p, U, pts_run, w=640, h=260) -> str:
    cur = p.get("_curve")
    if not cur:
        return ""
    xs = [x for x, _ in cur] + [x for x, _ in pts_run]
    ys = [y for _, y in cur] + [y for _, y in pts_run]
    X, Y, o = _frame(w, h, 58, 12, 14, 34, 0, max(xs) * 1.05, 0, max(ys) * 1.1, f"flow ({U['flow']})", f"head gain ({U['head']})")
    o.append(f'<polyline points="{" ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in cur)}" fill="none" stroke="#c2185b" stroke-width="2"/>')
    for x, y in p.get("_points", []):
        o.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="5" fill="#fff" stroke="#c2185b" stroke-width="2"><title>curve point {g(x)}, {g(y)}</title></circle>')
    for x, y in pts_run:
        o.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="3" fill="#1f5fa8" opacity=".7"/>')
    leg = '<div class="legend"><span><i class="ln" style="border-color:#c2185b"></i>pump curve</span>' + (
        '<span><i class="mk" style="background:#1f5fa8;border-radius:50%"></i>operating points from .out</span>' if pts_run else "") + "</div>"
    return f'<h4>Pump {_h.escape(p["id"])} - {_h.escape(p.get("kind", ""))}</h4><svg class="plot" viewBox="0 0 {w} {h}">' + "".join(o) + "</svg>" + leg


def _pcolor(p, lo, hi):
    if p < 0:
        return "#c62828"
    t = max(0.0, min(1.0, (p - lo) / ((hi - lo) or 1)))
    stops = [(0, (253, 174, 97)), (0.5, (171, 217, 233)), (1, (44, 123, 182))]
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t <= t1:
            f = (t - t0) / (t1 - t0)
            return "#%02x%02x%02x" % tuple(round(a + (b - a) * f) for a, b in zip(c0, c1))
    return "#2c7bb6"


def svg_map(cfg, results=False, w=900) -> str:
    m, nodes, links = cfg["model"], cfg["nodes"], cfg["links"]
    C = m["coords"]
    pts = [p for p in list(C.values()) + [q for v in m["vertices"].values() for q in v] if p[0] is not None and p[1] is not None]
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
    out = cfg.get("out") if results else None
    lr, nr = (out or {}).get("_link", {}), (out or {}).get("_node", {})
    many = len(nodes) > 120
    ds = [l.get("diameter") for l in links.values() if l.get("diameter")]
    dmax = max(ds) if ds else 1
    vmax = max((v["vmax"] for v in lr.values()), default=1) or 1
    o = [f'<svg class="map" viewBox="0 0 {w:.0f} {h:.0f}">']
    for lk, l in links.items():
        a, b = C.get(l["n1"]), C.get(l["n2"])
        if not a or not b or None in a or None in b:
            continue
        P = [a] + m["vertices"].get(lk, []) + [b]
        wd = 1 + 3 * (l.get("diameter") or dmax * 0.3) / dmax if l["type"] in ("pipe", "CV pipe") else 2
        col = LINK_STYLE.get(l["type"], "#6b7a90")
        tip = f"{l['type']} {lk}: {l['n1']} -> {l['n2']}" + (f", D {g(l.get('diameter'))}" if l.get("diameter") else "")
        if results and lk in lr:
            v = lr[lk]["vmax"]
            if l["type"] in ("pipe", "CV pipe"):
                t = v / vmax
                col = "#%02x%02x%02x" % (round(44 + (215 - 44) * t), round(123 + (48 - 123) * t), round(182 + (39 - 182) * t))
            tip += f"; max velocity {g(v)}, max |flow| {g(lr[lk]['qmax'])}"
        o.append(f'<polyline points="{" ".join(f"{X(p[0]):.1f},{Y(p[1]):.1f}" for p in P)}" fill="none" stroke="{col}" '
                 f'stroke-width="{wd:.1f}" stroke-linecap="round"><title>{_h.escape(tip)}</title></polyline>')
        if l["type"] not in ("pipe", "CV pipe"):
            mx, my = (X(a[0]) + X(b[0])) / 2, (Y(a[1]) + Y(b[1])) / 2
            if l["type"] == "pump":
                o.append(f'<circle cx="{mx:.1f}" cy="{my:.1f}" r="7" fill="#fff" stroke="{col}" stroke-width="2"><title>{_h.escape(tip)}</title></circle>'
                         f'<polygon points="{mx - 3:.1f},{my - 4:.1f} {mx + 4:.1f},{my:.1f} {mx - 3:.1f},{my + 4:.1f}" fill="{col}"/>')
            else:
                o.append(f'<polygon points="{mx - 7:.1f},{my - 5:.1f} {mx + 7:.1f},{my + 5:.1f} {mx + 7:.1f},{my - 5:.1f} {mx - 7:.1f},{my + 5:.1f}" '
                         f'fill="{col}"><title>{_h.escape(tip)}</title></polygon>')
    pv = [v["pmin"] for v in nr.values()]
    plo, phi = (min(pv), max(v["pmax"] for v in nr.values())) if pv else (0, 1)
    for n, d in nodes.items():
        p = C.get(n)
        if not p or None in p:
            continue
        x, y = X(p[0]), Y(p[1])
        tip = f"{d['type']} {n}: elevation/head {g(d['elev'])}"
        if d["type"] == "junction":
            col, r = "#1f5fa8", 3 if many else 5
            if results and n in nr:
                v = nr[n]
                col = _pcolor(v["pmin"], max(0, plo), phi)
                tip += f"; pressure {g(v['pmin'])} .. {g(v['pmax'])}"
                if v["pmin"] < 0:
                    o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r + 5}" fill="none" stroke="#c62828" stroke-width="2.5"/>')
            o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{col}" stroke="#fff"><title>{_h.escape(tip)}</title></circle>')
        elif d["type"] == "reservoir":
            o.append(f'<path d="M{x - 9:.1f},{y - 6:.1f} L{x - 9:.1f},{y + 6:.1f} L{x + 9:.1f},{y + 6:.1f} L{x + 9:.1f},{y - 6:.1f}" '
                     f'fill="#90caf9" stroke="#1565c0" stroke-width="2"><title>{_h.escape(tip)}</title></path>')
        else:
            o.append(f'<rect x="{x - 7:.1f}" y="{y - 10:.1f}" width="14" height="20" fill="#b3e5fc" stroke="#0277bd" stroke-width="2">'
                     f'<title>{_h.escape(tip)}</title></rect>')
        if not many:
            o.append(f'<text x="{x + 9:.1f}" y="{y - 8:.1f}" class="nlab">{_h.escape(n)}</text>')
    for x, y, t in m["labels"]:
        anchor = "end" if X(x) > 0.8 * w else "start"
        o.append(f'<text x="{X(x):.1f}" y="{Y(y):.1f}" class="lbl" text-anchor="{anchor}">{_h.escape(t)}</text>')
    o.append("</svg>")
    leg = ['<div class="legend"><span><b style="color:#1f5fa8">●</b> junction</span>',
           '<span><b style="color:#1565c0">⊔</b> reservoir</span><span><b style="color:#0277bd">▮</b> tank</span>',
           '<span><b style="color:#c2185b">▶</b> pump</span><span><b style="color:#00897b">⧗</b> valve</span>',
           '<span>pipe width = diameter</span>']
    if results and nr:
        leg.append('<span><i class="sw" style="background:linear-gradient(90deg,#fdae61,#abd9e9,#2c7bb6)"></i>minimum pressure, low .. high</span>'
                   '<span><i class="mk" style="border:2px solid #c62828;border-radius:50%"></i>negative pressure</span>'
                   '<span><i class="sw" style="background:linear-gradient(90deg,#2c7bb6,#d73027)"></i>pipe max velocity</span>')
    leg.append("</div>")
    return "".join(o) + "".join(leg)


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #
EXTRA_CSS = """
svg.map{width:100%;max-width:1000px;height:auto;display:block;margin:8px 0;background:var(--bg);border:1px solid var(--line);border-radius:6px}
svg.map .nlab{font-size:10px;fill:var(--ink)} svg.map .lbl{font-size:12px;font-style:italic;fill:var(--muted)}
svg.plot{width:100%;max-width:980px;height:auto;display:block;margin:8px 0;font-size:11px}
svg.plot text{fill:var(--muted)} svg.plot .lab{fill:var(--ink)} svg.plot .gl{stroke:var(--line)} svg.plot .frame{fill:none;stroke:var(--line)}
svg.plot .bar{fill:#5c9bd5}
.legend{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:12px;color:var(--muted);margin-bottom:8px}
.legend i{display:inline-block;vertical-align:middle;margin-right:5px}
.legend .mk{width:10px;height:10px} .legend .ln{width:22px;height:0;border-top:3px solid} .legend .sw{width:110px;height:10px;border-radius:2px}
td.msg{white-space:normal}
.desc{white-space:pre-wrap;color:var(--muted);font-size:13px;max-height:260px;overflow:auto;border:1px solid var(--line);border-radius:6px;padding:8px 12px;margin:6px 0}
"""


def _rows_html(cols, rows, limit=500):
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

    if cfg.get("error"):
        out.append(f'<p class="warn">{_e(cfg["error"])}</p>')
        return "".join(out), present
    if cfg.get("only_out"):
        add("out", "Binary output", _out_html(cfg))
        return "".join(out), present
    m, U = cfg["model"], cfg["model"]["units"]
    if m["title"]:
        out.append(f'<div class="desc">{_e(chr(10).join(m["title"]))}</div>')
    add("sum", "Summary", _two(list(cfg["summary"].items())))
    add("map", "Network map", svg_map(cfg), f"{len(m['coords'])} nodes with coordinates")
    orows = []
    for k in OPTION_KEYS:
        if k in m["options"]:
            orows.append([f"<b>{_e(k)}</b>", _e(m["options"][k]), (_e(OPTION_HELP.get(k, "")), "msg")])
        elif k in OPTION_DEFAULTS:
            orows.append([_e(k), f'<span class="muted">{_e(OPTION_DEFAULTS[k])} (default)</span>', (_e(OPTION_HELP.get(k, "")), "msg")])
    trows = [[(f"<b>{_e(k)}</b>" if k in m["times"] else _e(k)),
              _e(m["times"][k]) if k in m["times"] else f'<span class="muted">{_e(TIME_DEFAULTS.get(k, "1/10 of hydraulic step"))} (default)</span>']
             for k in TIME_KEYS]
    add("opt", "Options & times", '<div class="grid2"><div>' + _table(["Option", "Value", "Meaning"], orows, sortable=False)
        + "</div><div>" + _table(["Time", "Value"], trows, sortable=False) + "</div></div>"
        + (("<h4>Energy</h4><pre class='desc'>" + _e("\n".join(m["energy"])) + "</pre>") if m["energy"] else "")
        + (("<h4>Reactions</h4><pre class='desc'>" + _e("\n".join(m["reactions"])) + "</pre>") if m["reactions"] else ""))
    for title, secs in (("Nodes", ("JUNCTIONS", "RESERVOIRS", "TANKS", "DEMANDS", "EMITTERS")),
                        ("Links", ("PIPES", "PUMPS", "VALVES", "STATUS", "LEAKAGE")),
                        ("Water quality", ("QUALITY", "SOURCES", "MIXING"))):
        body, n = "", 0
        for sec in secs:
            cols, rows = section_rows(cfg, sec)
            if rows:
                n += len(rows)
                body += _blk(f"[{sec}]", _rows_html(cols, rows), f"{len(rows)} rows", open_=sec == secs[0])
        if title == "Nodes" and cfg["tanks"]:
            body += "<h4>Tanks</h4>" + _rows_html(["Tank", "Base elev", "Levels min/init/max", "Diameter", f"Working volume ({U['volume']})", "Volume curve"],
                                                  [[t["id"], g(t["base"]), t["levels"], g(t["diameter"]), g(t["working volume"]), t["curve"]] for t in cfg["tanks"]])
        if title == "Links" and cfg["pipe_diameters"]:
            body += "<h4>Pipe length by diameter</h4>" + _rows_html([f"Diameter ({U['diameter']})", f"Length ({U['length']})"],
                                                                    [[g(d), g(l)] for d, l in cfg["pipe_diameters"]])
        add(title.lower().replace(" ", ""), title, body, f"{n} rows", open_=title != "Water quality")
    step_h = (parse_seconds(m["times"].get("PATTERN TIMESTEP", "1:00")) or 3600) / 3600
    body = "".join(svg_pattern(n, f, step_h) for n, f in list(m["patterns"].items())[:12])
    o = cfg.get("out") or {}
    for p in cfg["pumps"]:
        body += svg_pump(p, U, o.get("_pump_points", {}).get(p["id"], []))
    if m["curves"]:
        body += "<h4>Curves</h4>" + _rows_html(["Curve", "Type", "Points", "Comment"],
                                               [[n, c["type"], " ".join(f"({g(x)}, {g(y)})" for x, y in c["points"][:8]), c.get("comment", "")]
                                                for n, c in m["curves"].items()])
    add("pat", "Patterns, curves & pumps", body, open_=True)
    body = ""
    if m["controls"]:
        body += f"<h4>Simple controls ({len(m['controls'])})</h4><pre class='desc'>{_e(chr(10).join(m['controls']))}</pre>"
    if m["rules_text"]:
        body += f"<h4>Rule-based controls ({len(m['rules'])})</h4><pre class='desc'>{_e(chr(10).join(m['rules_text']))}</pre>"
    add("ctl", "Controls", body)
    add("secs", "Sections in the file", _table(["Section", "Rows"], [[_e(s["section"]), _e(s["rows"])] for s in cfg["sections"]]),
        f"{len(cfg['sections'])} sections", open_=False)
    add("rpt", "Report file (.rpt)", _rpt_html(cfg))
    add("out", "Binary output (.out)", _out_html(cfg))
    return "".join(out), present


def _rpt_html(cfg):
    r = cfg.get("rpt")
    if not r:
        return '<p class="muted">No report file found (model not run, or saved elsewhere - use --rpt).</p>'
    body = _two([(k, r.get(k)) for k in ("file", "engine", "analysis begun", "analysis ended", "trials per step", "status messages")])
    body += _table(["Level", "Message"], [['<span class="warn">ERROR</span>', (_e(x), "msg")] for x in r["errors"]]
                   + [["warning", (_e(x), "msg")] for x in r["warnings"]], sortable=False)
    if r["pump/valve changes"]:
        body += "<h4>Pump / valve status changes</h4><pre class='desc'>" + _e("\n".join(r["pump/valve changes"])) + "</pre>"
    if r.get("energy"):
        body += "<h4>Energy usage</h4>" + _rows_html(["Pump", "Usage %", "Avg effic %", "kWh/unit", "Avg kW", "Peak kW", "Cost/day"], r["energy"])
    return body


def _out_html(cfg):
    o = cfg.get("out")
    if not o:
        return '<p class="muted">No binary output file found.</p>'
    if o.get("error"):
        return f'<p class="warn">{_e(o["file"])}: {_e(o["error"])}</p>'
    U = cfg["model"]["units"] if cfg.get("model") else units_for(str(o.get("flow units", "GPM")))
    body = _two([(k, v) for k, v in o.items() if not k.startswith("_")])
    t = o["_times_h"]
    if len(t) > 1:
        body += "<h4>System demand and lowest junction pressure</h4>"
        body += svg_lines([("total junction demand", "#1f5fa8", list(zip(t, o["_system_demand"])))], "hours", f"demand ({o['flow units']})", h=200)
        body += svg_lines([("lowest junction pressure", "#c62828", list(zip(t, o["_min_pressure"])))], "hours",
                          f"pressure ({o['pressure units']})", h=200, hlines=[(0, "0", "#555")])
        if o["_tank_level"]:
            body += "<h4>Tank levels</h4>" + svg_lines([(n, None, list(zip(t, v))) for n, v in o["_tank_level"].items()], "hours",
                                                        f"level ({U['length']})", h=220)
        if o["_reservoir_out"]:
            body += "<h4>Reservoir outflow</h4>" + svg_lines([(n, None, list(zip(t, v))) for n, v in o["_reservoir_out"].items()],
                                                              "hours", f"flow ({o['flow units']})", h=200)
    if cfg.get("model"):
        body += "<h4>Map - minimum pressure and maximum velocity</h4>" + svg_map(cfg, results=True)
    if o["_energy"]:
        body += "<h4>Pump energy (from .out)</h4>" + _rows_html(list(o["_energy"][0].keys()),
                                                              [[e["pump"]] + [g(v) for k, v in e.items() if k != "pump"] for e in o["_energy"]])
    nm = sorted(o["_node"].items(), key=lambda kv: kv[1]["pmin"])
    body += _blk("Junction pressures", _rows_html(["Junction", "Min pressure", "Time of min", "Max pressure", "Max demand", "Quality min", "Quality max"],
                 [[n, g(v["pmin"]), hhmm(v["pmin_t"]), g(v["pmax"]), g(v["dmax"]), g(v["qmin"]), g(v["qmax"])] for n, v in nm]),
                 f"{len(nm)} junctions, lowest first", open_=bool(o.get("_negative")))
    lk = sorted(o["_link"].items(), key=lambda kv: -kv[1]["vmax"])
    body += _blk("Link maxima", _rows_html(["Link", "Type", "Max |flow|", "Time", "Max velocity", "Max head loss", "Closed periods"],
                 [[n, v["type"], g(v["qmax"]), hhmm(v["qmax_t"]), g(v["vmax"]), g(v["hlmax"]), v["closed_periods"]] for n, v in lk]),
                 f"{len(lk)} links, fastest first", open_=False)
    return body


def render_html(cfg) -> str:
    nav, body = [], []
    for i, m in enumerate(cfg["models"]):
        sec = [f'<section id="m{i}"><h2>{_e(m["name"])}']
        if m.get("error"):
            sec.append(' <span class="badge w">not readable</span>')
        elif m.get("only_out"):
            sec.append(' <span class="badge o">output file only</span>')
        elif m.get("rpt") or m.get("out"):
            bad = (m.get("rpt") or {}).get("errors") or (m.get("rpt") or {}).get("warnings") or (m.get("out") or {}).get("_negative")
            sec.append(f' <span class="badge{" w" if bad else ""}">{"run with problems" if bad else "run OK"}</span>')
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
           " · ".join(x for x in (m0["summary"].get("Units", "").split(" - ")[0], m0["summary"].get("Simulation")) if x))
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
    ap = argparse.ArgumentParser(prog=prog, description="Show the configuration of an EPANET 2.x model.")
    ap.add_argument("model", help=".inp file, a folder of .inp files, or a binary .out file")
    ap.add_argument("--rpt", help="report file (default: <model>.rpt)")
    ap.add_argument("--out", help="binary output file (default: <model>.out or .bin)")
    ap.add_argument("--all", action="store_true", help="terminal: every row of every section, and the rules")
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
