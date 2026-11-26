"""
modelvue.rasvue - Read a HEC-RAS project and show its configuration.

Parses the project file (.prj) and the files it lists - geometry (.gNN: river
reaches, cross sections, bridges, culverts, inline and lateral structures,
junctions, storage areas, 2D flow areas, SA/2D connections, boundary condition
lines, break lines, pump stations), steady flow (.fNN), unsteady flow (.uNN:
boundary and initial conditions, hydrographs, gate openings, meteorology), plans
(.pNN: simulation window, intervals, solver settings) and the RAS Mapper file
(.rasmap) - checks the data, and reports what a run left behind (computation
messages, volume accounting, plan HDF results).

The text files are read with the standard library. With h5py installed the
geometry HDF (.gNN.hdf) adds 2D mesh statistics and the plan HDF (.pNN.hdf) a
results summary. HEC-DSS files (.dss) are binary and are not read.

Usage
-----
    rasvue  path/to/project.prj
    rasvue  path/to/project_folder          # every HEC-RAS project in the folder (and below)
    rasvue  folder --project Muncie         # pick one project when several are found
    rasvue  project.prj --all               # terminal: every cross section and boundary
    rasvue  project.prj --json config.json  # also dump to JSON
    rasvue  project.prj --html              # webpage only, opens in the browser
    rasvue  project.prj --html report.html --no-browser
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from ._common import _blk, _e, _kv, _table, _two, h1, h2, row, strip_private
from ._report import MapFrame, PALETTE, find_ci, g, legend, num, page, read_text, svg_lines, text_table, write_page, frame

NODE_TYPES = {1: "Cross section", 2: "Culvert", 3: "Bridge", 4: "Multiple opening", 5: "Inline structure", 6: "Lateral structure",
              7: "Lateral structure"}
FILE_KINDS = [("Geom File", "geometry", "g"), ("Flow File", "steady flow", "f"), ("Unsteady File", "unsteady flow", "u"),
              ("QuasiSteady File", "quasi-unsteady flow", "q"), ("Plan File", "plan", "p"), ("Sediment File", "sediment", "s"),
              ("Water Quality File", "water quality", "w"), ("Hydro File", "hydrologic", "h")]
STEADY_BC = {"0": "none", "1": "known WS", "2": "critical depth", "3": "normal depth", "4": "rating curve"}
UNSTEADY_BC = ["Flow Hydrograph", "Stage Hydrograph", "Stage and Flow Hydrograph", "Lateral Inflow Hydrograph",
               "Uniform Lateral Inflow Hydrograph", "Ground Water Interflow", "Precipitation Hydrograph", "Rating Curve",
               "Friction Slope", "Gate Openings", "Elev Controlled Gate", "Navigation Dam", "Rule Operation",
               "Observed Stage and Flow Hydrograph", "IB Stage/Flow"]
BC_COUNTED = {"Flow Hydrograph": 1, "Stage Hydrograph": 1, "Lateral Inflow Hydrograph": 1, "Uniform Lateral Inflow Hydrograph": 1,
              "Precipitation Hydrograph": 1, "Ground Water Interflow": 1, "Gate Openings": 1, "Rating Curve": 2,
              "Stage and Flow Hydrograph": 1}
LOC_FIELDS = ["river", "reach", "rs", "rs_dn", "connection", "area", "pump", "bc_line", "extra"]
OBJECT_KEYS = {"Junct Name", "Storage Area", "Connection", "BC Line Name", "BreakLine Name", "Pump Station"}
BC_LABEL = {"Friction Slope": "Normal depth", "Gate Openings": "Gate openings (time series)", "Elev Controlled Gate": "Elevation-controlled gate",
            "Rule Operation": "Rules", "IB Stage/Flow": "Internal stage/flow"}
INTERVAL_SEC = {"SEC": 1, "MIN": 60, "HOUR": 3600, "DAY": 86400, "WEEK": 604800}


# --------------------------------------------------------------------------- #
# Low-level helpers
# --------------------------------------------------------------------------- #
def kv(line: str):
    if "=" not in line:
        return None, None
    k, v = line.split("=", 1)
    return k.strip(), v


def csv(v: str) -> list[str]:
    return [x.strip() for x in (v or "").split(",")]


def fw(lines, i, count, width, per_line=None):
    """Read ``count`` fixed-width numbers (``width`` characters each) from lines[i:]; returns (values, next i).

    Blank fields are skipped - unless ``per_line`` is given: then every line holds exactly that many fields
    (RAS trims trailing blanks) and a blank field is kept as None (e.g. the #Mann table).
    """
    vals = []
    while i < len(lines) and len(vals) < count:
        ln = lines[i].rstrip("\r\n")
        if "=" in ln and not re.match(r"^[\s\d.eE+-]*$", ln.split("=")[0]):
            break
        if per_line:
            ln = ln.ljust(per_line * width)
            for j in range(per_line):
                c = ln[j * width:(j + 1) * width].strip()
                vals.append(num(c) if c else None)
        else:
            for j in range(0, len(ln), width):
                c = ln[j:j + width].strip()
                if c:
                    vals.append(num(c))
        i += 1
    return vals[:count], i


def pairs(vals):
    return [(vals[k], vals[k + 1]) for k in range(0, len(vals) - 1, 2) if vals[k] is not None and vals[k + 1] is not None]


def count_of(v):
    x = num(csv(v)[0])
    return int(x) if x is not None else 0


def parse_interval(s):
    m = re.match(r"\s*(\d+(?:\.\d+)?)\s*([A-Za-z]+)", s or "")
    if not m:
        return None
    unit = m.group(2).upper()
    for k, sec in INTERVAL_SEC.items():
        if unit.startswith(k):
            return float(m.group(1)) * sec
    return None


def ras_datetime(d, t="0000"):
    try:
        t = (t or "0000").strip().replace(":", "")
        base = _dt.datetime.strptime(d.strip().upper(), "%d%b%Y")
        return base + _dt.timedelta(hours=int(t[:2] or 0), minutes=int(t[2:4] or 0))
    except (ValueError, AttributeError):
        return None


def _h5py():
    try:
        import h5py  # noqa: F401
        return h5py
    except ImportError:
        return None


def _dec(x):
    if isinstance(x, bytes):
        return x.decode("utf-8", "replace").strip()
    if hasattr(x, "tolist"):
        x = x.tolist()
        if isinstance(x, list):
            return [_dec(y) for y in x] if len(x) <= 8 else f"array of {len(x)}"
    if isinstance(x, float):
        return round(x, 6)
    return x


def poly_len(P):
    return sum(math.hypot(P[k + 1][0] - P[k][0], P[k + 1][1] - P[k][1]) for k in range(len(P) - 1))


def poly_area(P):
    if len(P) < 3:
        return 0.0
    return abs(sum(P[k][0] * P[k + 1][1] - P[k + 1][0] * P[k][1] for k in range(-1, len(P) - 1))) / 2


# --------------------------------------------------------------------------- #
# Project file
# --------------------------------------------------------------------------- #
def is_ras_project(p: Path) -> bool:
    try:
        with open(p, "rb") as f:
            head = f.read(200).decode("latin-1", "replace")
        return head.lstrip().startswith("Proj Title=")
    except OSError:
        return False


def read_prj(path: Path) -> dict:
    L = read_text(path).splitlines()
    d = {"title": None, "current_plan": None, "units": None, "files": {k: [] for _, k, _ in FILE_KINDS}, "settings": {},
         "description": ""}
    in_desc = False
    for ln in L:
        s = ln.strip()
        if s == "BEGIN DESCRIPTION:":
            in_desc = True
            continue
        if s == "END DESCRIPTION:":
            in_desc = False
            continue
        if in_desc:
            d["description"] += ln + "\n"
            continue
        if s in ("English Units", "SI Units"):
            d["units"] = s.replace(" Units", "")
            continue
        k, v = kv(ln)
        if k is None:
            continue
        v = v.strip()
        if k == "Proj Title":
            d["title"] = v
        elif k == "Current Plan":
            d["current_plan"] = v
        else:
            kind = next((name for key, name, _ in FILE_KINDS if key == k), None)
            if kind:
                d["files"][kind].append(v)
            elif v:
                d["settings"][k] = v
    d["description"] = d["description"].strip()
    return d


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def read_geometry(path: Path, warnings: list) -> dict:
    L = read_text(path).splitlines()
    G = {"file": path.name, "title": None, "version": None, "reaches": [], "nodes": [], "junctions": [], "storage": [],
         "connections": [], "bclines": [], "breaklines": [], "pumps": [], "rasters": [], "lateral": []}
    cur = reach = None
    n = len(L)
    i = 0
    while i < n:
        ln = L[i]
        i += 1
        s = ln.strip()
        if "=" not in ln:
            # gate tables: header line "IW Gate Name Wd,H,..." then a data line starting with the gate name
            if cur is not None and re.match(r"^(IW |LW |Conn |)Gate Name", s):
                cur.setdefault("gates", []).append(csv(L[i])[0] if i < n else "")
            continue
        k, v = kv(ln)
        if k == "Geom Title":
            G["title"] = v.strip()
        elif k == "Program Version":
            G["version"] = v.strip()
        elif k == "River Reach":
            f = csv(v)
            reach = {"river": f[0], "reach": f[1] if len(f) > 1 else "", "_xy": []}
            G["reaches"].append(reach)
            cur = reach
        elif k == "Reach XY" and reach is not None:
            vals, i = fw(L, i, 2 * count_of(v), 16)
            reach["_xy"] = pairs(vals)
        elif k == "Type RM Length L Ch R":
            f = csv(v)
            t = int(num(f[0], 0) or 0)
            cur = {"type": t, "kind": NODE_TYPES.get(t, f"type {t}"), "river": reach["river"] if reach else "",
                   "reach": reach["reach"] if reach else "", "rs": f[1] if len(f) > 1 else "",
                   "len_lob": num(f[2]) if len(f) > 2 else None, "len_ch": num(f[3]) if len(f) > 3 else None,
                   "len_rob": num(f[4]) if len(f) > 4 else None}
            G["nodes"].append(cur)
        elif cur is None and k not in OBJECT_KEYS:  # keys before the first object (pure-2D geometries start with a storage area)
            continue
        elif k == "Node Name":
            cur["name"] = v.strip()
        elif k == "XS GIS Cut Line":
            vals, i = fw(L, i, 2 * count_of(v), 16)
            cur["_cut"] = pairs(vals)
        elif k == "#Sta/Elev":
            vals, i = fw(L, i, 2 * count_of(v), 8)
            cur["_se"] = pairs(vals)
        elif k == "#Mann":
            nm = count_of(v)
            vals, i = fw(L, i, 3 * nm, 8, per_line=9)
            cur["_mann"] = [(vals[k3], vals[k3 + 1]) for k3 in range(0, len(vals) - 1, 3)]
        elif k == "Bank Sta":
            f = csv(v)
            cur["banks"] = (num(f[0]), num(f[1]) if len(f) > 1 else None)
        elif k == "Exp/Cntr":
            cur["exp_cntr"] = v.strip()
        elif k == "Levee":
            cur["levee"] = True
        elif k == "#XS Ineff":
            cur["ineffective"] = count_of(v)
        elif k == "#Block Obstruct":
            cur["blocked"] = count_of(v)
        elif k in ("Culvert", "Multiple Barrel Culv", "Connection Culv", "IW Culv", "LW Culv"):
            cur.setdefault("culverts", []).append(csv(v)[12] if len(csv(v)) > 12 and csv(v)[12] else csv(v)[0])
        elif k == "Pier Skew, UpSta & Num, DnSta & Num":
            cur["piers"] = cur.get("piers", 0) + 1
        elif k == "Lateral Weir Pos":
            cur["lw_pos"] = v.strip()
        elif k == "Lateral Weir End":
            f = csv(v)
            cur["lw_to"] = " / ".join(x for x in f if x)
        elif k == "Lateral Weir WD":
            cur["weir_width"] = num(v)
        elif k == "Lateral Weir Coef":
            cur["weir_coef"] = num(v)
        elif k == "Lateral Weir Distance":
            cur["lw_dist"] = num(v)
        elif k == "Lateral Weir Centerline":
            vals, i = fw(L, i, 2 * count_of(v), 16)
            cur["_line"] = pairs(vals)
        elif k == "Junct Name":
            cur = {"name": v.strip(), "up": [], "dn": []}
            G["junctions"].append(cur)
            reach = None
        elif k == "Junct X Y & Text X Y":
            f = csv(v)
            cur["x"], cur["y"] = num(f[0]), num(f[1]) if len(f) > 1 else None
        elif k == "Up River,Reach":
            cur.setdefault("up", []).append(" / ".join(csv(v)[:2]))
        elif k == "Dn River,Reach":
            cur.setdefault("dn", []).append(" / ".join(csv(v)[:2]))
        elif k == "Storage Area":
            f = csv(v)
            cur = {"name": f[0], "x": num(f[1]) if len(f) > 1 else None, "y": num(f[2]) if len(f) > 2 else None, "is2d": False}
            G["storage"].append(cur)
            reach = None
        elif k == "Storage Area Surface Line":
            m = count_of(v)
            pts = []
            for _ in range(m):
                if i >= n:
                    break
                a, b = num(L[i][0:16]), num(L[i][16:32])
                i += 1
                if a is not None and b is not None:
                    pts.append((a, b))
            cur["_poly"] = pts
        elif k == "Storage Area Type":
            cur["sa_type"] = v.strip()
        elif k == "Storage Area Area":
            cur["sa_area"] = v.strip()
        elif k == "Storage Area Min Elev":
            cur["min_elev"] = num(v)
        elif k == "Storage Area Is2D":
            cur["is2d"] = v.strip() in ("-1", "1", "True")
        elif k == "Storage Area Point Generation Data":
            f = csv(v)
            cur["spacing"] = " x ".join(x for x in f[2:4] if x)
        elif k == "Storage Area 2D Points":
            cur["cells_text"] = count_of(v)
        elif k == "Storage Area Mannings":
            cur["mannings"] = num(v)
        elif k == "Storage Area Vol Elev":
            cur["vol_elev_points"] = count_of(v)
        elif k == "Connection":
            f = csv(v)
            cur = {"name": f[0], "x": num(f[1]) if len(f) > 1 else None, "y": num(f[2]) if len(f) > 2 else None}
            G["connections"].append(cur)
            reach = None
        elif k == "Connection Desc":
            cur["description"] = v.strip()
        elif k == "Connection Line":
            vals, i = fw(L, i, 2 * count_of(v), 16)
            cur["_line"] = pairs(vals)
        elif k == "Connection Up SA":
            cur["up"] = v.strip()
        elif k == "Connection Dn SA":
            cur["dn"] = v.strip()
        elif k == "Conn Routing Type":
            cur["routing"] = v.strip()
        elif k == "Conn Weir WD":
            cur["weir_width"] = num(v)
        elif k == "Conn Weir Coef":
            cur["weir_coef"] = num(v)
        elif k == "BC Line Name":
            cur = {"name": v.strip()}
            G["bclines"].append(cur)
            reach = None
        elif k == "BC Line Storage Area":
            cur["area"] = v.strip()
        elif k == "BC Line Arc":
            vals, i = fw(L, i, 2 * count_of(v), 16)
            cur["_line"] = pairs(vals)
        elif k == "BreakLine Name":
            cur = {"name": v.strip()}
            G["breaklines"].append(cur)
            reach = None
        elif k == "BreakLine CellSize Min":
            cur["cell_min"] = v.strip()
        elif k == "BreakLine CellSize Max":
            cur["cell_max"] = v.strip()
        elif k == "BreakLine Polyline":
            vals, i = fw(L, i, 2 * count_of(v), 16)
            cur["_line"] = pairs(vals)
        elif k == "Pump Station":
            f = csv(v)
            cur = {"name": f[0], "x": num(f[1]) if len(f) > 1 else None, "y": num(f[2]) if len(f) > 2 else None}
            G["pumps"].append(cur)
            reach = None
        elif k == "Pump Station From":
            cur["from"] = " / ".join(x for x in csv(v) if x)
        elif k == "Pump Station To":
            cur["to"] = " / ".join(x for x in csv(v)[:4] if x)
        elif k == "Geom Raster":
            G["rasters"].append(csv(v)[0])

    # ---- derived values and checks -------------------------------------
    for nd in G["nodes"]:
        se = nd.get("_se") or []
        if se:
            st = [p[0] for p in se]
            zs = [p[1] for p in se]
            nd["points"] = len(se)
            nd["min_elev"] = min(zs)
            nd["max_elev"] = max(zs)
            nd["width"] = st[-1] - st[0]
            if any(st[k + 1] < st[k] for k in range(len(st) - 1)):
                warnings.append(f"{G['file']}: XS {nd['river']} / {nd['reach']} RS {nd['rs']}: stations decrease along the section")
            b = nd.get("banks")
            if b and b[0] is not None and b[1] is not None:
                if not (st[0] <= b[0] <= b[1] <= st[-1]):
                    warnings.append(f"{G['file']}: XS {nd['river']} / {nd['reach']} RS {nd['rs']}: bank stations {g(b[0])}, {g(b[1])} "
                                    f"outside the section ({g(st[0])} .. {g(st[-1])})")
                il = min(range(len(st)), key=lambda q: abs(st[q] - b[0]))
                ir = min(range(len(st)), key=lambda q: abs(st[q] - b[1]))
                nd["bank_elev"] = (zs[il], zs[ir])
        mann = [m for m in (nd.get("_mann") or []) if m[0] is not None and m[1] is not None]
        if se and mann:
            # drop zero-width segments (n change points that coincide, e.g. at a bank station on the section edge)
            ends = [m[0] for m in mann[1:]] + [se[-1][0]]
            mann = [m for m, e in zip(mann, ends) if e - m[0] > 1e-9] or mann
            nd["_mann"] = mann
        if mann:
            ns = [m[1] for m in mann]
            if min(ns) <= 0 or max(ns) > 0.3:
                warnings.append(f"{G['file']}: XS {nd['river']} / {nd['reach']} RS {nd['rs']}: unusual Manning's n "
                                f"({g(min(ns))} .. {g(max(ns))}; expected between 0 and 0.3)")

            def n_at(x):  # n of the segment that starts at or before station x
                cand = [m[1] for m in mann if m[0] <= x + 1e-9]
                return cand[-1] if cand else mann[0][1]
            b = nd.get("banks")
            if b and b[0] is not None and b[1] is not None and se:
                if b[0] > se[0][0]:
                    nd["n_lob"] = n_at(se[0][0])
                nd["n_ch"] = n_at((b[0] + b[1]) / 2)
                if b[1] < se[-1][0]:
                    nd["n_rob"] = n_at(b[1] + 1e-6)
            elif len(mann) == 1:
                nd["n_ch"] = mann[0][1]
    for r in G["reaches"]:
        xs = [nd for nd in G["nodes"] if nd["river"] == r["river"] and nd["reach"] == r["reach"]]
        r["n_xs"] = sum(1 for nd in xs if nd["type"] == 1)
        r["n_struct"] = sum(1 for nd in xs if nd["type"] != 1)
        r["length_ch"] = sum(nd.get("len_ch") or 0 for nd in xs if nd["type"] == 1)
        r["length_xy"] = poly_len(r["_xy"]) if r["_xy"] else None
        rss = [nd["rs"] for nd in xs if nd["type"] == 1]
        r["rs_range"] = f"{rss[0]} .. {rss[-1]}" if rss else ""
        if r["n_xs"] < 2:
            warnings.append(f"{G['file']}: reach {r['river']} / {r['reach']} has fewer than 2 cross sections")
    for sa in G["storage"]:
        sa["polygon_area"] = poly_area(sa.get("_poly") or []) or None
    G["counts"] = {
        "River reaches": len(G["reaches"]), "Cross sections": sum(1 for nd in G["nodes"] if nd["type"] == 1),
        "Bridges": sum(1 for nd in G["nodes"] if nd["type"] == 3), "Culverts": sum(1 for nd in G["nodes"] if nd["type"] == 2),
        "Multiple openings": sum(1 for nd in G["nodes"] if nd["type"] == 4),
        "Inline structures": sum(1 for nd in G["nodes"] if nd["type"] == 5),
        "Lateral structures": sum(1 for nd in G["nodes"] if nd["type"] in (6, 7)),
        "Junctions": len(G["junctions"]), "Storage areas": sum(1 for s in G["storage"] if not s["is2d"]),
        "2D flow areas": sum(1 for s in G["storage"] if s["is2d"]), "SA/2D connections": len(G["connections"]),
        "BC lines": len(G["bclines"]), "Break lines": len(G["breaklines"]), "Pump stations": len(G["pumps"]),
    }
    G["counts"] = {k: v for k, v in G["counts"].items() if v}
    return G


def geometry_hdf(path: Path) -> dict | None:
    h5py = _h5py()
    if h5py is None or not path.exists():
        return None
    out = {"file": path.name, "areas": []}
    try:
        with h5py.File(path, "r") as f:
            out["attrs"] = {k: _dec(v) for k, v in f.attrs.items()}
            if "Geometry" in f:
                out["geometry"] = {k: _dec(v) for k, v in f["Geometry"].attrs.items()}
            base = "Geometry/2D Flow Areas"
            if base in f and "Attributes" in f[base]:
                at = f[base]["Attributes"][()]
                names = at.dtype.names or ()
                for rec in at:
                    d = {nm: _dec(rec[nm]) for nm in names}
                    name = d.get("Name")
                    grp = f[base].get(name) if name else None
                    if grp is not None:
                        d.update({f"mesh {k.lower()}": _dec(v) for k, v in grp.attrs.items() if k.startswith("Cell")})
                        if "Cells Center Coordinate" in grp:
                            d["cells incl. ghost"] = grp["Cells Center Coordinate"].shape[0]
                        if "Faces FacePoint Indexes" in grp:
                            d["faces"] = grp["Faces FacePoint Indexes"].shape[0]
                        if "Perimeter" in grp:
                            d["_perimeter"] = [tuple(map(float, p)) for p in grp["Perimeter"][()]]
                    out["areas"].append(d)
            lc = "Geometry/Land Cover (Manning's n)/Calibration Table"
            if lc in f:
                t = f[lc][()]
                out["landcover"] = [{nm: _dec(r[nm]) for nm in t.dtype.names} for r in t]
    except Exception as e:  # damaged / unexpected file
        out["error"] = str(e)
    return out


# --------------------------------------------------------------------------- #
# Flow files
# --------------------------------------------------------------------------- #
def read_steady(path: Path) -> dict:
    L = read_text(path).splitlines()
    d = {"file": path.name, "title": None, "profiles": [], "flows": [], "boundaries": [], "storage": []}
    nprof = 1
    i, n = 0, len(L)
    bnd = None
    while i < n:
        k, v = kv(L[i])
        i += 1
        if k is None:
            continue
        if k == "Flow Title":
            d["title"] = v.strip()
        elif k == "Program Version":
            d["version"] = v.strip()
        elif k == "Number of Profiles":
            nprof = count_of(v) or 1
        elif k == "Profile Names":
            d["profiles"] = csv(v)
        elif k == "River Rch & RM":
            f = csv(v)
            vals, i = fw(L, i, nprof, 8)
            d["flows"].append({"river": f[0], "reach": f[1] if len(f) > 1 else "", "rs": f[2] if len(f) > 2 else "", "q": vals})
        elif k == "Boundary for River Rch & Prof#":
            f = csv(v)
            bnd = {"river": f[0], "reach": f[1] if len(f) > 1 else "", "profile": f[2] if len(f) > 2 else ""}
            d["boundaries"].append(bnd)
        elif bnd is not None and k in ("Up Type", "Dn Type"):
            bnd[k] = STEADY_BC.get(v.strip(), v.strip())
        elif bnd is not None and k in ("Up Slope", "Dn Slope", "Up Known WS", "Dn Known WS", "Up Rating Curve", "Dn Rating Curve"):
            bnd[k] = v.strip()
        elif k == "Storage Area Elev":
            f = csv(v)
            vals, i = fw(L, i, nprof, 8)
            d["storage"].append({"name": f[0], "elev": vals})
    if not d["profiles"]:
        d["profiles"] = [f"PF {k + 1}" for k in range(nprof)]
    return d


def read_unsteady(path: Path) -> dict:
    L = read_text(path).splitlines()
    d = {"file": path.name, "title": None, "boundaries": [], "initial": [], "met": {}, "settings": {}}
    i, n = 0, len(L)
    b = None
    while i < n:
        k, v = kv(L[i])
        i += 1
        if k is None:
            continue
        val = v.strip()
        if k == "Flow Title":
            d["title"] = val
        elif k == "Program Version":
            d["version"] = val
        elif k == "Boundary Location":
            f = csv(v)
            loc = {LOC_FIELDS[j]: f[j] for j in range(min(len(f), len(LOC_FIELDS))) if f[j]}
            b = {"location": loc, "type": None, "settings": {}}
            d["boundaries"].append(b)
        elif k in ("Initial Flow Loc", "Initial Storage Elev", "Initial RRR Elev"):
            d["initial"].append({"kind": k.replace("Initial ", ""), "value": " / ".join(x for x in csv(v) if x)})
        elif k.startswith("Met BC") or k in ("Precipitation Mode", "Wind Mode", "Air Density Mode", "Met Point Raster Parameters"):
            if k.startswith("Met BC"):
                var, _, prop = val.partition("|")
                prop, _, pv = prop.partition("=")
                if pv:
                    d["met"].setdefault(var, {})[prop] = pv
            elif val:
                d["met"][k] = val
        elif b is not None and (k in UNSTEADY_BC or k.startswith("Gate Name")):
            if k.startswith("Gate Name"):
                b.setdefault("gates", []).append(val)
                b["type"] = b["type"] or "Gate Openings"
                continue
            b["type"] = b["type"] or k
            if k in BC_COUNTED:
                m = count_of(v)
                vals, i = fw(L, i, BC_COUNTED[k] * m, 8)
                if k == "Rating Curve":
                    b["_curve"] = pairs(vals)
                    b["points"] = m
                elif k == "Gate Openings":
                    b.setdefault("_gates", []).append(vals)
                    b["points"] = m
                else:
                    b["_values"] = vals
                    b["points"] = m
            elif k == "Friction Slope":
                b["slope"] = csv(v)[0]
            else:
                b["settings"][k] = val
        elif b is not None and k in ("Interval", "DSS Path", "Use DSS", "DSS File", "Use Fixed Start Time", "Fixed Start Date/Time",
                                     "Is Critical Boundary", "Gate DSS Path", "Gate Use DSS", "Gate Time Interval", "Flow Hydrograph QMult",
                                     "Flow Hydrograph Slope", "Stage Hydrograph TW Check", "Flow Hydrograph Min"):
            if val and val not in (",",):
                b["settings"][k] = val
        elif b is None and val and k not in ("Restart Filename",):
            d["settings"][k] = val
        elif k == "Use Restart" or k == "Restart Filename":
            d["settings"][k] = val
    for bc in d["boundaries"]:
        vals = bc.get("_values")
        if vals:
            fin = [x for x in vals if x is not None]
            if fin:
                bc["min"], bc["max"] = min(fin), max(fin)
        st = bc["settings"]
        if st.get("Use DSS", "").lower() == "true":
            bc["source"] = "DSS " + st.get("DSS Path", "")
        elif bc.get("points"):
            bc["source"] = f"table ({bc['points']} values" + (f", every {st['Interval']}" if st.get("Interval") else "") + ")"
        elif bc.get("slope"):
            bc["source"] = f"slope {bc['slope']}"
        bc["where"] = loc_text(bc["location"])
        if bc["type"]:
            bc["type"] = BC_LABEL.get(bc["type"], bc["type"])
    return d


def loc_text(loc: dict) -> str:
    parts = []
    if loc.get("river"):
        parts.append(f"{loc['river']} / {loc.get('reach', '')} RS {loc.get('rs', '')}" + (f" .. {loc['rs_dn']}" if loc.get("rs_dn") else ""))
    if loc.get("area"):
        parts.append(f"area {loc['area']}")
    if loc.get("connection"):
        parts.append(f"connection {loc['connection']}")
    if loc.get("pump"):
        parts.append(f"pump station {loc['pump']}")
    if loc.get("bc_line"):
        parts.append(f"BC line {loc['bc_line']}")
    return ", ".join(parts)


# --------------------------------------------------------------------------- #
# Plans and run output
# --------------------------------------------------------------------------- #
def read_plan(path: Path) -> dict:
    L = read_text(path).splitlines()
    d = {"file": path.name, "ext": path.suffix[1:].lower(), "title": None, "short_id": None, "regime": None, "settings": {}}
    for ln in L:
        s = ln.strip()
        if s in ("Subcritical Flow", "Supercritical Flow", "Mixed Flow"):
            d["regime"] = s.replace(" Flow", "").lower()
            continue
        k, v = kv(ln)
        if k is None:
            continue
        val = v.strip()
        if k == "Plan Title":
            d["title"] = val
        elif k == "Short Identifier":
            d["short_id"] = val
        elif k == "Program Version":
            d["version"] = val
        elif k == "Geom File":
            d["geom"] = val
        elif k == "Flow File":
            d["flow"] = val
        elif k == "Simulation Date":
            f = csv(v)
            if len(f) >= 4:
                d["start"], d["end"] = f"{f[0]} {f[1]}", f"{f[2]} {f[3]}"
                t0, t1 = ras_datetime(f[0], f[1]), ras_datetime(f[2], f[3])
                if t0 and t1:
                    d["_t0"] = t0
                    d["duration_h"] = round((t1 - t0).total_seconds() / 3600, 3)
        elif val and not k.startswith(("Stage Flow Hydrograph", "Profile ", "Encroach", "Global ", "Lateral Flow Hydrograph")):
            d["settings"][k] = val
    st = d["settings"]
    flow = d.get("flow", "")
    d["kind"] = {"u": "unsteady", "f": "steady", "q": "quasi-unsteady"}.get(flow[:1].lower(), "?")
    if st.get("Run Sediment", "0").strip() not in ("0", ""):
        d["kind"] += " + sediment"
    d["key"] = {k: st.get(k) for k in ("Computation Interval", "Output Interval", "Instantaneous Interval", "Mapping Interval") if st.get(k)}
    for k in ("UNET D2 Equation", "UNET Theta", "UNET D2 Theta", "UNET D2 Cores", "UNET D1 Cores", "UNET MxIter", "UNET ZTol", "UNET D2 Z Tol",
              "Std Step Tol", "Critical Tol", "Friction Slope Method"):
        if st.get(k):
            d["key"][k] = st[k]
    d["runs"] = {k.replace("Run ", ""): st[k] for k in ("Run HTab", "Run UNet", "Run Sediment", "Run PostProcess", "Run WQNet", "Run RASMapper")
                 if k in st and st[k] not in ("0", "")}
    return d


def read_compute_messages(folder: Path, stem: str, plan_ext: str) -> dict | None:
    for name in (f"{stem}.{plan_ext}.computeMsgs.txt", f"{stem}.{plan_ext}.comp_msgs.txt"):
        p = find_ci(folder, name)
        if not p:
            continue
        lines = []
        for ln in read_text(p).splitlines():
            if ln.startswith("Progress |"):
                continue
            lines.append(re.sub(r"^Message\s*\|\s?", "", ln).rstrip())
        text = "\n".join(lines)
        out = {"file": p.name, "_text": re.sub(r"\n{3,}", "\n\n", text)[-30000:]}
        m = re.search(r"Simulation started at:\s*(.+)", text)
        out["started"] = m.group(1).strip() if m else None
        out["finished"] = bool(re.search(r"Finished (Unsteady Flow Simulation|Steady Flow Simulation|Post Processing|Simulation)|"
                                         r"Complete Process", text))
        out["errors"] = sorted({ln.strip() for ln in lines if re.search(r"\berror\b|\bfailed\b|unstable|went unstable", ln, re.I)
                                and not re.search(r"Max(imum)? error|ERROR\s*$|error tol", ln, re.I)})[:20]
        m = re.search(r"Overall Volume Accounting Error.*?([-\d.Ee+]+)\s*%?", text)
        out["volume_error"] = m.group(0).strip() if m else None
        m = re.search(r"Complete Process\s+([^\n]+)", text)
        out["run_time"] = m.group(1).strip() if m else None
        out["max_iterations"] = len(re.findall(r"Maximum iterations", text))
        return out
    return None


def plan_hdf(path: Path) -> dict | None:
    h5py = _h5py()
    if not path.exists():
        return None
    if h5py is None:
        return {"file": path.name, "note": "install h5py to read the plan results HDF"}
    out = {"file": path.name}
    try:
        with h5py.File(path, "r") as f:
            out["attrs"] = {k: _dec(v) for k, v in f.attrs.items()}
            for grp in ("Plan Data/Plan Information", "Results/Unsteady/Summary", "Results/Unsteady/Summary/Volume Accounting",
                        "Results/Steady/Summary"):
                if grp in f:
                    out[grp] = {k: _dec(v) for k, v in f[grp].attrs.items()}
            for ds in ("Results/Summary/Compute Messages (text)", "Results/Unsteady/Summary/Compute Messages (text)",
                       "Results/Steady/Summary/Compute Messages (text)"):
                if ds in f:
                    raw = f[ds][()]
                    raw = raw.tobytes() if hasattr(raw, "tobytes") else raw
                    out["_messages"] = (raw if isinstance(raw, bytes) else str(raw).encode()).decode("utf-8", "replace").strip("\x00 ")[-20000:]
                    break
            base = "Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas"
            if base in f:
                areas = []
                for name, grp in f[base].items():
                    if "Maximum Water Surface" not in grp:
                        continue
                    arr = grp["Maximum Water Surface"][()]
                    wse = arr[0] if getattr(arr, "ndim", 1) == 2 else arr
                    rec = {"area": name, "cells": int(len(wse))}
                    vals = [float(x) for x in wse if math.isfinite(float(x))]
                    if vals:
                        rec["max WSE range"] = f"{g(min(vals))} .. {g(max(vals))}"
                    zp = f"Geometry/2D Flow Areas/{name}/Cells Minimum Elevation"
                    if zp in f:
                        z = f[zp][()]
                        dep = [float(a) - float(b) for a, b in zip(wse, z) if math.isfinite(float(a)) and math.isfinite(float(b))]
                        wet = [x for x in dep if x > 0.01]
                        rec["wet cells (max depth > 0.01)"] = len(wet)
                        if wet:
                            rec["max depth"] = g(max(wet))
                    areas.append(rec)
                out["areas_2d"] = areas
            listing = []

            def visit(name, obj):
                if len(listing) < 80 and isinstance(obj, h5py.Dataset) and name.startswith("Results"):
                    listing.append((name, "x".join(map(str, obj.shape))))
            f.visititems(visit)
            out["datasets"] = listing
    except Exception as e:
        out["error"] = str(e)
    return out


def read_rasmap(path: Path) -> dict | None:
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError) as e:
        return {"file": path.name, "error": str(e)}
    out = {"file": path.name, "version": root.findtext("Version")}
    pr = root.find("RASProjectionFilename")
    out["projection"] = pr.get("Filename") if pr is not None else None
    groups = {}
    for grp in root:
        if grp.tag in ("Version", "RASProjectionFilename"):
            continue
        items = [(lay.get("Name") or lay.get("Type"), lay.get("Type"), lay.get("Filename")) for lay in grp.findall("Layer")]
        items = [x for x in items if x[0]]
        if items:
            groups[grp.tag] = items
    out["groups"] = groups
    t = root.find("Units")
    out["units"] = t.text if t is not None else None
    return out


# --------------------------------------------------------------------------- #
# Project
# --------------------------------------------------------------------------- #
def find_projects(arg: str) -> list[Path]:
    p = Path(arg).expanduser().resolve()
    if not p.exists():
        sys.exit(f"Path not found: {p}")
    if p.is_file():
        if p.suffix.lower() == ".prj" and is_ras_project(p):
            return [p]
        cands = [q for q in sorted(p.parent.glob("*.prj")) + sorted(p.parent.glob("*.PRJ")) if is_ras_project(q)]
        stem = [q for q in cands if q.stem.lower() == p.name.split(".")[0].lower()]
        if stem or cands:
            return (stem or cands)[:1]
        sys.exit(f"{p.name} is not a HEC-RAS project file (.prj starting with 'Proj Title=')")
    fs = sorted({q for q in list(p.rglob("*.prj")) + list(p.rglob("*.PRJ")) if is_ras_project(q)})
    if not fs:
        sys.exit(f"No HEC-RAS project (.prj) found under {p}")
    return fs


def build_project(path: Path) -> dict:
    folder, stem = path.parent, path.stem
    warnings = []
    prj = read_prj(path)

    def load(ext, reader):
        p = find_ci(folder, f"{stem}.{ext}")
        if not p:
            warnings.append(f"{stem}.{ext} listed in the project but not found")
            return {"file": f"{stem}.{ext}", "ext": ext, "missing": True}
        try:
            d = reader(p)
        except Exception as e:
            warnings.append(f"{p.name}: could not be read ({e})")
            d = {"file": p.name, "error": str(e)}
        d["ext"] = ext
        return d

    geoms = [load(e, lambda p: read_geometry(p, warnings)) for e in prj["files"]["geometry"]]
    for gm in geoms:
        if gm.get("missing") or gm.get("error"):
            continue
        hp = find_ci(folder, f"{stem}.{gm['ext']}.hdf")
        hd = geometry_hdf(hp) if hp else None
        if hd:
            gm["hdf"] = hd
    steady = [load(e, read_steady) for e in prj["files"]["steady flow"]]
    unsteady = [load(e, read_unsteady) for e in prj["files"]["unsteady flow"]]
    plans = [load(e, read_plan) for e in prj["files"]["plan"]]
    exts = {x["ext"] for x in geoms + steady + unsteady if not x.get("missing")}
    exts |= {e for e in prj["files"]["quasi-unsteady flow"] if find_ci(folder, f"{stem}.{e}")}
    for pl in plans:
        if pl.get("missing") or pl.get("error"):
            continue
        for key in ("geom", "flow"):
            if pl.get(key) and pl[key].lower() not in {e.lower() for e in exts}:
                warnings.append(f"plan {pl['ext']} ({pl.get('title')}): {key} file {pl[key]} is not in the project")
        msgs = read_compute_messages(folder, stem, pl["ext"])
        if msgs:
            pl["messages"] = msgs
            for e in msgs["errors"][:3]:
                warnings.append(f"plan {pl['ext']} computation: {e}")
        hp = find_ci(folder, f"{stem}.{pl['ext']}.hdf")
        if hp:
            pl["results_hdf"] = plan_hdf(hp)
        old = find_ci(folder, f"{stem}.O{pl['ext'][1:]}")
        if old:
            pl["binary_output"] = old.name
        pl["computed"] = bool(msgs or hp or old)
    others = {}
    for kind in ("quasi-unsteady flow", "sediment", "water quality", "hydrologic"):
        for e in prj["files"][kind]:
            p = find_ci(folder, f"{stem}.{e}")
            title = None
            if p:
                m = re.search(r"^(?:Flow|Sediment|Water Quality|WQ|Hydro)\s*Title=(.*)$", read_text(p), re.M)
                title = m.group(1).strip() if m else None
            others.setdefault(kind, []).append({"ext": e, "file": p.name if p else f"{stem}.{e}", "title": title, "missing": not p})
    rm = find_ci(folder, f"{stem}.rasmap")
    rasmap = read_rasmap(rm) if rm else None
    if prj["current_plan"] and prj["current_plan"].lower() not in {pl["ext"] for pl in plans}:
        warnings.append(f"current plan {prj['current_plan']} is not in the plan list")
    cur = next((pl for pl in plans if pl.get("ext") == (prj["current_plan"] or "").lower()), None)
    tot = {}
    for gm in geoms:
        for k, v in gm.get("counts", {}).items():
            tot[k] = max(tot.get(k, 0), v)
    summary = {
        "Project title": prj["title"], "Units": prj["units"], "Current plan": f"{prj['current_plan']} - {cur.get('title')}" if cur else prj["current_plan"],
        "Geometries": len(geoms), "Steady flow files": len(steady), "Unsteady flow files": len(unsteady), "Plans": len(plans),
        "Computed plans": sum(1 for p in plans if p.get("computed")),
        "Largest geometry": ", ".join(f"{k}: {v}" for k, v in tot.items()) or None,
        "RAS version (plans)": ", ".join(sorted({p.get("version") for p in plans if p.get("version")})) or None,
        "DSS file": prj["settings"].get("DSS File"),
    }
    files = sorted(q.name for q in folder.iterdir() if q.is_file() and q.name.lower().startswith(stem.lower() + "."))
    dss = [{"file": q.name, "size_kb": round(q.stat().st_size / 1024, 1)} for q in sorted(folder.glob("*.dss"))]
    return {"name": stem, "path": str(path), "folder": str(folder), "title": prj["title"], "units": prj["units"],
            "description": prj["description"], "current_plan": prj["current_plan"], "settings": prj["settings"], "summary": summary,
            "geometries": geoms, "steady": steady, "unsteady": unsteady, "plans": plans, "other_files": others, "rasmap": rasmap,
            "dss": dss, "files": files, "warnings": warnings}


def build_config(arg: str, only: str | None = None) -> dict:
    files = find_projects(arg)
    if only:
        sel = [f for f in files if only.lower() in (f.stem.lower(), f.parent.name.lower())]
        files = sel or files
    root = Path(arg).expanduser().resolve()
    return {"root": str(root if root.is_dir() else root.parent), "name": root.stem if root.is_file() else root.name,
            "projects": [build_project(f) for f in files]}


# --------------------------------------------------------------------------- #
# Terminal
# --------------------------------------------------------------------------- #
def show_project(p, full=False):
    h1(f"HEC-RAS PROJECT: {p['name']}")
    row("File", p["path"])
    h2("Summary")
    for k, v in p["summary"].items():
        row(k, v)
    if p["description"]:
        print("  | " + p["description"].replace("\n", "\n  | ")[:800])
    if p["plans"]:
        h2("Plans")
        text_table(["Plan", "Title", "Kind", "Geometry", "Flow", "Simulation window", "Computed"],
                   [[pl.get("ext"), pl.get("title", ""), pl.get("kind", ""), pl.get("geom", ""), pl.get("flow", ""),
                     f"{pl.get('start', '')} - {pl.get('end', '')}" if pl.get("start") else "", "yes" if pl.get("computed") else "-"]
                    for pl in p["plans"]], maxw=34)
        for pl in p["plans"]:
            if full or pl.get("ext") == (p["current_plan"] or "").lower():
                print(f"\n  Plan {pl.get('ext')}: {pl.get('title')}  ({pl.get('kind')}, {pl.get('regime') or ''})")
                for k, v in pl.get("key", {}).items():
                    row(k, v, indent=4)
                row("runs", pl.get("runs"), indent=4)
                ms = pl.get("messages")
                if ms:
                    row("computation", ("finished" if ms["finished"] else "did not finish") + (f", {ms['run_time']}" if ms.get("run_time") else ""), indent=4)
                    row("volume accounting", ms.get("volume_error"), indent=4)
                    for e in ms["errors"][:5]:
                        print(f"      ! {e}")
                rh = pl.get("results_hdf")
                if rh:
                    row("results HDF", rh.get("note") or rh.get("error") or rh["file"], indent=4)
                    for a in rh.get("areas_2d", []):
                        row(f"2D area {a['area']}", {k: v for k, v in a.items() if k != "area"}, indent=6)
    for gm in p["geometries"]:
        h2(f"Geometry {gm.get('ext')}: {gm.get('title') or ''} ({gm.get('file')})")
        if gm.get("missing") or gm.get("error"):
            print(f"  ! {gm.get('error') or 'file missing'}")
            continue
        row("Contents", gm["counts"])
        if gm["reaches"]:
            text_table(["River", "Reach", "XS", "Structures", "RS range", "Channel length"],
                       [[r["river"], r["reach"], r["n_xs"], r["n_struct"], r["rs_range"], g(r["length_ch"])] for r in gm["reaches"]], indent=4)
        st = [nd for nd in gm["nodes"] if nd["type"] != 1]
        if st:
            print("    Structures")
            text_table(["Type", "River / reach", "RS", "Name", "Detail"],
                       [[nd["kind"], f"{nd['river']} / {nd['reach']}", nd["rs"], nd.get("name", ""), struct_detail(nd)] for nd in st],
                       indent=6, limit=None if full else 20, maxw=36)
        if full:
            xs = [nd for nd in gm["nodes"] if nd["type"] == 1]
            print("    Cross sections")
            text_table(["River", "Reach", "RS", "Points", "Min elev", "Banks", "n (LOB/Ch/ROB)", "L Ch"],
                       [[nd["river"], nd["reach"], nd["rs"], nd.get("points", ""), g(nd.get("min_elev")),
                         "/".join(g(b) for b in (nd.get("banks") or ())), "/".join(g(x) for x in (nd.get("n_lob"), nd.get("n_ch"), nd.get("n_rob")) if x is not None),
                         g(nd.get("len_ch"))] for nd in xs], indent=6)
        for sa in gm["storage"]:
            hd = next((a for a in (gm.get("hdf") or {}).get("areas", []) if a.get("Name") == sa["name"]), {})
            cells = hd.get("Cell Count") or sa.get("cells_text")
            row(("2D flow area " if sa["is2d"] else "storage area ") + sa["name"],
                (f"{cells} cells" if sa["is2d"] and cells else "") + (f", spacing {sa['spacing']}" if sa.get("spacing") else "")
                + (f", n {g(sa.get('mannings'))}" if sa["is2d"] and sa.get("mannings") else ""), indent=4)
        for c in gm["connections"]:
            row(f"connection {c['name']}", f"{c.get('up', '')} -> {c.get('dn', '')}", indent=4)
    for f_ in p["steady"]:
        h2(f"Steady flow {f_.get('ext')}: {f_.get('title') or ''}")
        if f_.get("missing") or f_.get("error"):
            continue
        row("profiles", f_["profiles"])
        text_table(["River", "Reach", "RS"] + f_["profiles"][:6],
                   [[q["river"], q["reach"], q["rs"]] + [g(x) for x in q["q"][:6]] for q in f_["flows"]], indent=4, limit=None if full else 15)
        for b in f_["boundaries"]:
            row(f"{b['river']} / {b['reach']} prof {b['profile']}",
                f"up: {b.get('Up Type', '')} {b.get('Up Slope', '') or b.get('Up Known WS', '')}  dn: {b.get('Dn Type', '')} "
                f"{b.get('Dn Slope', '') or b.get('Dn Known WS', '')}", indent=4)
    for u in p["unsteady"]:
        h2(f"Unsteady flow {u.get('ext')}: {u.get('title') or ''}")
        if u.get("missing") or u.get("error"):
            continue
        text_table(["Type", "Location", "Source", "Range"],
                   [[b["type"] or "(none)", b["where"], b.get("source", ""), f"{g(b['min'])} .. {g(b['max'])}" if "min" in b else ""]
                    for b in u["boundaries"]], indent=4, maxw=48, limit=None if full else 25)
        if u["initial"]:
            row("initial conditions", f"{len(u['initial'])} ({', '.join(sorted({x['kind'] for x in u['initial']}))})", indent=4)
        if u["met"].get("Precipitation Mode") and u["met"]["Precipitation Mode"] != "Disable":
            row("precipitation", u["met"]["Precipitation Mode"], indent=4)
    if p["warnings"]:
        h2("Warnings")
        for w in p["warnings"][:60 if not full else None]:
            print(f"  ! {w}")
        if len(p["warnings"]) > 60 and not full:
            print(f"  ... {len(p['warnings']) - 60} more (use --all)")


def struct_detail(nd) -> str:
    bits = []
    if nd.get("culverts"):
        bits.append(f"{len(nd['culverts'])} culvert(s)")
    if nd.get("gates"):
        bits.append(f"{len(nd['gates'])} gate group(s)")
    if nd.get("piers"):
        bits.append(f"{nd['piers']} pier(s)")
    if nd.get("weir_width"):
        bits.append(f"weir {g(nd['weir_width'])} wide" + (f", Cd {g(nd['weir_coef'])}" if nd.get("weir_coef") else ""))
    if nd.get("lw_to"):
        bits.append(f"to {nd['lw_to']}")
    return ", ".join(bits)


def show(cfg, full=False):
    for p in cfg["projects"]:
        show_project(p, full)
    print()


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #
EXTRA_CSS = """
svg.map .rch{fill:none;stroke:#1565c0;stroke-width:2.6} svg.map .xs{fill:none;stroke:#2e7d32;stroke-width:1;opacity:.75}
svg.map .xs.st{stroke:#c62828;stroke-width:2.5;opacity:1} svg.map .sa{fill:#ffcc80;fill-opacity:.45;stroke:#e65100;stroke-width:1.2}
svg.map .a2d{fill:#90caf9;fill-opacity:.35;stroke:#0d47a1;stroke-width:1.4} svg.map .conn{fill:none;stroke:#6a1b9a;stroke-width:2.5}
svg.map .bcl{fill:none;stroke:#d32f2f;stroke-width:4} svg.map .brk{fill:none;stroke:#455a64;stroke-width:1.2;stroke-dasharray:5 3}
svg.map .lw{fill:none;stroke:#ef6c00;stroke-width:2.5} svg.map .jn{fill:#6a1b9a;stroke:#fff}
svg.plot .bank{stroke:#c62828;stroke-dasharray:3 3} svg.plot .grd{fill:#d7ccc8;fill-opacity:.6;stroke:#6d4c41;stroke-width:1.4}
"""


def geometry_map(gm) -> str:
    hd_per = {a.get("Name"): a.get("_perimeter") for a in (gm.get("hdf") or {}).get("areas", [])}
    pts = []
    for r in gm["reaches"]:
        pts += r["_xy"]
    for nd in gm["nodes"]:
        pts += nd.get("_cut") or []
    for sa in gm["storage"]:
        pts += sa.get("_poly") or hd_per.get(sa["name"]) or []
    for c in gm["connections"] + gm["bclines"] + gm["breaklines"]:
        pts += c.get("_line") or []
    if len(pts) < 2:
        return '<p class="muted">The geometry holds no georeferenced coordinates - no map.</p>'
    # schematic-only models use tiny coordinates (0..1); still draw them
    F = MapFrame(pts)
    o = []
    for sa in gm["storage"]:
        poly = sa.get("_poly") or hd_per.get(sa["name"])
        if poly:
            o.append(f'<polygon points="{F.pts(poly)}" class="{"a2d" if sa["is2d"] else "sa"}"><title>{"2D flow area" if sa["is2d"] else "storage area"} '
                     f'{_e(sa["name"])}</title></polygon>')
        elif sa.get("x") is not None:
            o.append(f'<rect x="{F.X(sa["x"]) - 6:.1f}" y="{F.Y(sa["y"]) - 6:.1f}" width="12" height="12" class="sa"><title>storage area {_e(sa["name"])}</title></rect>')
    for b in gm["breaklines"]:
        if b.get("_line"):
            o.append(f'<polyline points="{F.pts(b["_line"])}" class="brk"><title>break line {_e(b["name"])}</title></polyline>')
    for nd in gm["nodes"]:
        cut = nd.get("_cut")
        if cut:
            cls = "xs" if nd["type"] == 1 else "xs st"
            o.append(f'<polyline points="{F.pts(cut)}" class="{cls}"><title>{_e(nd["kind"])} {_e(nd["river"])} / {_e(nd["reach"])} RS {_e(nd["rs"])}'
                     f'{" " + _e(nd["name"]) if nd.get("name") else ""}</title></polyline>')
        if nd.get("_line"):
            o.append(f'<polyline points="{F.pts(nd["_line"])}" class="lw"><title>lateral structure RS {_e(nd["rs"])}</title></polyline>')
    for r in gm["reaches"]:
        if r["_xy"]:
            o.append(f'<polyline points="{F.pts(r["_xy"])}" class="rch"><title>{_e(r["river"])} / {_e(r["reach"])}: {r["n_xs"]} XS</title></polyline>')
    for c in gm["connections"]:
        if c.get("_line"):
            o.append(f'<polyline points="{F.pts(c["_line"])}" class="conn"><title>connection {_e(c["name"])}: {_e(c.get("up", ""))} -> {_e(c.get("dn", ""))}</title></polyline>')
    for b in gm["bclines"]:
        if b.get("_line"):
            o.append(f'<polyline points="{F.pts(b["_line"])}" class="bcl"><title>BC line {_e(b["name"])} ({_e(b.get("area", ""))})</title></polyline>')
    for j in gm["junctions"]:
        if j.get("x") is not None and j.get("y") is not None:
            o.append(f'<circle cx="{F.X(j["x"]):.1f}" cy="{F.Y(j["y"]):.1f}" r="5" class="jn"><title>junction {_e(j["name"])}</title></circle>')
    for pmp in gm["pumps"]:
        if pmp.get("x") is not None:
            o.append(f'<rect x="{F.X(pmp["x"]) - 5:.1f}" y="{F.Y(pmp["y"]) - 5:.1f}" width="10" height="10" fill="#c2185b"><title>pump station {_e(pmp["name"])}</title></rect>')
    items = [("border-color:#1565c0", "river reach", "ln"), ("border-color:#2e7d32;border-top-width:1px", "cross section", "ln")]
    if any(nd["type"] != 1 and nd.get("_cut") for nd in gm["nodes"]):
        items.append(("border-color:#c62828", "structure", "ln"))
    if any(s["is2d"] for s in gm["storage"]):
        items.append(("background:#90caf9;border:1px solid #0d47a1", "2D flow area", "bx"))
    if any(not s["is2d"] for s in gm["storage"]):
        items.append(("background:#ffcc80;border:1px solid #e65100", "storage area", "bx"))
    if gm["connections"]:
        items.append(("border-color:#6a1b9a", "SA/2D connection", "ln"))
    if gm["bclines"]:
        items.append(("border-color:#d32f2f", "BC line", "ln"))
    if gm["breaklines"]:
        items.append(("border-color:#455a64;border-top-style:dashed;border-top-width:1px", "break line", "ln"))
    if any(nd.get("_line") for nd in gm["nodes"]):
        items.append(("border-color:#ef6c00", "lateral structure", "ln"))
    if gm["junctions"]:
        items.append(("background:#6a1b9a", "junction", "dot"))
    return F.svg(o) + legend(items)


def xs_plot(nd, units, w=420, h=210) -> str:
    se = nd.get("_se")
    if not se or len(se) < 2:
        return ""
    st = [p[0] for p in se]
    zs = [p[1] for p in se]
    z0, z1 = min(zs), max(zs)
    pad = (z1 - z0) * 0.08 or 1
    X, Y, o = frame(w, h, 52, 10, 10, 32, st[0], st[-1], z0 - pad, z1 + pad, f"station ({_len(units)})", f"elevation ({_len(units)})")
    P = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in se)
    o.append(f'<polygon points="{X(st[0]):.1f},{Y(z0 - pad):.1f} {P} {X(st[-1]):.1f},{Y(z0 - pad):.1f}" class="grd"/>')
    b = nd.get("banks")
    if b:
        for x in b:
            if x is not None and st[0] <= x <= st[-1]:
                o.append(f'<line x1="{X(x):.1f}" x2="{X(x):.1f}" y1="10" y2="{h - 32}" class="bank"/>')
    title = f"{nd['river']} / {nd['reach']} RS {nd['rs']}" + (f" ({nd['name']})" if nd.get("name") else "")
    n = "/".join(g(x) for x in (nd.get("n_lob"), nd.get("n_ch"), nd.get("n_rob")) if x is not None)
    return (f'<div><b style="font-size:12px">{_e(title)}</b><span class="muted" style="font-size:12px"> · n {n}</span>'
            f'<svg class="plot" viewBox="0 0 {w} {h}">' + "".join(o) + "</svg></div>")


def _len(units):
    return "m" if units == "SI" else "ft"


def profile_plot(gm, r, units) -> str:
    xs = [nd for nd in gm["nodes"] if nd["river"] == r["river"] and nd["reach"] == r["reach"] and nd["type"] == 1 and nd.get("min_elev") is not None]
    if len(xs) < 2:
        return ""
    dist, d = [], 0.0
    for k, nd in enumerate(xs):
        dist.append(d)
        d += nd.get("len_ch") or 0
    tot = dist[-1]
    xd = [tot - x for x in dist]  # distance from the downstream end
    series = [("thalweg", "#1565c0", xd, [nd["min_elev"] for nd in xs])]
    lb = [nd.get("bank_elev", (None, None))[0] for nd in xs]
    rb = [nd.get("bank_elev", (None, None))[1] for nd in xs]
    if any(v is not None for v in lb):
        series.append(("left bank", "#2e7d32", xd, lb))
    if any(v is not None for v in rb):
        series.append(("right bank", "#c62828", xd, rb))
    return svg_lines(series, f"channel distance from downstream end ({_len(units)})", f"elevation ({_len(units)})", h=240)


def _geometry_html(gm, units) -> str:
    if gm.get("missing") or gm.get("error"):
        return f'<p class="warn">{_e(gm.get("error") or "file not found")}</p>'
    out = [f'<p class="muted">{_e(gm.get("title") or "")} · <code>{_e(gm["file"])}</code>'
           + (f' · RAS {_e(gm["version"])}' if gm.get("version") else "") + "</p>",
           '<div class="stats">' + "".join(f'<div class="stat"><b>{v}</b><span>{_e(k)}</span></div>' for k, v in gm["counts"].items()) + "</div>",
           geometry_map(gm)]
    if gm["reaches"]:
        out.append("<h4>River reaches</h4>" + _table(["River", "Reach", "Cross sections", "Structures", "RS range", f"Channel length ({_len(units)})",
                                                      f"Centerline length ({_len(units)})"],
                   [[_e(r["river"]), _e(r["reach"]), r["n_xs"], r["n_struct"], _e(r["rs_range"]), g(r["length_ch"]), g(r["length_xy"])] for r in gm["reaches"]]))
        prof = "".join(f"<p><b>{_e(r['river'])} / {_e(r['reach'])}</b></p>" + profile_plot(gm, r, units) for r in gm["reaches"][:12])
        out.append(_blk("Longitudinal profiles (thalweg and bank elevations)", prof, open_=True))
        sample = []
        for r in gm["reaches"][:6]:
            xs = [nd for nd in gm["nodes"] if nd["river"] == r["river"] and nd["reach"] == r["reach"] and nd["type"] == 1 and nd.get("_se")]
            if xs:
                pick = {0, len(xs) // 2, len(xs) - 1}
                sample += [xs[k] for k in sorted(pick)]
        if sample:
            out.append(_blk("Cross sections (upstream, middle, downstream of each reach)",
                            '<div class="grid2">' + "".join(xs_plot(nd, units) for nd in sample[:18]) + "</div>", open_=False))
        xs_all = [nd for nd in gm["nodes"] if nd["type"] == 1]
        rows = [[_e(nd["river"]), _e(nd["reach"]), _e(nd["rs"]), _e(nd.get("name", "")), nd.get("points", ""), g(nd.get("min_elev")),
                 "/".join(g(b) for b in (nd.get("banks") or ())), "/".join(g(x) for x in (nd.get("n_lob"), nd.get("n_ch"), nd.get("n_rob")) if x is not None),
                 g(nd.get("len_lob")), g(nd.get("len_ch")), g(nd.get("len_rob")), _e(nd.get("exp_cntr", "")),
                 ("levee " if nd.get("levee") else "") + (f"{nd['ineffective']} ineff. " if nd.get("ineffective") else "")
                 + (f"{nd['blocked']} blocked" if nd.get("blocked") else "")] for nd in xs_all[:3000]]
        out.append(_blk("Cross sections", '<input class="pfilter" placeholder="filter..." style="margin:6px 0;padding:4px 8px">' +
                        _table(["River", "Reach", "RS", "Name", "Points", "Min elev", "Banks", "n LOB/Ch/ROB", "L LOB", "L Ch", "L ROB",
                                "Exp/Cntr", "Other"], rows), f"{len(xs_all)}" + (" (first 3000)" if len(xs_all) > 3000 else ""), open_=False))
    st = [nd for nd in gm["nodes"] if nd["type"] != 1]
    if st:
        out.append(_blk("Structures", _table(["Type", "River", "Reach", "RS", "Name", "Detail"],
                        [[_e(nd["kind"]), _e(nd["river"]), _e(nd["reach"]), _e(nd["rs"]), _e(nd.get("name", "")),
                          f'<span class="wrap">{_e(struct_detail(nd))}</span>'] for nd in st]), f"{len(st)}"))
    if gm["junctions"]:
        out.append(_blk("Junctions", _table(["Name", "Upstream reaches", "Downstream reaches"],
                        [[_e(j["name"]), _e(", ".join(j["up"])), _e(", ".join(j["dn"]))] for j in gm["junctions"]]), f"{len(gm['junctions'])}"))
    if gm["storage"]:
        hd = {a.get("Name"): a for a in (gm.get("hdf") or {}).get("areas", [])}
        rows = []
        for sa in gm["storage"]:
            a = hd.get(sa["name"], {})
            rows.append([("2D flow area" if sa["is2d"] else "storage area"), _e(sa["name"]),
                         a.get("Cell Count") or sa.get("cells_text") or "", _e(sa.get("spacing", "")), g(sa.get("mannings")) if sa["is2d"] else "",
                         g(a.get("mesh cell average size")), g(sa.get("polygon_area")), g(sa.get("min_elev")),
                         _e(sa.get("sa_type", "")) if not sa["is2d"] else ""])
        out.append(_blk("Storage areas and 2D flow areas", _table(["Kind", "Name", "Cells", "Cell spacing", "Manning's n", "Avg cell area",
                                                                    "Polygon area", "Min elev", "SA type"], rows),
                        "cell counts from the geometry HDF" if hd else ""))
    if gm["connections"]:
        out.append(_blk("SA/2D connections", _table(["Name", "From", "To", "Routing", "Weir width", "Weir coef", "Culverts", "Gates"],
                        [[_e(c["name"]), _e(c.get("up", "")), _e(c.get("dn", "")), _e(c.get("routing", "")), g(c.get("weir_width")),
                          g(c.get("weir_coef")), len(c.get("culverts", [])), len(c.get("gates", []))] for c in gm["connections"]])))
    if gm["bclines"] or gm["breaklines"]:
        out.append(_blk("BC lines and break lines", _table(["Kind", "Name", "2D area / cell size", "Length"],
                        [["BC line", _e(b["name"]), _e(b.get("area", "")), g(poly_len(b.get("_line") or []))] for b in gm["bclines"]]
                        + [["break line", _e(b["name"]), _e(f"{b.get('cell_min', '')} .. {b.get('cell_max', '')}".strip(" .")),
                            g(poly_len(b.get("_line") or []))] for b in gm["breaklines"]])))
    if gm["pumps"]:
        out.append(_blk("Pump stations", _table(["Name", "From", "To"], [[_e(x["name"]), _e(x.get("from", "")), _e(x.get("to", ""))] for x in gm["pumps"]])))
    hd = gm.get("hdf")
    if hd:
        body = _kv(list((hd.get("geometry") or {}).items()) + list((hd.get("attrs") or {}).items()))
        if hd.get("areas"):
            cols = [c for c in hd["areas"][0] if not c.startswith("_")]
            body += _table(cols, [[_e(g(a.get(c)) if isinstance(a.get(c), float) else a.get(c)) for c in cols] for a in hd["areas"]])
        if hd.get("landcover"):
            cols = list(hd["landcover"][0])
            body += "<h4>Land cover Manning's n</h4>" + _table(cols, [[_e(r.get(c)) for c in cols] for r in hd["landcover"]])
        if hd.get("error"):
            body += f'<p class="warn">{_e(hd["error"])}</p>'
        out.append(_blk(f"Geometry HDF · {_e(hd['file'])}", body, open_=False))
    if gm["rasters"]:
        out.append(_blk("Background layers", " ".join(f"<code>{_e(x)}</code>" for x in gm["rasters"]), open_=False))
    return "".join(out)


def _steady_html(f_) -> str:
    if f_.get("missing") or f_.get("error"):
        return f'<p class="warn">{_e(f_.get("error") or "file not found")}</p>'
    out = [_kv([("File", f_["file"]), ("Title", f_.get("title")), ("Profiles", ", ".join(f_["profiles"]))])]
    out.append("<h4>Flows</h4>" + _table(["River", "Reach", "RS"] + f_["profiles"],
                                        [[_e(q["river"]), _e(q["reach"]), _e(q["rs"])] + [g(x) for x in q["q"]] for q in f_["flows"]]))
    if f_["boundaries"]:
        out.append("<h4>Boundary conditions</h4>" + _table(["River", "Reach", "Profile", "Upstream", "Downstream"],
                   [[_e(b["river"]), _e(b["reach"]), _e(b["profile"]), _e(f"{b.get('Up Type', '')} {b.get('Up Slope', '') or b.get('Up Known WS', '')}"),
                     _e(f"{b.get('Dn Type', '')} {b.get('Dn Slope', '') or b.get('Dn Known WS', '')}")] for b in f_["boundaries"]]))
    if f_["storage"]:
        out.append("<h4>Storage area elevations</h4>" + _table(["Storage area"] + f_["profiles"],
                   [[_e(s["name"])] + [g(x) for x in s["elev"]] for s in f_["storage"]]))
    return "".join(out)


def _unsteady_html(u, plans) -> str:
    if u.get("missing") or u.get("error"):
        return f'<p class="warn">{_e(u.get("error") or "file not found")}</p>'
    out = [_two([("File", u["file"]), ("Title", u.get("title")), ("RAS version", u.get("version")),
                 ("Used by plans", ", ".join(f"{p['ext']}" for p in plans if (p.get("flow") or "").lower() == u["ext"].lower()))]
                + [(k, v) for k, v in u["settings"].items() if not k.startswith("Non-Newtonian")][:12])]
    def extra(b):  # settings not already shown in the Source column, and not at their defaults
        skip = {"Interval", "DSS Path", "Use DSS", "DSS File", "Gate DSS Path", "Gate Use DSS"}
        dflt = {"Use Fixed Start Time": "False", "Fixed Start Date/Time": ",", "Is Critical Boundary": "False", "Stage Hydrograph TW Check": "0"}
        return ", ".join(f"{k}: {v}" for k, v in b["settings"].items() if k not in skip and dflt.get(k) != v)

    out.append("<h4>Boundary conditions</h4>" + _table(["Type", "Location", "Source", "Min", "Max", "Other settings"],
               [[_e(b["type"] or "(none)"), f'<span class="wrap">{_e(b["where"])}</span>',
                 f'<span class="wrap">{_e(b.get("source", ""))}</span>', g(b.get("min")), g(b.get("max")),
                 f'<span class="wrap">{_e(extra(b))}</span>'] for b in u["boundaries"]]))
    t0 = next((p.get("_t0") for p in plans if (p.get("flow") or "").lower() == u["ext"].lower() and p.get("_t0")), None)
    plots = []
    for k, b in enumerate([b for b in u["boundaries"] if b.get("_values")][:12]):
        dt = parse_interval(b["settings"].get("Interval", "")) or 3600
        hrs = [j * dt / 3600 for j in range(len(b["_values"]))]
        unit = "stage" if "Stage" in (b["type"] or "") else ("precipitation" if "Precip" in (b["type"] or "") else "flow")
        plots.append(f"<p><b>{_e(b['type'])}</b> · {_e(b['where'])}</p>" +
                     svg_lines([(b["type"], PALETTE[k % len(PALETTE)], hrs, b["_values"])],
                               "hours from simulation start" + (f" ({t0:%d%b%Y %H:%M})" if t0 else ""), unit, h=200,
                               step="Precip" in (b["type"] or "")))
    if plots:
        out.append(_blk("Hydrographs (tables in the file)", "".join(plots), open_=False))
    for b in u["boundaries"]:
        if b.get("_curve"):
            c = b["_curve"]
            out.append(_blk(f"Rating curve · {_e(b['where'])}", svg_lines([("rating curve", "#6a1b9a", [p[1] for p in c], [p[0] for p in c])], "flow", "stage", h=200), open_=False))
    if u["initial"]:
        out.append(_blk("Initial conditions", _table(["Kind", "Location / value"], [[_e(x["kind"]), _e(x["value"])] for x in u["initial"]]),
                        f"{len(u['initial'])}", open_=False))
    met = {k: v for k, v in u["met"].items() if not isinstance(v, dict)}
    metv = {k: v for k, v in u["met"].items() if isinstance(v, dict) and any(x not in ("0", "Nearest", "DSS") for x in v.values())}
    if met or metv:
        out.append(_blk("Meteorology", _kv(list(met.items()) + [(k, ", ".join(f"{a}: {b}" for a, b in v.items())) for k, v in metv.items()]), open_=False))
    return "".join(out)


def _plan_html(pl, p) -> str:
    if pl.get("missing") or pl.get("error"):
        return f'<p class="warn">{_e(pl.get("error") or "file not found")}</p>'
    out = [_two([("File", pl["file"]), ("Title", pl.get("title")), ("Short identifier", pl.get("short_id")), ("Kind", pl.get("kind")),
                 ("Flow regime", pl.get("regime")), ("Geometry", pl.get("geom")), ("Flow file", pl.get("flow")),
                 ("Simulation start", pl.get("start")), ("Simulation end", pl.get("end")), ("Duration (h)", g(pl.get("duration_h"))),
                 ("RAS version", pl.get("version")), ("Computations", ", ".join(pl.get("runs", {})))] + list(pl.get("key", {}).items()))]
    ms = pl.get("messages")
    if ms:
        body = _kv([("File", ms["file"]), ("Started", ms.get("started")), ("Finished", "yes" if ms["finished"] else "no"),
                    ("Run time", ms.get("run_time")), ("Volume accounting", ms.get("volume_error")),
                    ("'Maximum iterations' messages", ms.get("max_iterations") or None)])
        if ms["errors"]:
            body += '<ul class="warns">' + "".join(f"<li>{_e(e)}</li>" for e in ms["errors"]) + "</ul>"
        body += f'<pre class="log">{_e(ms["_text"])}</pre>'
        out.append(_blk("Computation messages", body, "finished" if ms["finished"] else "not finished", open_=bool(ms["errors"])))
    rh = pl.get("results_hdf")
    if rh:
        body = ""
        if rh.get("note") or rh.get("error"):
            body += f'<p class="{"warn" if rh.get("error") else "muted"}">{_e(rh.get("note") or rh.get("error"))}</p>'
        for grp in ("Plan Data/Plan Information", "Results/Unsteady/Summary", "Results/Unsteady/Summary/Volume Accounting", "Results/Steady/Summary"):
            if rh.get(grp):
                body += f"<h4>{_e(grp)}</h4>" + _kv(list(rh[grp].items()))
        if rh.get("areas_2d"):
            cols = list(rh["areas_2d"][0])
            body += "<h4>2D flow areas - maximum water surface</h4>" + _table(cols, [[_e(a.get(c, "")) for c in cols] for a in rh["areas_2d"]])
        if rh.get("_messages"):
            body += f'<h4>Compute messages (HDF)</h4><pre class="log">{_e(rh["_messages"])}</pre>'
        if rh.get("datasets"):
            body += _blk("Result datasets", _table(["Dataset", "Shape"], [[_e(a), _e(b)] for a, b in rh["datasets"]]), open_=False)
        out.append(_blk(f"Results HDF · {_e(rh['file'])}", body))
    if pl.get("binary_output"):
        out.append(f'<p class="muted">Binary output file <code>{_e(pl["binary_output"])}</code> (not read).</p>')
    if not pl.get("computed"):
        out.append('<p class="muted">No computation messages or results files for this plan in the project folder.</p>')
    out.append(_blk("All plan settings", '<input class="pfilter" placeholder="filter..." style="margin:6px 0;padding:4px 8px">' +
                    _table(["Key", "Value"], [[_e(k), f'<span class="wrap">{_e(v)}</span>'] for k, v in pl["settings"].items()]),
                    f"{len(pl['settings'])}", open_=False))
    return "".join(out)


def _project_html(p, pi) -> tuple[str, list]:
    present, out = [], []

    def add(bid, title, body, note=""):
        if body:
            present.append((bid, title))
            out.append(_blk(title, body, note, bid=bid))

    summ = _two(list(p["summary"].items()))
    if p["description"]:
        summ += f'<div class="desc">{_e(p["description"])}</div>'
    if p["plans"]:
        summ += "<h4>Plans</h4>" + _table(["Plan", "Title", "Kind", "Geometry", "Flow", "Start", "End", "Computed"],
                                          [[("<b>" + _e(pl.get("ext")) + "</b>") if pl.get("ext") == (p["current_plan"] or "").lower() else _e(pl.get("ext")),
                                            _e(pl.get("title", "")), _e(pl.get("kind", "")), _e(pl.get("geom", "")), _e(pl.get("flow", "")),
                                            _e(pl.get("start", "")), _e(pl.get("end", "")), "yes" if pl.get("computed") else ""] for pl in p["plans"]])
    add(f"p{pi}-sum", "Summary", summ)
    for k, pl in enumerate(p["plans"]):
        add(f"p{pi}-pl{k}", f"Plan {pl.get('ext')} · {pl.get('title') or ''}", _plan_html(pl, p), pl.get("kind", ""))
    for k, gm in enumerate(p["geometries"]):
        add(f"p{pi}-g{k}", f"Geometry {gm.get('ext')} · {gm.get('title') or ''}", _geometry_html(gm, p["units"]), gm.get("file", ""))
    for k, f_ in enumerate(p["steady"]):
        add(f"p{pi}-f{k}", f"Steady flow {f_.get('ext')} · {f_.get('title') or ''}", _steady_html(f_), f_.get("file", ""))
    for k, u in enumerate(p["unsteady"]):
        add(f"p{pi}-u{k}", f"Unsteady flow {u.get('ext')} · {u.get('title') or ''}", _unsteady_html(u, p["plans"]), u.get("file", ""))
    oth = [[_e(kind), _e(x["ext"]), _e(x["file"]), _e(x.get("title") or ""), "missing" if x["missing"] else ""] for kind, lst in p["other_files"].items() for x in lst]
    if oth:
        add(f"p{pi}-oth", "Other model files", _table(["Kind", "Ext", "File", "Title", ""], oth))
    rm = p.get("rasmap")
    if rm:
        body = _kv([("File", rm["file"]), ("Version", rm.get("version")), ("Projection", rm.get("projection")), ("Units", rm.get("units"))])
        for grp, items in rm.get("groups", {}).items():
            body += f"<h4>{_e(grp)}</h4>" + _table(["Name", "Type", "File"], [[_e(a), _e(b), _e(c or "")] for a, b, c in items])
        if rm.get("error"):
            body += f'<p class="warn">{_e(rm["error"])}</p>'
        add(f"p{pi}-map", "RAS Mapper", body, rm["file"])
    files = '<p class="wrap">' + " ".join(f"<code>{_e(f)}</code>" for f in p["files"]) + "</p>"
    if p["dss"]:
        files += "<h4>HEC-DSS files (not read)</h4>" + _table(["File", "Size (kB)"], [[_e(d["file"]), g(d["size_kb"])] for d in p["dss"]])
    if p["settings"]:
        files += _blk("Other project settings", _kv(list(p["settings"].items())), open_=False)
    add(f"p{pi}-files", "Project files", files, f"{len(p['files'])}")
    return "".join(out), present


def render_html(cfg) -> str:
    nav, body = [], []
    for i, p in enumerate(cfg["projects"]):
        comp = sum(1 for pl in p["plans"] if pl.get("computed"))
        badge = f'<span class="badge{"" if comp else " o"}">{comp}/{len(p["plans"])} plans computed</span>'
        if p["warnings"]:
            badge += f' <span class="badge w">{len(p["warnings"])} warning(s)</span>'
        sec = [f'<section id="p{i}"><h2>{_e(p["name"])} {badge}</h2><p class="muted">{_e(p["title"] or "")} · <code>{_e(p["path"])}</code></p>']
        if p["warnings"]:
            ws = p["warnings"]
            sec.append('<ul class="warns">' + "".join(f"<li>{_e(w)}</li>" for w in ws[:40]) +
                       (f"<li>... {len(ws) - 40} more</li>" if len(ws) > 40 else "") + "</ul>")
        h, present = _project_html(p, i)
        sec.append(h + "</section>")
        body.append("".join(sec))
        nav.append(f'<a href="#p{i}"><b>{_e(p["name"])}</b></a>')
        if len(cfg["projects"]) <= 12:
            nav += [f'<a class="sub" href="#{b}">{_e(t)}</a>' for b, t in present]
    p0 = cfg["projects"][0]
    sub = f'{len(cfg["projects"])} projects' if len(cfg["projects"]) > 1 else \
        " · ".join(x for x in (p0["title"], p0["units"] and f'{p0["units"]} units', p0["summary"].get("Current plan")) if x)
    return page("HEC-RAS", "rasvue", cfg["name"], sub, nav, body, cfg["root"], EXTRA_CSS)


def write_html(cfg, out: Path, open_browser=True):
    return write_page(render_html(cfg), out, open_browser)


# --------------------------------------------------------------------------- #
def main(argv=None, prog=None):
    ap = argparse.ArgumentParser(prog=prog, description="Show the configuration of a HEC-RAS project.")
    ap.add_argument("model", help="HEC-RAS project file (.prj), any file of the project, or a folder containing projects")
    ap.add_argument("--project", dest="only", metavar="NAME", help="when several projects are found, show only this one (file stem or folder name)")
    ap.add_argument("--all", action="store_true", help="terminal: every plan, cross section, boundary and warning")
    ap.add_argument("--json", metavar="FILE", help="also write the configuration to a JSON file")
    ap.add_argument("--html", nargs="?", const="", metavar="FILE",
                    help="write the configuration as a webpage and open it in the browser, instead of printing it "
                         "(default file: <project>_config.html next to the .prj)")
    ap.add_argument("--no-browser", action="store_true", help="with --html, write the page but don't open it")
    a = ap.parse_args(argv)

    cfg = build_config(a.model, a.only)
    if a.html is None:
        show(cfg, a.all)
    else:
        p0 = cfg["projects"][0]
        out = Path(a.html) if a.html else (Path(p0["folder"]) / f"{p0['name']}_config.html" if len(cfg["projects"]) == 1
                                           else Path(cfg["root"]) / f"{cfg['name']}_config.html")
        write_html(cfg, out, open_browser=not a.no_browser)
    if a.json:
        Path(a.json).write_text(json.dumps(strip_private(cfg), indent=2, default=str), encoding="utf-8")
        print(f"Configuration written to {a.json}")


if __name__ == "__main__":
    main()
