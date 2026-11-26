"""
modelvue.hmsvue - Read a HEC-HMS project and show its configuration.

Parses the project file (.hms) and the files it lists - basin models (.basin),
meteorologic models (.met), control specifications (.control), simulation runs
(.run), time-series gages (.gage), paired data (.pdata), grids (.grid) - checks
the element network, and, if the runs have been computed, reads the run results
(results/RUN_<name>.results: peak flow, volumes, observed-flow statistics) and the
run logs (notes, warnings, errors).

Time series stored in HEC-DSS files (.dss) are binary and are not read; their
file names and pathnames are listed. Standard library only.

Usage
-----
    hmsvue  path/to/project.hms
    hmsvue  path/to/project_folder          # every .hms in the folder (and below)
    hmsvue  project.hms --all               # terminal: every parameter of every element
    hmsvue  project.hms --json config.json  # also dump to JSON
    hmsvue  project.hms --html              # webpage only, opens in the browser
    hmsvue  project.hms --html report.html --no-browser
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from struct import error as struct_error

from ._common import _blk, _e, _kv, _table, _two, h1, h2, row, strip_private
from ._report import MapFrame, find_ci, g, legend, num, page, read_shp, read_text, shp_kind, svg_bars, text_table, write_page

ELEMENT_TYPES = ["Subbasin", "Reach", "Junction", "Reservoir", "Diversion", "Source", "Sink"]
PLURAL = {"Subbasin": "subbasins", "Reach": "reaches", "Junction": "junctions", "Reservoir": "reservoirs",
          "Diversion": "diversions", "Source": "sources", "Sink": "sinks"}
ELEMENT_STYLE = {"Subbasin": ("#2e7d32", "bx"), "Reach": ("#1565c0", "ln"), "Junction": ("#6a1b9a", "dot"),
                 "Reservoir": ("#00838f", "bx"), "Diversion": ("#ef6c00", "dot"), "Source": ("#558b2f", "dot"),
                 "Sink": ("#c62828", "dot")}
META_KEYS = {"Last Modified Date", "Last Modified Time", "Canvas X", "Canvas Y", "Label X", "Label Y", "From Canvas X",
             "From Canvas Y", "Downstream", "Description", "File", "Version", "Filepath Separator", "Allow Simultaneous Precip Et"}
SUBBASIN_METHODS = ["Canopy", "Surface", "LossRate", "Transform", "Baseflow"]
STAT_LABELS = [("Outflow Maximum", "Peak outflow"), ("Outflow Maximum Time", "Time of peak"), ("Outflow Volume", "Outflow volume"),
               ("Outflow Depth", "Outflow depth"), ("Precipitation Total", "Precipitation"), ("Loss Total", "Loss"),
               ("Excess Total", "Excess"), ("Baseflow Volume", "Baseflow volume")]
OBS_LABELS = [("Observed Flow Maximum", "Observed peak"), ("Observed Flow Maximum Time", "Observed peak time"),
              ("Observed Flow Volume", "Observed volume"), ("Observed Flow Nash Sutcliffe", "Nash-Sutcliffe"),
              ("Observed Flow Percent Bias", "Percent bias"), ("Observed Flow RMSE Stdev", "RMSE / StDev")]


# --------------------------------------------------------------------------- #
# HMS text format: "Type: name" blocks closed by "End:", "     Key: value" lines,
# nested sub-blocks "Key: ..." ... "End Key:" (conduits, pumps, variants, layers)
# --------------------------------------------------------------------------- #
_LINE = re.compile(r"^([^:]+?):\s?(.*)$")


def parse_blocks(text: str) -> list[dict]:
    blocks, cur, lines = [], None, text.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i].rstrip()
        i += 1
        if not ln.strip():
            continue
        if cur is None:
            if ln[0].isspace():
                continue
            m = _LINE.match(ln)
            if not m or m.group(1).strip() == "End":
                continue
            cur = {"type": m.group(1).strip(), "name": m.group(2).strip(), "items": [], "sub": []}
            end_own = f"End {cur['type']}:"
            body = []
            while i < len(lines) and not lines[i].startswith("End:") and not lines[i].startswith(end_own):
                body.append(lines[i])
                i += 1
            i += 1
            kv = [(m2.group(1).strip(), m2.group(2).strip()) for b in body if (m2 := _LINE.match(b.strip()))]

            def opens(k, v, c):  # "Conduit: Culvert" .. "End Conduit:";  "Map Layer:" (no value) .. "End Layer:"
                return k == c or (k.endswith(" " + c) and v == "")

            # "End X:" closes a sub-block only when a key opens one (so not "End Date:", "End Time:")
            closers = {k[4:] for k, _ in kv if k.startswith("End ")}
            closers = {c for c in closers if any(opens(k, v, c) for k, v in kv if not k.startswith("End "))}

            def opener(k, v):
                return next((c for c in closers if opens(k, v, c)), None)

            stack = []
            for b in body:
                s = b.strip()
                m2 = _LINE.match(s)
                if not s or not m2:
                    continue
                k, v = m2.group(1).strip(), m2.group(2).strip()
                if k.startswith("End ") and k[4:] in closers:
                    if stack:
                        stack.pop()
                    continue
                c = opener(k, v)
                if c and not (stack and stack[-1]["_closer"] == c):  # "Layer:0" inside "Map Layer:" is a value
                    sb = {"type": k, "name": v, "items": [], "_closer": c}
                    cur["sub"].append(sb)  # nested sub-blocks are flattened
                    stack.append(sb)
                    continue
                (stack[-1] if stack else cur)["items"].append((k, v))
            for sb in cur["sub"]:
                sb.pop("_closer", None)
            blocks.append(cur)
            cur = None
    return blocks


def bget(block, key, default=None):
    for k, v in block["items"]:
        if k == key:
            return v
    return default


def bgetall(block, key):
    return [v for k, v in block["items"] if k == key]


def bdict(block, skip=META_KEYS):
    d = {}
    for k, v in block["items"]:
        if k in skip:
            continue
        d[k] = v if k not in d else f"{d[k]}; {v}"
    return d


def of_type(blocks, t):
    return [b for b in blocks if b["type"] == t]


def _date(d, t=None):
    for fmt in ("%d %B %Y", "%d %b %Y", "%d%b%Y"):
        try:
            x = _dt.datetime.strptime(d.strip(), fmt)
            break
        except (ValueError, AttributeError):
            x = None
    if x is None:
        return None
    if t:
        hh, _, mm = t.strip().partition(":")
        try:
            x += _dt.timedelta(hours=int(hh), minutes=int(mm or 0))
        except ValueError:
            pass
    return x


def _path(folder: Path, name: str | None) -> Path | None:
    if not name:
        return None
    name = name.replace("\\", "/")
    p = find_ci(folder, name) if "/" not in name else folder / name
    return p if p and p.exists() else None


# --------------------------------------------------------------------------- #
# Components
# --------------------------------------------------------------------------- #
def read_basin(path: Path, folder: Path, warnings: list) -> dict:
    blocks = parse_blocks(read_text(path))
    hdr = (of_type(blocks, "Basin") or [{"items": [], "name": path.stem}])[0]
    elems = {}
    for b in blocks:
        if b["type"] not in ELEMENT_TYPES:
            continue
        e = {"name": b["name"], "type": b["type"], "downstream": bget(b, "Downstream"),
             "description": bget(b, "Description"), "x": num(bget(b, "Canvas X")), "y": num(bget(b, "Canvas Y")),
             "fx": num(bget(b, "From Canvas X")), "fy": num(bget(b, "From Canvas Y")), "params": bdict(b),
             "parts": [{"type": s["type"], "name": s["name"], **{k: v for k, v in s["items"]}} for s in b["sub"]]}
        if b["type"] == "Subbasin":
            e["area"] = num(bget(b, "Area"))
            e["methods"] = {m: bget(b, m) for m in SUBBASIN_METHODS if bget(b, m) not in (None, "", "None")}
        if b["type"] in ("Reach", "Reservoir"):
            e["route"] = bget(b, "Route")
        if b["type"] == "Junction" and bget(b, "Observed Hydrograph Gage"):
            e["observed_gage"] = bget(b, "Observed Hydrograph Gage")
        for k in ("Observed Hydrograph Gage", "Observed Pool Elevation Gage", "Observed Storage Gage"):
            if bget(b, k):
                e.setdefault("gages", []).append(bget(b, k))
        tabs = [v for k, v in b["items"] if k.endswith("Table Name") or k.endswith("Table") and not k.startswith("Use")]
        tabs += [v for s in b["sub"] for k, v in s["items"] if k.endswith("Table Name")]
        e["tables"] = [t for t in tabs if t]
        elems[b["name"]] = e
    # network checks
    for e in elems.values():
        d = e["downstream"]
        if d and d not in elems:
            warnings.append(f"basin {hdr['name']}: {e['type'].lower()} '{e['name']}' drains to unknown element '{d}'")
    upstream = {}
    for e in elems.values():
        if e["downstream"] in elems:
            upstream.setdefault(e["downstream"], []).append(e["name"])
    for e in elems.values():
        e["upstream"] = upstream.get(e["name"], [])
        if e["type"] == "Subbasin" and (e.get("area") is None or e["area"] <= 0):
            warnings.append(f"basin {hdr['name']}: subbasin '{e['name']}' has no positive area")
        if e["type"] in ("Reach", "Junction", "Reservoir", "Diversion", "Sink") and not e["upstream"]:
            warnings.append(f"basin {hdr['name']}: {e['type'].lower()} '{e['name']}' receives no inflow from any element")
    outlets = [e["name"] for e in elems.values() if not e["downstream"]]
    # map layers
    layers = []
    for b in of_type(blocks, "Basin Layer Properties"):
        for s in b["sub"]:
            fn = dict(s["items"]).get("Filename")
            if s["type"] == "Map Layer" and fn:
                p = _path(folder, fn)
                layers.append({"file": fn, "name": dict(s["items"]).get("Display name", Path(fn.replace("\\", "/")).stem),
                               "shown": dict(s["items"]).get("Layer shown"), "_path": str(p) if p else None})
    area = sum(e.get("area") or 0 for e in elems.values() if e["type"] == "Subbasin")
    counts = {t: sum(1 for e in elems.values() if e["type"] == t) for t in ELEMENT_TYPES}
    return {"name": hdr["name"], "file": path.name, "description": bget(hdr, "Description"),
            "settings": bdict(hdr, META_KEYS - {"Description"}), "unit_system": bget(hdr, "Unit System"),
            "elements": elems, "counts": {k: v for k, v in counts.items() if v}, "outlets": outlets,
            "total_area": area, "layers": layers,
            "computation_points": [b["name"] for b in of_type(blocks, "Computation Point")]}


def read_met(path: Path, warnings: list) -> dict:
    blocks = parse_blocks(read_text(path))
    hdr = (of_type(blocks, "Meteorology") or [{"items": [], "name": path.stem}])[0]
    methods = {k: v for k, v in hdr["items"] if k.endswith("Method") and v not in ("None", "")}
    subs = {}
    for b in of_type(blocks, "Subbasin"):
        pairs, cur = [], None
        for k, v in b["items"]:
            if k == "Gage":
                cur = [v, None, None]
                pairs.append(cur)
            elif k == "Volume Weight" and cur is not None:
                cur[1] = num(v)
            elif k == "Temporal Weight" and cur is not None:
                cur[2] = num(v)
        other = {k: v for k, v in b["items"] if k not in ("Gage", "Volume Weight", "Temporal Weight") and k not in META_KEYS}
        subs[b["name"]] = {"gages": pairs, "other": other}
        vw = [p[1] for p in pairs if p[1] is not None]
        if methods.get("Precipitation Method") == "Weighted Gages" and vw and not 0.95 <= sum(vw) <= 1.05:
            warnings.append(f"meteorologic model {hdr['name']}: volume weights of subbasin '{b['name']}' sum to {g(sum(vw))}")
    gages = [{"name": b["name"], **bdict(b)} for b in of_type(blocks, "Gage")]
    params = [{"name": b["name"], **bdict(b)} for b in blocks if b["type"].endswith("Method Parameters")]
    return {"name": hdr["name"], "file": path.name, "description": bget(hdr, "Description"),
            "unit_system": bget(hdr, "Unit System"), "methods": methods, "basins": bgetall(hdr, "Use Basin Model"),
            "subbasins": subs, "gages": gages, "method_parameters": params,
            "other_blocks": sorted({b["type"] for b in blocks} - {"Meteorology", "Subbasin", "Gage"} - {p["name"] for p in params})}


def read_control(path: Path) -> dict:
    blocks = parse_blocks(read_text(path))
    b = (of_type(blocks, "Control") or [{"items": [], "name": path.stem}])[0]
    t0 = _date(bget(b, "Start Date", ""), bget(b, "Start Time"))
    t1 = _date(bget(b, "End Date", ""), bget(b, "End Time"))
    dt = num(bget(b, "Time Interval"))
    out = {"name": b["name"], "file": path.name, "description": bget(b, "Description"),
           "start": t0.strftime("%Y-%m-%d %H:%M") if t0 else f"{bget(b, 'Start Date')} {bget(b, 'Start Time')}",
           "end": t1.strftime("%Y-%m-%d %H:%M") if t1 else f"{bget(b, 'End Date')} {bget(b, 'End Time')}",
           "interval_min": dt}
    if t0 and t1:
        hrs = (t1 - t0).total_seconds() / 3600
        out["duration_h"] = round(hrs, 3)
        if dt:
            out["steps"] = int(round(hrs * 60 / dt))
    return out


def read_results(path: Path) -> dict | None:
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError) as e:
        return {"file": path.name, "error": str(e)}
    res = {"file": path.name, "run": root.findtext("RunName"), "executed": root.findtext("ExecutionTime"),
           "start": root.findtext("StartTime"), "end": root.findtext("EndTime"), "elements": []}
    for el in root.findall("BasinElement"):
        stats = {}
        for s in el.iter("StatisticMeasure"):
            stats[s.get("type")] = (s.get("value"), s.get("units"))
        da = el.find("DrainageArea")
        res["elements"].append({"name": el.get("name"), "type": el.get("type"),
                                "drainage_area": (da.get("area"), da.get("units")) if da is not None else None,
                                "upstream": [f.get("upstreamName") for f in el.findall("FlowLink")],
                                "observed_gage": (el.find("ObservedFlowGage").get("name") if el.find("ObservedFlowGage") is not None else None),
                                "stats": stats, "series": len(el.findall(".//TimeSeries"))})
    return res


def read_log(path: Path) -> dict:
    text = read_text(path)
    kinds = {"NOTE": 0, "WARNING": 0, "ERROR": 0}
    seen, msgs = set(), []
    for ln in text.splitlines():
        m = re.match(r"^(NOTE|WARNING|ERROR)\s+(\d+):\s+(.*)$", ln.strip())
        if not m:
            continue
        kinds[m.group(1)] += 1
        if m.group(1) != "NOTE" and m.group(3) not in seen:
            seen.add(m.group(3))
            msgs.append(f"{m.group(1)} {m.group(2)}: {m.group(3)}")
    return {"file": path.name, "counts": kinds, "messages": msgs,
            "finished": "Finished computing" in text, "_text": text[-20000:]}


def read_manager(path: Path, block_type: str) -> list[dict]:
    if not path or not path.exists():
        return []
    out = []
    for b in of_type(parse_blocks(read_text(path)), block_type):
        d = {"name": b["name"], **bdict(b, META_KEYS - {"Description"})}
        d.pop(block_type, None)
        for s in b["sub"]:
            if s["type"] == "Variant":
                it = dict(s["items"])
                d["period"] = f"{it.get('Start Time', '')} .. {it.get('End Time', '')}".strip(" .")
        out.append(d)
    return out


def dss_catalog(folder: Path) -> list[dict]:
    out = []
    for p in sorted(folder.glob("*.dss")):
        d = {"file": p.name, "size_kb": round(p.stat().st_size / 1024, 1)}
        dsc = find_ci(folder, p.stem + ".dsc")
        if dsc:
            paths = re.findall(r"^\s*\d+\s+\S+\s+(/.*/)\s*$", read_text(dsc), re.M)
            d["records"] = len(paths)
            parts = {}
            for pn in paths:
                f = pn.split("/")
                if len(f) > 3:
                    parts.setdefault(f[3], set()).add(f[2])
            d["parameters"] = {k: len(v) for k, v in sorted(parts.items())}
        out.append(d)
    return out


# --------------------------------------------------------------------------- #
# Project
# --------------------------------------------------------------------------- #
def find_projects(arg: str) -> list[Path]:
    p = Path(arg).expanduser().resolve()
    if not p.exists():
        sys.exit(f"Path not found: {p}")
    if p.is_file():
        if p.suffix.lower() != ".hms":
            hms = sorted(p.parent.glob("*.hms"))
            if hms:
                return hms[:1]
        return [p]
    fs = sorted(p.rglob("*.hms"))
    if not fs:
        sys.exit(f"No .hms project file found under {p}")
    return fs


def build_project(path: Path) -> dict:
    folder = path.parent
    warnings = []
    blocks = parse_blocks(read_text(path))
    hdr = (of_type(blocks, "Project") or [{"items": [], "name": path.stem}])[0]

    def comp(btype, reader, key="Filename"):
        out = []
        for b in of_type(blocks, btype):
            fn = bget(b, key) or bget(b, "FileName") or bget(b, "Filename")
            p = _path(folder, fn)
            if not p:
                warnings.append(f"{btype.lower()} '{b['name']}': file {fn} not found")
                out.append({"name": b["name"], "file": fn, "missing": True, "description": bget(b, "Description")})
                continue
            try:
                d = reader(p)
            except Exception as e:  # keep going on a damaged file
                warnings.append(f"{btype.lower()} '{b['name']}' ({p.name}): {e}")
                d = {"name": b["name"], "file": p.name, "error": str(e)}
            d.setdefault("description", bget(b, "Description"))
            out.append(d)
        return out

    basins = comp("Basin", lambda p: read_basin(p, folder, warnings))
    mets = comp("Precipitation", lambda p: read_met(p, warnings))
    controls = comp("Control", read_control)

    runs = []
    run_file = find_ci(folder, path.stem + ".run")
    if run_file:
        for b in of_type(parse_blocks(read_text(run_file)), "Run"):
            r = {"name": b["name"], "description": bget(b, "Description"), "basin": bget(b, "Basin"),
                 "met": bget(b, "Precip"), "control": bget(b, "Control"), "dss": bget(b, "DSS File"),
                 "log": bget(b, "Log File"), "executed": " ".join(x for x in (bget(b, "Last Execution Date"), bget(b, "Last Execution Time")) if x),
                 "settings": {k: v for k, v in bdict(b).items() if k not in ("Basin", "Precip", "Control", "DSS File", "Log File",
                                                                           "Last Execution Date", "Last Execution Time")}}
            for kind, lst in (("basin", basins), ("met", mets), ("control", controls)):
                if r[kind] and r[kind] not in {x["name"] for x in lst}:
                    warnings.append(f"run '{r['name']}': {kind} '{r[kind]}' is not in the project")
            res = None
            if (folder / "results").is_dir():
                res = find_ci(folder / "results", f"RUN_{r['name'].replace(' ', '_')}.results") or \
                    find_ci(folder / "results", f"RUN_{r['name']}.results")
            if res:
                r["results"] = read_results(res)
            lg = _path(folder, r["log"]) if r["log"] else find_ci(folder, f"{r['name']}.log")
            if lg:
                r["log_summary"] = read_log(lg)
                if r["log_summary"]["counts"]["ERROR"]:
                    warnings.append(f"run '{r['name']}': {r['log_summary']['counts']['ERROR']} error(s) in {lg.name}")
            runs.append(r)

    gages = read_manager(find_ci(folder, path.stem + ".gage"), "Gage")
    tables = read_manager(find_ci(folder, path.stem + ".pdata"), "Table")
    grids = read_manager(find_ci(folder, path.stem + ".grid"), "Grid")
    gage_names = {x["name"] for x in gages}
    table_names = {x["name"] for x in tables}
    for bs in basins:
        for e in bs.get("elements", {}).values():
            for gn in e.get("gages", []):
                if gage_names and gn not in gage_names:
                    warnings.append(f"basin {bs['name']}: {e['name']} uses gage '{gn}' that is not in the gage manager")
            for tn in e.get("tables", []):
                if table_names and tn not in table_names:
                    warnings.append(f"basin {bs['name']}: {e['name']} uses paired-data table '{tn}' that is not in the project")
    for m in mets:
        for bn in m.get("basins", []):
            bs = next((b for b in basins if b["name"] == bn and "elements" in b), None)
            if bs and m.get("methods", {}).get("Precipitation Method") not in ("Gridded Precipitation", None):
                subs = {e["name"] for e in bs["elements"].values() if e["type"] == "Subbasin"}
                miss = sorted(subs - set(m["subbasins"]))
                if miss and m["subbasins"]:
                    warnings.append(f"meteorologic model {m['name']}: no data for subbasin(s) {', '.join(miss[:6])} of basin {bn}")

    project = {k: v for k, v in hdr["items"]}
    summary = {
        "Project": hdr["name"], "Description": project.get("Description"), "HMS version": project.get("Version"),
        "Time zone": project.get("Time Zone ID"), "Default DSS file": project.get("DSS File Name"),
        "Basin models": len(basins), "Meteorologic models": len(mets), "Control specifications": len(controls),
        "Simulation runs": len(runs), "Computed runs": sum(1 for r in runs if r.get("results")),
        "Time-series gages": len(gages), "Paired-data tables": len(tables), "Grids": len(grids),
    }
    files = sorted({q.name for q in folder.iterdir() if q.is_file()})
    return {"name": hdr["name"], "path": str(path), "folder": str(folder), "project": project, "summary": summary,
            "basins": basins, "mets": mets, "controls": controls, "runs": runs, "gages": gages, "tables": tables,
            "grids": grids, "dss": dss_catalog(folder), "files": files, "warnings": warnings}


def build_config(arg: str, only: str | None = None) -> dict:
    files = find_projects(arg)
    if only:
        files = [f for f in files if f.stem.lower() == only.lower() or f.parent.name.lower() == only.lower()] or files
    root = Path(arg).expanduser().resolve()
    return {"root": str(root if root.is_dir() else root.parent), "name": root.stem if root.is_file() else root.name,
            "projects": [build_project(f) for f in files]}


# --------------------------------------------------------------------------- #
# Terminal
# --------------------------------------------------------------------------- #
def _stat(stats, key):
    v = stats.get(key)
    if not v:
        return ""
    val, unit = v
    x = num(val)
    return (g(x) if x is not None else val) + (f" {unit}" if unit and x is not None else "")


def show_project(p, full=False):
    h1(f"HEC-HMS PROJECT: {p['name']}")
    row("File", p["path"])
    h2("Summary")
    for k, v in p["summary"].items():
        row(k, v)
    for b in p["basins"]:
        h2(f"Basin model: {b['name']} ({b.get('file')})")
        if b.get("missing") or b.get("error"):
            print(f"  ! {b.get('error') or 'file missing'}")
            continue
        row("Description", b.get("description"))
        row("Unit system", b.get("unit_system"))
        row("Elements", b["counts"])
        row("Total subbasin area", g(b["total_area"]))
        row("Outlets", b["outlets"])
        subs = [e for e in b["elements"].values() if e["type"] == "Subbasin"]
        if subs:
            print("  Subbasins")
            text_table(["Name", "Area", "Downstream", "Loss", "Transform", "Baseflow"],
                       [[e["name"], g(e.get("area")), e["downstream"], e["methods"].get("LossRate", ""),
                         e["methods"].get("Transform", ""), e["methods"].get("Baseflow", "")] for e in subs], indent=4,
                       limit=None if full else 20)
        others = [e for e in b["elements"].values() if e["type"] != "Subbasin"]
        if others:
            print("  Other elements")
            text_table(["Type", "Name", "Downstream", "Method / upstream"],
                       [[e["type"], e["name"], e["downstream"] or "(outlet)", e.get("route") or ", ".join(e["upstream"])]
                        for e in others], indent=4, limit=None if full else 25, maxw=40)
        if full:
            for e in b["elements"].values():
                print(f"    [{e['type']}] {e['name']}")
                for k, v in e["params"].items():
                    row(k, v, indent=8)
                for part in e["parts"]:
                    row(part["type"], {k: v for k, v in part.items() if k != "type"}, indent=8)
    for m in p["mets"]:
        h2(f"Meteorologic model: {m['name']} ({m.get('file')})")
        if m.get("missing") or m.get("error"):
            print(f"  ! {m.get('error') or 'file missing'}")
            continue
        row("Description", m.get("description"))
        for k, v in m["methods"].items():
            row(k, v)
        row("Basin models", m["basins"])
        for gg in m["gages"]:
            row(f"gage {gg['name']}", {k: v for k, v in gg.items() if k != "name"})
        for sn, s in m["subbasins"].items():
            if s["gages"]:
                row(sn, ", ".join(f"{gn} {g(w) if w is not None else ''}".strip() for gn, w, _ in s["gages"]), indent=4)
            elif s["other"] and full:
                row(sn, s["other"], indent=4)
    if p["controls"]:
        h2("Control specifications")
        text_table(["Name", "Start", "End", "Interval (min)", "Duration (h)", "Steps"],
                   [[c["name"], c.get("start", ""), c.get("end", ""), g(c.get("interval_min")), g(c.get("duration_h")), c.get("steps", "")]
                    for c in p["controls"]])
    if p["runs"]:
        h2("Simulation runs")
        text_table(["Run", "Basin", "Met", "Control", "Last computed", "Results"],
                   [[r["name"], r["basin"], r["met"], r["control"], r["executed"], "yes" if r.get("results") else "-"] for r in p["runs"]],
                   maxw=26)
        for r in p["runs"]:
            res = r.get("results")
            if res and res.get("elements"):
                print(f"\n  Run '{r['name']}' results ({res.get('start')} - {res.get('end')})")
                text_table(["Element", "Type", "Peak outflow", "Time of peak", "Volume"],
                           [[e["name"], e["type"], _stat(e["stats"], "Outflow Maximum"), _stat(e["stats"], "Outflow Maximum Time"),
                             _stat(e["stats"], "Outflow Volume")] for e in res["elements"]], indent=4, limit=None if full else 15, maxw=26)
                obs = [e for e in res["elements"] if "Observed Flow Nash Sutcliffe" in e["stats"]]
                for e in obs:
                    print(f"    observed at {e['name']}: peak {_stat(e['stats'], 'Observed Flow Maximum')}, "
                          f"NSE {_stat(e['stats'], 'Observed Flow Nash Sutcliffe')}, bias {_stat(e['stats'], 'Observed Flow Percent Bias')} %")
            lg = r.get("log_summary")
            if lg:
                c = lg["counts"]
                print(f"    log {lg['file']}: {c['NOTE']} notes, {c['WARNING']} warnings, {c['ERROR']} errors")
                for msg in lg["messages"][:8]:
                    print(f"      {msg}")
    if p["gages"]:
        h2(f"Time-series gages ({len(p['gages'])})")
        text_table(["Name", "Type", "Source", "File", "Pathname"],
                   [[x["name"], x.get("Gage Type", ""), x.get("Data Source Type", ""), x.get("Filename", ""), x.get("Pathname", "")]
                    for x in p["gages"]], maxw=46, limit=None if full else 20)
    if p["tables"]:
        h2(f"Paired-data tables ({len(p['tables'])})")
        text_table(["Name", "Type", "X units", "Y units", "Pathname"],
                   [[x["name"], x.get("Table Type", ""), x.get("X-Units", ""), x.get("Y-Units", ""), x.get("Pathname", "")]
                    for x in p["tables"]], maxw=46)
    if p["grids"]:
        h2(f"Grids ({len(p['grids'])})")
        text_table(["Name", "Type", "Source", "File", "Pathname"],
                   [[x["name"], x.get("Grid Type", ""), x.get("Data Source Type", ""), x.get("Filename", ""), x.get("Pathname", "")]
                    for x in p["grids"]], maxw=46)
    if p["dss"]:
        h2("HEC-DSS files (not read)")
        for d in p["dss"]:
            row(d["file"], f"{d['size_kb']} kB" + (f", {d['records']} records in catalog" if d.get("records") else ""))
    if p["warnings"]:
        h2("Warnings")
        for w in p["warnings"]:
            print(f"  ! {w}")


def show(cfg, full=False):
    for p in cfg["projects"]:
        show_project(p, full)
    print()


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #
EXTRA_CSS = """
svg.map .shp{fill:#e8f1e4;fill-opacity:.6;stroke:#9bbf8f;stroke-width:1} svg.map .shpl{fill:none;stroke:#90caf9;stroke-width:1.5}
svg.map .lnk{stroke:#7b8794;stroke-width:1.2;stroke-dasharray:4 3} svg.map .rch{fill:none;stroke:#1565c0;stroke-width:3}
svg.map .el{stroke:#fff;stroke-width:1}
"""


def basin_map(b, res=None) -> str:
    E = b["elements"]
    pts = [(e["x"], e["y"]) for e in E.values() if e["x"] is not None] + \
          [(e["fx"], e["fy"]) for e in E.values() if e["fx"] is not None]
    if len(pts) < 2:
        return '<p class="muted">The basin model holds no element coordinates - no schematic.</p>'
    ex0, ex1 = min(p[0] for p in pts), max(p[0] for p in pts)
    ey0, ey1 = min(p[1] for p in pts), max(p[1] for p in pts)
    shapes = []
    for lay in b.get("layers", []):
        if not lay.get("_path"):
            continue
        try:
            sh = read_shp(Path(lay["_path"]))
            kind = shp_kind(Path(lay["_path"]))
        except (OSError, ValueError, struct_error):
            continue
        xs = [x for s in sh for part in s for x, _ in part]
        ys = [y for s in sh for part in s for _, y in part]
        if kind != "point" and xs and not (max(xs) < ex0 or min(xs) > ex1 or max(ys) < ey0 or min(ys) > ey1):  # same CRS
            shapes.append((lay, sh, kind))
    allp = pts + [pt for _, sh, _ in shapes for s in sh for part in s for pt in part[::max(1, len(part) // 50)]]
    F = MapFrame(allp)
    o = []
    for lay, sh, kind in sorted(shapes, key=lambda t: t[2] != "polygon"):  # polygons below lines
        for s in sh:
            for part in s:
                if kind == "polygon":
                    o.append(f'<polygon points="{F.pts(part)}" class="shp"><title>{_e(lay["name"])}</title></polygon>')
                else:
                    o.append(f'<polyline points="{F.pts(part)}" class="shpl"><title>{_e(lay["name"])}</title></polyline>')
    peaks = {}
    if res and res.get("elements"):
        for re_ in res["elements"]:
            peaks[re_["name"]] = _stat(re_["stats"], "Outflow Maximum")
    for e in E.values():
        d = E.get(e["downstream"]) if e["downstream"] else None
        if e["type"] == "Reach" and e["fx"] is not None and e["x"] is not None:
            o.append(f'<polyline points="{F.pts([(e["fx"], e["fy"]), (e["x"], e["y"])])}" class="rch"><title>reach {_e(e["name"])}'
                     f'{" - peak " + peaks[e["name"]] if e["name"] in peaks else ""}</title></polyline>')
        elif d and e["x"] is not None and d["x"] is not None:
            tx, ty = (d["fx"], d["fy"]) if d["type"] == "Reach" and d["fx"] is not None else (d["x"], d["y"])
            o.append(f'<line x1="{F.X(e["x"]):.1f}" y1="{F.Y(e["y"]):.1f}" x2="{F.X(tx):.1f}" y2="{F.Y(ty):.1f}" class="lnk"/>')
    many = len(E) > 60
    for e in E.values():
        if e["x"] is None or e["type"] == "Reach":
            continue
        col, kind = ELEMENT_STYLE.get(e["type"], ("#555", "dot"))
        x, y = F.X(e["x"]), F.Y(e["y"])
        tip = f"{e['type'].lower()} {e['name']}" + (f", area {g(e.get('area'))}" if e.get("area") else "") + \
              (f" - peak outflow {peaks[e['name']]}" if e["name"] in peaks else "")
        if kind == "bx":
            o.append(f'<rect x="{x - 6:.1f}" y="{y - 6:.1f}" width="12" height="12" fill="{col}" class="el"><title>{_e(tip)}</title></rect>')
        else:
            o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5.5" fill="{col}" class="el"><title>{_e(tip)}</title></circle>')
        if not many:
            right = x > F.w - 8 * len(e["name"]) - 20  # keep labels near the right edge inside the map
            anchor = ' text-anchor="end"' if right else ""
            o.append(f'<text x="{x - 8 if right else x + 8:.1f}" y="{y - 6:.1f}"{anchor}>{_e(e["name"])}</text>')
    items = [(f"border-color:{ELEMENT_STYLE[t][0]};background:{ELEMENT_STYLE[t][0]}", t.lower(), ELEMENT_STYLE[t][1])
             for t in ELEMENT_TYPES if b["counts"].get(t)]
    items.append(("border-color:#7b8794;border-top-style:dashed", "drains to", "ln"))
    if shapes:
        items.append(("background:#e8f1e4;border:1px solid #9bbf8f", "map layer: " + ", ".join(lay["name"] for lay, _, _ in shapes), "bx"))
    return F.svg(o) + legend(items)


def _basin_html(b, runs) -> str:
    if b.get("missing") or b.get("error"):
        return f'<p class="warn">{_e(b.get("error") or "file " + str(b.get("file")) + " not found")}</p>'
    res = next((r.get("results") for r in runs if r["basin"] == b["name"] and r.get("results")), None)
    out = [f'<p class="muted">{_e(b.get("description") or "")} <code>{_e(b["file"])}</code></p>',
           '<div class="stats">' + "".join(f'<div class="stat"><b>{v}</b><span>{PLURAL.get(k, k.lower())}</span></div>' for k, v in b["counts"].items())
           + f'<div class="stat"><b>{g(b["total_area"])}</b><span>subbasin area</span></div></div>',
           basin_map(b, res)]
    subs = [e for e in b["elements"].values() if e["type"] == "Subbasin"]
    if subs:
        meths = sorted({m for e in subs for m in e["methods"]}, key=SUBBASIN_METHODS.index)
        out.append("<h4>Subbasins</h4>" + _table(["Name", "Area", "Downstream"] + meths,
                                                 [[_e(e["name"]), g(e.get("area")), _e(e["downstream"])] + [_e(e["methods"].get(m, "")) for m in meths]
                                                  for e in subs]))
    others = [e for e in b["elements"].values() if e["type"] != "Subbasin"]
    if others:
        out.append("<h4>Reaches, junctions and other elements</h4>" +
                   _table(["Type", "Name", "Downstream", "Method", "Upstream"],
                          [[e["type"], _e(e["name"]), _e(e["downstream"] or "(outlet)"), _e(e.get("route") or ""),
                            f'<span class="wrap">{_e(", ".join(e["upstream"]))}</span>'] for e in others]))
    det = []
    for e in b["elements"].values():
        pairs = list(e["params"].items())
        parts = "".join(f"<p><b>{_e(pt['type'])}</b> {_e(pt['name'])}</p>" + _kv([(k, v) for k, v in pt.items() if k not in ("type", "name")])
                        for pt in e["parts"])
        det.append(_blk(f"{_e(e['type'])} · {_e(e['name'])}", _kv(pairs) + parts, open_=False))
    out.append(_blk("All element parameters", "".join(det), f"{len(det)} elements", open_=False))
    if b["settings"]:
        out.append(_blk("Basin model options", _kv(list(b["settings"].items())), open_=False))
    return "".join(out)


def _met_html(m) -> str:
    if m.get("missing") or m.get("error"):
        return f'<p class="warn">{_e(m.get("error") or "file not found")}</p>'
    out = [_two([("File", m["file"]), ("Description", m.get("description")), ("Unit system", m.get("unit_system")),
                 ("Basin models", ", ".join(m["basins"]))] + list(m["methods"].items()))]
    if m["gages"]:
        cols = sorted({k for gg in m["gages"] for k in gg if k != "name"})
        out.append("<h4>Gages</h4>" + _table(["Gage"] + cols, [[_e(gg["name"])] + [_e(gg.get(c, "")) for c in cols] for gg in m["gages"]]))
    rows = [[_e(sn), _e(", ".join(f"{gn} ({g(w)})" if w is not None else gn for gn, w, _ in s["gages"])),
             f'<span class="wrap">{_e(", ".join(f"{k}: {v}" for k, v in s["other"].items()))}</span>']
            for sn, s in m["subbasins"].items()]
    if rows:
        out.append("<h4>Subbasin data</h4>" + _table(["Subbasin", "Gages (volume weight)", "Other settings"], rows))
    for pm in m["method_parameters"]:
        out.append(_blk(_e(pm["name"]), _kv([(k, v) for k, v in pm.items() if k != "name"]), "method parameters", open_=False))
    return "".join(out)


def _run_html(r) -> str:
    out = [_two([("Description", r.get("description")), ("Basin model", r["basin"]), ("Meteorologic model", r["met"]),
                 ("Control specification", r["control"]), ("Output DSS file", r["dss"]), ("Log file", r["log"]),
                 ("Last computed", r["executed"])] + list(r["settings"].items()))]
    res = r.get("results")
    if res and res.get("elements"):
        out.append(f'<h4>Results <small class="muted">({_e(res.get("start"))} - {_e(res.get("end"))}, computed {_e(res.get("executed"))})</small></h4>')
        out.append(_table(["Element", "Type", "Drainage area"] + [lab for _, lab in STAT_LABELS],
                          [[_e(e["name"]), _e(e["type"]), _e(" ".join(e["drainage_area"]) if e["drainage_area"] else "")]
                           + [_e(_stat(e["stats"], k)) for k, _ in STAT_LABELS] for e in res["elements"]]))
        obs = [e for e in res["elements"] if any(k in e["stats"] for k, _ in OBS_LABELS)]
        if obs:
            out.append("<h4>Observed flow comparison</h4>" +
                       _table(["Element", "Gage", "Computed peak"] + [lab for _, lab in OBS_LABELS],
                              [[_e(e["name"]), _e(e.get("observed_gage") or ""), _e(_stat(e["stats"], "Outflow Maximum"))]
                               + [_e(_stat(e["stats"], k)) for k, _ in OBS_LABELS] for e in obs]))
        pk = [(e["name"], num(e["stats"].get("Outflow Maximum", (None,))[0]), e["stats"].get("Outflow Maximum", (None, ""))[1])
              for e in res["elements"]]
        pk = [t for t in pk if t[1] is not None]
        if len(pk) > 1:
            pk.sort(key=lambda t: -t[1])
            out.append("<h4>Peak outflow by element</h4>" + svg_bars([t[0] for t in pk], [t[1] for t in pk], pk[0][2] or ""))
    elif res and res.get("error"):
        out.append(f'<p class="warn">results file {_e(res["file"])}: {_e(res["error"])}</p>')
    else:
        out.append('<p class="muted">No results file (results/RUN_&lt;name&gt;.results) - the run has not been computed in this folder.</p>')
    lg = r.get("log_summary")
    if lg:
        c = lg["counts"]
        body = (f'<p>{c["NOTE"]} notes, <span class="{"warn" if c["WARNING"] else ""}">{c["WARNING"]} warnings</span>, '
                f'<span class="{"warn" if c["ERROR"] else ""}">{c["ERROR"]} errors</span>'
                f'{" · finished" if lg["finished"] else ""}</p>')
        if lg["messages"]:
            body += '<ul class="warns">' + "".join(f"<li>{_e(x)}</li>" for x in lg["messages"][:40]) + "</ul>"
        body += f'<pre class="log">{_e(lg["_text"])}</pre>'
        out.append(_blk(f"Log · {_e(lg['file'])}", body, open_=bool(c["ERROR"] or c["WARNING"])))
    return "".join(out)


def _project_html(p, pi) -> tuple[str, list]:
    present, out = [], []

    def add(bid, title, body, note=""):
        if body:
            present.append((bid, title))
            out.append(_blk(title, body, note, bid=bid))

    add(f"p{pi}-sum", "Summary", _two(list(p["summary"].items())))
    for bi, b in enumerate(p["basins"]):
        add(f"p{pi}-b{bi}", f"Basin model · {b['name']}", _basin_html(b, p["runs"]), b.get("file") or "")
    for mi, m in enumerate(p["mets"]):
        add(f"p{pi}-m{mi}", f"Meteorologic model · {m['name']}", _met_html(m), m.get("file") or "")
    if p["controls"]:
        add(f"p{pi}-c", "Control specifications",
            _table(["Name", "Description", "Start", "End", "Interval (min)", "Duration (h)", "Steps", "File"],
                   [[_e(c["name"]), _e(c.get("description") or ""), _e(c.get("start", "")), _e(c.get("end", "")),
                     g(c.get("interval_min")), g(c.get("duration_h")), c.get("steps", ""), _e(c.get("file"))] for c in p["controls"]]))
    for ri, r in enumerate(p["runs"]):
        add(f"p{pi}-r{ri}", f"Run · {r['name']}", _run_html(r), "computed" if r.get("results") else "not computed")
    ts = []
    if p["gages"]:
        ts.append("<h4>Time-series gages</h4>" + _table(["Name", "Type", "Source", "File", "Pathname", "Period"],
                  [[_e(x["name"]), _e(x.get("Gage Type", "")), _e(x.get("Data Source Type", "")), _e(x.get("Filename", "")),
                    _e(x.get("Pathname", "")), _e(x.get("period", ""))] for x in p["gages"]]))
    if p["tables"]:
        ts.append("<h4>Paired-data tables</h4>" + _table(["Name", "Type", "X units", "Y units", "Description", "Pathname"],
                  [[_e(x["name"]), _e(x.get("Table Type", "")), _e(x.get("X-Units", "")), _e(x.get("Y-Units", "")),
                    _e(x.get("Description", "")), _e(x.get("Pathname", ""))] for x in p["tables"]]))
    if p["grids"]:
        ts.append("<h4>Grids</h4>" + _table(["Name", "Type", "Source", "File", "Pathname"],
                  [[_e(x["name"]), _e(x.get("Grid Type", "")), _e(x.get("Data Source Type", "")), _e(x.get("Filename", "")),
                    _e(x.get("Pathname", ""))] for x in p["grids"]]))
    if p["dss"]:
        ts.append("<h4>HEC-DSS files</h4><p class='muted'>DSS files are binary and are not read; the record count comes from the .dsc catalog when present.</p>"
                  + _table(["File", "Size (kB)", "Records", "Parameters (locations)"],
                           [[_e(d["file"]), g(d["size_kb"]), d.get("records", ""),
                             f'<span class="wrap">{_e(", ".join(f"{k} ({v})" for k, v in d.get("parameters", {}).items()))}</span>']
                            for d in p["dss"]]))
    add(f"p{pi}-ts", "Time series, paired data & grids", "".join(ts))
    add(f"p{pi}-f", "Files in the project folder", '<p class="wrap">' + " ".join(f"<code>{_e(f)}</code>" for f in p["files"]) + "</p>",
        f"{len(p['files'])} files", )
    return "".join(out), present


def render_html(cfg) -> str:
    nav, body = [], []
    for i, p in enumerate(cfg["projects"]):
        computed = sum(1 for r in p["runs"] if r.get("results"))
        badge = (f'<span class="badge">{computed}/{len(p["runs"])} runs computed</span>' if p["runs"] else '<span class="badge o">no runs</span>')
        if p["warnings"]:
            badge += f' <span class="badge w">{len(p["warnings"])} warning(s)</span>'
        sec = [f'<section id="p{i}"><h2>{_e(p["name"])} {badge}</h2><p class="muted"><code>{_e(p["path"])}</code></p>']
        if p["warnings"]:
            sec.append('<ul class="warns">' + "".join(f"<li>{_e(w)}</li>" for w in p["warnings"]) + "</ul>")
        h, present = _project_html(p, i)
        sec.append(h + "</section>")
        body.append("".join(sec))
        nav.append(f'<a href="#p{i}"><b>{_e(p["name"])}</b></a>')
        nav += [f'<a class="sub" href="#{b}">{_e(t)}</a>' for b, t in present]
    p0 = cfg["projects"][0]
    sub = f'{len(cfg["projects"])} projects' if len(cfg["projects"]) > 1 else \
        " · ".join(x for x in (p0["summary"].get("Description"), f'HMS {p0["summary"].get("HMS version")}' if p0["summary"].get("HMS version") else "") if x)
    return page("HEC-HMS", "hmsvue", cfg["name"], sub, nav, body, cfg["root"], EXTRA_CSS)


def write_html(cfg, out: Path, open_browser=True):
    return write_page(render_html(cfg), out, open_browser)


# --------------------------------------------------------------------------- #
def main(argv=None, prog=None):
    ap = argparse.ArgumentParser(prog=prog, description="Show the configuration of a HEC-HMS project.")
    ap.add_argument("model", help=".hms project file, or a folder containing projects")
    ap.add_argument("--project", dest="only", metavar="NAME", help="when several projects are found, show only this one")
    ap.add_argument("--all", action="store_true", help="terminal: every parameter of every element")
    ap.add_argument("--json", metavar="FILE", help="also write the configuration to a JSON file")
    ap.add_argument("--html", nargs="?", const="", metavar="FILE",
                    help="write the configuration as a webpage and open it in the browser, instead of printing it "
                         "(default file: <project>_config.html next to the .hms)")
    ap.add_argument("--no-browser", action="store_true", help="with --html, write the page but don't open it")
    a = ap.parse_args(argv)

    cfg = build_config(a.model, a.only)
    if a.html is None:
        show(cfg, a.all)
    else:
        p0 = cfg["projects"][0]
        out = Path(a.html) if a.html else (Path(p0["folder"]) / f"{Path(p0['path']).stem}_config.html" if len(cfg["projects"]) == 1
                                           else Path(cfg["root"]) / f"{cfg['name']}_config.html")
        write_html(cfg, out, open_browser=not a.no_browser)
    if a.json:
        Path(a.json).write_text(json.dumps(strip_private(cfg), indent=2, default=str), encoding="utf-8")
        print(f"Configuration written to {a.json}")


if __name__ == "__main__":
    main()
