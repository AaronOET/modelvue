"""
modelvue.dfmvue - Read a Delft3D FM (D-Flow FM / D-HYDRO) 1D2D model and print
its configuration.

Parses the MDU and the files it points to - UGRID network/mesh file (*_net.nc),
new-style INI inputs (external forcing .ext, boundary .bc, cross sections,
roughness, storage nodes, structures, initial fields), polylines (.pli) - and,
if the model has been run, the diagnostics (.dia) and map / his NetCDF output.

No extra packages needed: classic NetCDF-3 files (what D-Flow FM / Delta Shell
write by default) are read with a built-in reader. NetCDF-4 files need netCDF4
or h5py. numpy is optional (faster).

Usage
-----
    dfmvue  path/to/FM_model.mdu
    dfmvue  path/to/project_folder          # finds the .mdu (incl. *.dsproj_data)
    dfmvue  path/to/Project.dsproj
    dfmvue  folder --model FM_model         # pick one MDU when several exist
    dfmvue  model.mdu --all                 # terminal: every MDU keyword
    dfmvue  model.mdu --json config.json    # also dump to JSON
    dfmvue  model.mdu --html                # webpage only, opens in browser
    dfmvue  model.mdu --html report.html --no-browser
"""
from __future__ import annotations

import argparse
import datetime as _dt
import html as _h
import json
import re
import sys
from pathlib import Path


from ._common import _blk, _e, h1, h2, _kv, row, strip_private, _table, _two, wrap_pairs
from ._dflowfm import (
    bdict, bget, bgetall, BRANCH_GUI_TYPE, decode, FILE_KEYS, find_ci, find_mdus, fmt_dur, g,
    MATERIAL, MDU, mdu_groups, parse_ext_old, parse_ini, parse_initial_fields, parse_pli, ramp,
    RAMP_BED, RAMP_DEPTH, read_network, ref_datetime, results, sections, svg_lines, text_table,
    _time_unit, to_float, tunit_sec,
)


# --------------------------------------------------------------------------- #
# Forcing files
# --------------------------------------------------------------------------- #


def parse_bc(path: Path) -> list[dict]:
    out = []
    for b in parse_ini(path):
        if b["section"].lower() != "forcing":
            continue
        qs, us = bgetall(b, "quantity"), bgetall(b, "unit")
        rows = []
        for ln in b["data"]:
            vals = [to_float(t) for t in ln.split()]
            if vals and all(v is not None for v in vals):
                rows.append(vals)
        f = {"name": bget(b, "name"), "function": bget(b, "function"),
             "interpolation": bget(b, "timeInterpolation"),
             "quantities": [f"{q} [{u}]" if u else q for q, u in zip(qs, us + [""] * len(qs))],
             "rows": len(rows)}
        fn = (f["function"] or "").lower()
        if fn == "timeseries" and rows and qs and qs[0].lower() == "time":
            fac, ref = _time_unit(us[0] if us else "")
            t = [r[0] * fac for r in rows]
            f["reference"] = ref
            f["time_span"] = f"{g(t[0])} .. {g(t[-1])} s ({(t[-1] - t[0]) / 3600:.2f} h)"
            f["_series"] = []
            for c in range(1, len(rows[0])):
                vals = [r[c] for r in rows if len(r) > c]
                f.setdefault("value_range", []).append(f"{g(min(vals))} .. {g(max(vals))}")
                f["_series"].append((qs[c] if c < len(qs) else f"col{c + 1}", [x / 3600 for x in t], vals))
            f["value_range"] = "; ".join(f["value_range"])
        elif rows:
            flat = [v for r in rows for v in r]
            f["value_range"] = f"{g(min(flat))} .. {g(max(flat))}"
        out.append(f)
    return out


def parse_ext_new(path: Path, folder: Path, warnings: list) -> dict:
    out = {"file": path.name, "boundaries": [], "laterals": [], "meteo": [], "sources": []}
    for b in parse_ini(path):
        sec = b["section"].lower()
        d = bdict(b)
        if sec == "boundary":
            e = {"quantity": bget(b, "quantity"), "location": bget(b, "locationFile") or bget(b, "nodeId"),
                 "forcing_file": ";".join(bgetall(b, "forcingFile")), "return_time": bget(b, "returnTime")}
            if bget(b, "locationFile"):
                lf = find_ci(folder, bget(b, "locationFile"))
                if lf is None:
                    warnings.append(f"boundary locationFile '{bget(b, 'locationFile')}' not found")
                else:
                    try:
                        e["polylines"] = [{k: v for k, v in p.items() if k != "_xy"} for p in parse_pli(lf)]
                        e["_pli"] = parse_pli(lf)
                    except OSError as ex:
                        warnings.append(str(ex))
            e["forcings"] = []
            for ff in bgetall(b, "forcingFile"):
                fp = find_ci(folder, ff)
                if fp is None:
                    warnings.append(f"boundary forcingFile '{ff}' not found")
                    continue
                try:
                    bc = parse_bc(fp)
                except OSError as ex:
                    warnings.append(str(ex))
                    continue
                pref = [x["name"] for x in e.get("polylines", [])]
                e["forcings"] += [f for f in bc if not pref or any((f["name"] or "").startswith(p) for p in pref)
                                  or f["name"] == e["location"]] or bc
            out["boundaries"].append(e)
        elif sec == "lateral":
            e = {"id": d.get("id"), "name": d.get("name"), "type": d.get("type"),
                 "location": d.get("branchId") and f"{d.get('branchId')} @ {d.get('chainage')}" or
                             (d.get("nodeId") or (f"{d.get('numCoordinates')} coordinates" if d.get("numCoordinates") else None)),
                 "locationType": d.get("locationType"), "discharge": d.get("discharge")}
            fp = find_ci(folder, d.get("discharge", "")) if d.get("discharge") and to_float(d.get("discharge")) is None else None
            if fp and fp.suffix.lower() == ".bc":
                e["forcings"] = [f for f in parse_bc(fp) if f["name"] == d.get("id")] or parse_bc(fp)[:1]
            out["laterals"].append(e)
        elif sec == "meteo":
            out["meteo"].append({"quantity": d.get("quantity"), "forcingFile": d.get("forcingFile"),
                                 "forcingFileType": d.get("forcingFileType"), "method": d.get("interpolationMethod"),
                                 "operand": d.get("operand"),
                                 "exists": find_ci(folder, d.get("forcingFile", "")) is not None})
        elif sec == "sourcesink":
            out["sources"].append(d)
    return out


# --------------------------------------------------------------------------- #
# 1D inputs
# --------------------------------------------------------------------------- #
def parse_1d(mdu: MDU, warnings: list) -> dict:
    d = mdu.dir
    out = {}

    def load(key, sect="geometry"):
        res = []
        for fn in [x.strip() for x in (mdu(sect, key) or "").split(";") if x.strip()]:
            fp = find_ci(d, fn)
            if fp is None:
                warnings.append(f"{key} '{fn}' not found")
                continue
            try:
                res.append((fp, parse_ini(fp)))
            except OSError as e:
                warnings.append(str(e))
        return res

    # branch properties (Delta Shell branches.gui)
    gui = {}
    for fp, bl in load("BranchFile"):
        for b in sections(bl, "Branch"):
            gui[bget(b, "name")] = {"type": BRANCH_GUI_TYPE.get(bget(b, "branchType", ""), bget(b, "branchType")),
                                    "material": MATERIAL.get(bget(b, "material", ""), bget(b, "material")),
                                    "source": bget(b, "sourceCompartmentName"),
                                    "target": bget(b, "targetCompartmentName")}
    out["branch_gui"] = gui

    # cross sections
    defs = {}
    for fp, bl in load("CrossDefFile"):
        for b in sections(bl, "Definition"):
            t = (bget(b, "type") or "").lower()
            dims = ""
            if t == "circle":
                dims = f"D = {bget(b, 'diameter')}"
            elif t == "rectangle":
                dims = f"W = {bget(b, 'width')}, H = {bget(b, 'height')}, closed = {bget(b, 'closed', '-')}"
            elif t in ("yz", "xyz"):
                zs = [to_float(x) for x in (bget(b, "zCoordinates") or "").split()]
                ys = [to_float(x) for x in (bget(b, "yCoordinates") or "").split()]
                if zs and ys and None not in zs + ys:
                    dims = f"{len(zs)} pts, width {g(max(ys) - min(ys))}, z {g(min(zs))} .. {g(max(zs))}"
            elif t in ("zw", "tabulated", "zwriver"):
                lv = [to_float(x) for x in (bget(b, "levels") or "").split()]
                fw = [to_float(x) for x in (bget(b, "flowWidths") or bget(b, "totalWidths") or "").split()]
                if lv and None not in lv:
                    dims = f"{len(lv)} levels {g(min(lv))} .. {g(max(lv))}" + (f", max width {g(max(fw))}" if fw and None not in fw else "")
            defs[bget(b, "id")] = {"id": bget(b, "id"), "type": bget(b, "type"), "dims": dims,
                                   "friction": bget(b, "frictionIds") or bget(b, "frictionId") or "",
                                   "thalweg": bget(b, "thalweg")}
    locs = []
    for fp, bl in load("CrossLocFile"):
        for b in sections(bl, "CrossSection"):
            locs.append({"id": bget(b, "id"), "branch": bget(b, "branchId"),
                         "chainage": to_float(bget(b, "chainage")), "shift": bget(b, "shift"),
                         "definition": bget(b, "definitionId")})
    for L in locs:
        if L["definition"] and L["definition"] not in defs:
            warnings.append(f"cross section {L['id']}: definition '{L['definition']}' not in CrossDefFile")
    out["cross_section_definitions"] = list(defs.values())
    out["cross_section_locations"] = locs

    # roughness
    rough = []
    for fp, bl in load("FrictFile"):
        glob = sections(bl, "Global")
        br = sections(bl, "Branch")
        rough.append({"file": fp.name,
                      "global": [f"{bget(x, 'frictionId')}: {bget(x, 'frictionType')} {bget(x, 'frictionValue')}" for x in glob],
                      "branch_specific": len(br),
                      "branches": sorted({bget(x, "branchId") for x in br if bget(x, "branchId")})})
    out["roughness"] = rough

    # storage nodes / manholes
    nodes = []
    for fp, bl in load("StorageNodeFile"):
        hdr = sections(bl, "General")
        if hdr and bget(hdr[0], "useStreetStorage"):
            out["use_street_storage"] = bget(hdr[0], "useStreetStorage")
        for b in sections(bl, "StorageNode"):
            nodes.append({"id": bget(b, "id"), "node": bget(b, "nodeId"), "manhole": bget(b, "ManholeId"),
                          "type": bget(b, "storageType"), "bedLevel": bget(b, "bedLevel"), "area": bget(b, "area"),
                          "streetLevel": bget(b, "streetLevel"), "streetStorageArea": bget(b, "streetStorageArea"),
                          "useTable": bget(b, "useTable")})
    out["storage_nodes"] = nodes

    # structures
    structs = []
    for fp, bl in load("StructureFile"):
        for b in sections(bl, "Structure"):
            dct = bdict(b)
            loc = (f"{dct.get('branchId')} @ {dct.get('chainage')}" if dct.get("branchId") else
                   dct.get("polylinefile") or (f"{dct.get('numCoordinates')} coords" if dct.get("numCoordinates") else ""))
            main = {k: v for k, v in dct.items()
                    if k.lower() not in ("id", "name", "type", "branchid", "chainage", "polylinefile")}
            structs.append({"id": dct.get("id"), "name": dct.get("name"), "type": dct.get("type"),
                            "location": loc, "params": main})
    out["structures"] = structs
    return out


# --------------------------------------------------------------------------- #
# Model summary / MDU groups
# --------------------------------------------------------------------------- #


def time_info(mdu: MDU) -> dict:
    f = tunit_sec(mdu)
    t0 = to_float(mdu("time", "TStart"), 0) * f
    t1 = to_float(mdu("time", "TStop"), 0) * f
    ref = ref_datetime(mdu)
    out = {"RefDate": mdu("time", "RefDate"), "Tunit": mdu("time", "Tunit", "S"), "Tzone": mdu("time", "Tzone"),
           "TStart": fmt_dur(t0), "TStop": fmt_dur(t1), "duration": fmt_dur(t1 - t0),
           "DtUser": mdu("time", "DtUser"), "DtMax": mdu("time", "DtMax"), "DtInit": mdu("time", "DtInit"),
           "_t0": t0, "_t1": t1}
    if ref:
        out["period"] = f"{ref + _dt.timedelta(seconds=t0):%Y-%m-%d %H:%M} -> {ref + _dt.timedelta(seconds=t1):%Y-%m-%d %H:%M}"
    return out


def referenced_files(mdu: MDU) -> list[dict]:
    out = []
    for b in mdu.blocks:
        for k, v, _ in b["items"]:
            if FILE_KEYS.search(k) and v and k.lower() not in ("hisfile", "mapfile", "flowgeomfile", "classmapfile",
                                                               "restartfile" if not v else ""):
                for f in [x.strip() for x in v.split(";") if x.strip() and to_float(x) is None]:
                    out.append({"section": b["section"], "key": k, "file": f,
                                "exists": find_ci(mdu.dir, f) is not None})
    return out


def model_summary(mdu: MDU, net: dict, oned: dict, ext: dict, ti: dict) -> dict:
    br = net.get("branches", [])
    gui = oned.get("branch_gui", {})
    kinds = {}
    for b in br:
        k = gui.get(b["id"], {}).get("type") or ("Pipe" if "pipe" in b["id"].lower() else "Channel")
        kinds.setdefault(k, [0, 0.0])
        kinds[k][0] += 1
        kinds[k][1] += b["length_m"] or 0
    m2 = net.get("mesh2d", {})
    nb = sum(len(x.get("polylines", [])) or 1 for x in ext.get("boundaries", [])) if ext else 0
    out = {
        "Program": re.sub(r"^Deltares,\s*", "", mdu("General", "Version") or mdu("General", "Program") or ""),
        "Delta Shell GUI": mdu("General", "GuiVersion"),
        "MDU file version": mdu("General", "fileVersion"),
        "Simulation period": ti.get("period") or f"{ti['TStart']} -> {ti['TStop']}",
        "Duration": ti["duration"],
        "Time step": f"{decode(mdu, 'numerics', 'TimeStepType', '2')}; "
                     f"CFLMax {mdu('numerics', 'CFLMax', '0.7')}, DtMax {mdu('time', 'DtMax')} s, DtUser {mdu('time', 'DtUser')} s",
        "1D network": ", ".join(f"{n} {k.lower()}{'s' if n > 1 else ''} ({v / 1000:.2f} km)" for k, (n, v) in kinds.items()) or None,
        "1D mesh": f"{net['mesh1d']['nodes']} nodes, {net['mesh1d']['edges']} edges" if net.get("mesh1d") else None,
        "2D mesh": (f"{m2.get('faces')} cells ({m2.get('cell_types')}), {m2.get('approx_cell_size')}, {m2.get('total_area')}"
                    if m2 else None),
        "1D2D links": (f"{net['links_1d2d']['count']} (" + ", ".join(f"{k}: {v}" for k, v in net['links_1d2d']['by_type'].items()) + ")"
                       if net.get("links_1d2d") else None),
        "Projection": net.get("projection"),
        "Bed level type": decode(mdu, "geometry", "BedlevType", "3"),
        "Uniform friction": f"{decode(mdu, 'physics', 'UnifFrictType', '1')}, {mdu('physics', 'UnifFrictCoef')} (2D), 1D {mdu('physics', 'UnifFrictCoef1D', '-')}",
        "Open boundaries": f"{nb} ({', '.join(sorted({x['quantity'] for x in ext['boundaries'] if x['quantity']}))})" if ext and ext["boundaries"] else None,
        "Laterals": len(ext["laterals"]) if ext and ext["laterals"] else None,
        "Meteo forcing": ", ".join(m["quantity"] for m in ext["meteo"]) if ext and ext["meteo"] else None,
        "Storage nodes": len(oned.get("storage_nodes", [])) or None,
        "Structures": len(oned.get("structures", [])) or None,
        "Cross sections": f"{len(oned.get('cross_section_locations', []))} locations, {len(oned.get('cross_section_definitions', []))} definitions"
        if oned.get("cross_section_locations") else None,
        "Initial water level": mdu("geometry", "WaterLevIni"),
        "Infiltration": decode(mdu, "grw", "Infiltrationmodel", "0"),
        "Layers (Kmx)": mdu("geometry", "Kmx", "0") + (" (depth-averaged)" if mdu("geometry", "Kmx", "0") == "0" else ""),
        "Output": f"his every {fmt_dur(mdu('output', 'HisInterval')) or '-'}, map every {fmt_dur((mdu('output', 'MapInterval') or '').split()[0] if mdu('output', 'MapInterval') else None) or '-'}"
                  + (f", restart every {fmt_dur(mdu('output', 'RstInterval'))}" if to_float(mdu("output", "RstInterval"), 0) else ""),
    }
    return {k: v for k, v in out.items() if v not in (None, "", 0)}


# --------------------------------------------------------------------------- #
# Model discovery
# --------------------------------------------------------------------------- #


def build_model(path: Path) -> dict:
    warnings = []
    mdu = MDU(path)
    ti = time_info(mdu)
    m = {"name": path.stem, "mdu": str(path), "folder": str(path.parent)}
    netf = mdu("geometry", "NetFile")
    fp = find_ci(mdu.dir, netf or "")
    if not netf:
        warnings.append("no NetFile in [geometry]")
        net = {}
    elif fp is None:
        warnings.append(f"NetFile '{netf}' not found")
        net = {}
    else:
        net = read_network(fp, warnings)
    oned = parse_1d(mdu, warnings)
    ext = None
    nf = mdu("external forcing", "ExtForceFileNew")
    if nf:
        f = find_ci(mdu.dir, nf)
        if f is None:
            warnings.append(f"ExtForceFileNew '{nf}' not found")
        else:
            ext = parse_ext_new(f, f.parent, warnings)
    of = mdu("external forcing", "ExtForceFile")
    old = None
    if of:
        f = find_ci(mdu.dir, of)
        if f is None:
            warnings.append(f"ExtForceFile '{of}' not found")
        else:
            old = parse_ext_old(f, f.parent)
    m["time"] = ti
    m["network"] = net
    m["one_d"] = oned
    m["external_forcing"] = ext
    m["external_forcing_old"] = old
    m["initial_fields"] = parse_initial_fields(mdu, warnings)
    m["files"] = referenced_files(mdu)
    for f in m["files"]:
        if not f["exists"] and not any(f["file"] in w for w in warnings):
            warnings.append(f"{f['key']} '{f['file']}' not found")
    m["summary"] = model_summary(mdu, net, oned, ext, ti)
    m["mdu_parameters"] = mdu_groups(mdu)
    m["results"] = results(mdu, ti, warnings)
    # consistency checks
    gui_ids = set(oned.get("branch_gui", {}))
    net_ids = {b["id"] for b in net.get("branches", [])}
    if gui_ids - net_ids:
        warnings.append(f"branches.gui lists branches not in the net file: {', '.join(sorted(gui_ids - net_ids))}")
    crs_br = {L["branch"] for L in oned.get("cross_section_locations", [])}
    no_crs = sorted(net_ids - crs_br)
    if net_ids and no_crs:
        warnings.append(f"branches without a cross section (uniform profile used): {', '.join(no_crs)}")
    if ext:
        for b in ext["boundaries"]:
            if not b.get("forcings"):
                warnings.append(f"boundary {b['location']}: no matching forcing found")
            for fo in b.get("forcings", []):
                ts = fo.get("_series")
                if ts and ts[0][1] and ts[0][1][-1] * 3600 < ti["_t1"]:
                    warnings.append(f"forcing {fo['name']} ends before TStop")
    m["warnings"] = warnings
    return m


def build_config(arg: str, only=None) -> dict:
    mdus = find_mdus(arg)
    sel = [p for p in mdus if not only or p.stem.lower() == only.lower()]
    if not sel:
        sys.exit(f"Model '{only}' not found. MDUs: {', '.join(p.stem for p in mdus)}")
    root = Path(arg).expanduser().resolve()
    root = root if root.is_dir() else root.parent
    proj = next(iter(sorted(root.glob("*.dsproj"))), None) or next(
        (p.with_suffix("") for p in [sel[0]] for q in p.parents if q.name.endswith(".dsproj_data")
         for p in [q.with_name(q.name[:-5])] if p.exists()), None)
    return {"name": root.name if len(sel) > 1 or not proj else Path(str(proj)).stem,
            "root": str(root), "dsproj": str(proj) if proj else None,
            "models": [build_model(p) for p in sel]}


# --------------------------------------------------------------------------- #
# Terminal output
# --------------------------------------------------------------------------- #


def branch_rows(m):
    gui = m["one_d"].get("branch_gui", {})
    crs = {}
    for L in m["one_d"].get("cross_section_locations", []):
        crs[L["branch"]] = crs.get(L["branch"], 0) + 1
    return [[b["id"], gui.get(b["id"], {}).get("type") or ("Pipe" if "pipe" in b["id"].lower() else "Channel"), g(b["length_m"]), b["from"], b["to"],
             gui.get(b["id"], {}).get("material") or "", crs.get(b["id"], 0)]
            for b in m["network"].get("branches", [])]


def show_model(m, all_params=False):
    h1(f"D-FLOW FM 1D2D MODEL: {m['name']}")
    row("MDU file", m["mdu"])
    h2("Model summary")
    for k, v in m["summary"].items():
        row(k, v)
    net = m["network"]
    h2(f"Network & mesh ({net.get('file', '-')})")
    if net.get("error"):
        print(f"  ERROR: {net['error']}")
    row("Projection", net.get("projection"))
    for k, v in (net.get("mesh2d") or {}).items():
        row(f"2D {k}", v)
    for k, v in (net.get("mesh1d") or {}).items():
        row(f"1D mesh {k}", v)
    if net.get("links_1d2d"):
        row("1D2D links", f"{net['links_1d2d']['count']}  {net['links_1d2d']['by_type']}")
    if net.get("branches"):
        print("\n  1D branches")
        text_table(["Branch", "Type", "Length m", "From", "To", "Material", "Crs"], branch_rows(m), 4)
    o = m["one_d"]
    if o.get("cross_section_locations"):
        h2("Cross sections")
        text_table(["Location", "Branch", "Chainage", "Shift", "Definition"],
                   [[L["id"], L["branch"], g(L["chainage"]), L["shift"], L["definition"]] for L in o["cross_section_locations"]])
        print()
        text_table(["Definition", "Type", "Dimensions", "Friction"],
                   [[d["id"], d["type"], d["dims"], d["friction"]] for d in o["cross_section_definitions"]])
    if o.get("storage_nodes"):
        h2(f"Storage nodes / manholes{' (street storage on)' if o.get('use_street_storage') == '1' else ''}")
        text_table(["Id", "Manhole", "Type", "Bed", "Area", "Street lvl", "Street area"],
                   [[n["id"], n["manhole"], n["type"], n["bedLevel"], n["area"], n["streetLevel"], n["streetStorageArea"]]
                    for n in o["storage_nodes"]])
    if o.get("roughness"):
        h2("Roughness")
        for r in o["roughness"]:
            row(r["file"], "; ".join(r["global"]) + (f"  (+{r['branch_specific']} branch-specific)" if r["branch_specific"] else ""))
    h2("Structures")
    if o.get("structures"):
        for s in o["structures"]:
            print(f"  {s['type'] or '?':<12} {s['id']:<24} @ {s['location']}")
            wrap_pairs(list(s["params"].items())[:12], indent=6)
    else:
        print("  none")
    ext = m["external_forcing"]
    h2("Boundaries & forcing")
    if ext:
        for b in ext["boundaries"]:
            pl = ", ".join(f"{p['name']} ({p['points']} pts)" for p in b.get("polylines", []))
            print(f"  BND {b['quantity']:<18} @ {b['location']}  {pl}  <- {b['forcing_file']}")
            for f in b.get("forcings", []):
                print(f"      {f['name']}: {f['function']}, {f['rows']} rows, {f.get('time_span', '')}, value {f.get('value_range', '')}")
        for L in ext["laterals"]:
            print(f"  LAT {L['id']} ({L['type'] or ''}) @ {L['location']}  discharge = {L['discharge']}")
        for mt in ext["meteo"]:
            print(f"  MET {mt['quantity']} <- {mt['forcingFile']} ({mt['forcingFileType']}){'' if mt['exists'] else '  MISSING'}")
    for e in m["external_forcing_old"] or []:
        print(f"  OLD {e.get('QUANTITY')} <- {e.get('FILENAME')} ({e.get('FILETYPE')}, method {e.get('METHOD')})")
    if not ext and not m["external_forcing_old"]:
        print("  no external forcing file")
    if m["initial_fields"]:
        h2("Initial fields & parameters")
        for f in m["initial_fields"]:
            print(f"  [{f['block']}] {f['quantity']:<22} <- {f['dataFile']} ({f['dataFileType']}"
                  f"{', ' + f['method'] if f['method'] else ''}{', ' + f['locationType'] if f['locationType'] else ''})")
            if f.get("content"):
                print(f"      {f['content']}")
    t = m["time"]
    h2("Time")
    for k in ("RefDate", "Tunit", "TStart", "TStop", "period", "DtUser", "DtMax", "DtInit"):
        row(k, t.get(k))
    h2("MDU parameters" + ("" if all_params else " (non-empty; Wri* output switches summarised; --all for every key)"))
    for sec, items in m["mdu_parameters"].items():
        it = [e for e in items if all_params or e["value"] != ""]
        if not all_params:
            on = [e["key"] for e in it if re.match(r"wri(his|map|shp)_", e["key"], re.I) and e["value"] == "1"]
            it = [e for e in it if not re.match(r"wri(his|map|shp)_", e["key"], re.I)]
        if not it and not (not all_params and on):
            continue
        print(f"  [{sec}]")
        wrap_pairs([(e["key"], e["value"]) for e in it])
        if not all_params and on:
            print(f"    switched on: {', '.join(on)}")
    h2("Referenced files")
    for f in m["files"]:
        print(f"  {f['key']:<20} {f['file']:<40} {'ok' if f['exists'] else 'MISSING'}")
    h2("Results")
    r = m["results"]
    if not r:
        print("  No output folder found (model not run yet).")
    else:
        row("Output folder", r["output_dir"])
        d = r.get("dia")
        if d:
            row("Diagnostics", f"{d['file']} - {'finished' if d['finished'] else 'NOT finished'}, "
                               f"{len(d['errors'])} error(s), {len(d['warnings'])} warning(s)")
            for k in ("Computation started at", "Computation finished at", "total computation time (s)", "my model area (m2)"):
                row(k, d.get(k))
            row("Flow nodes / links", ", ".join(f"{k.replace('nr of ', '')}: {v}" for k, v in d["counts"].items()
                                                if "flow" in k or "1D" in k or "2D" in k))
            for x in d["errors"]:
                print(f"    ERROR   {x}")
            for x in d["warnings"]:
                print(f"    warning {x}")
        for kind in ("map", "his"):
            e = r.get(kind)
            if e:
                print(f"  {kind.upper()} {e['file']} ({e['size_MB']} MB)")
                for k, v in e.items():
                    if k not in ("file", "size_MB", "variables"):
                        row(k, v, indent=4)
                row("variables", f"{len(e.get('variables', []))}: {', '.join(e.get('variables', [])[:12])}"
                    + (" ..." if len(e.get("variables", [])) > 12 else ""), indent=4)
    if m["warnings"]:
        h2("Warnings")
        for w in m["warnings"]:
            print(f"  ! {w}")


def show(cfg, all_params=False):
    if len(cfg["models"]) > 1:
        h1(f"FOLDER: {cfg['root']}  -  {len(cfg['models'])} MDU files")
        for m in cfg["models"]:
            print(f"  {m['name']:<24} {m['summary'].get('Simulation period', '')}   {m['mdu']}")
    for m in cfg["models"]:
        show_model(m, all_params)
    print()


# --------------------------------------------------------------------------- #
# SVG graphics
# --------------------------------------------------------------------------- #


def svg_planview(m, values=None, ramp_stops=RAMP_BED, label="bed level (m)", w=880) -> str:
    """Plan view: 2D cells (coloured by values), 1D branches, nodes, 1D2D links, boundaries."""
    P = m["network"].get("_plot", {})
    faces = P.get("faces") or []
    xs, ys = [], []
    for px, py in faces:
        xs += px
        ys += py
    for _, bx, by in P.get("branches", []):
        xs += bx
        ys += by
    bnd = [pl for b in (m["external_forcing"] or {}).get("boundaries", []) for pl in b.get("_pli", [])]
    for pl in bnd:
        xs += [p[0] for p in pl["_xy"]]
        ys += [p[1] for p in pl["_xy"]]
    if not xs:
        return ""
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    span = max(x1 - x0, y1 - y0) or 1
    x0 -= span * 0.02
    x1 += span * 0.02
    y0 -= span * 0.02
    y1 += span * 0.02
    sc = (w - 20) / (x1 - x0)
    hgt = (y1 - y0) * sc + 20
    if hgt > 900:
        sc = 880 / (y1 - y0)
        hgt = 900
    X = lambda v: 10 + (v - x0) * sc
    Y = lambda v: 10 + (y1 - v) * sc
    o = [f'<svg class="plan" viewBox="0 0 {w} {hgt:.0f}">']
    vv = [v for v in (values or []) if v is not None]
    lo, hi = (min(vv), max(vv)) if vv else (0, 1)
    tips = len(faces) <= 4000
    for i, (px, py) in enumerate(faces[:60000]):
        v = values[i] if values and i < len(values) else None
        fill = ramp((v - lo) / ((hi - lo) or 1), ramp_stops) if v is not None else "none"
        pts = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in zip(px, py))
        t = f"<title>cell {i + 1}: {g(v)}</title>" if tips and v is not None else ""
        o.append(f'<polygon points="{pts}" fill="{fill}" class="cell">{t}</polygon>')
    # 1D2D links
    L = P.get("links")
    m1 = P.get("mesh1d_xy")
    fxy = P.get("face_xy")
    if L and m1 and fxy:
        s = P.get("links_start", 0)
        for i in range(0, len(L), 2):
            a, b = L[i] - s, L[i + 1] - s
            if 0 <= a < len(m1) and 0 <= b < len(fxy):
                o.append(f'<line x1="{X(m1[a][0]):.1f}" y1="{Y(m1[a][1]):.1f}" x2="{X(fxy[b][0]):.1f}" '
                         f'y2="{Y(fxy[b][1]):.1f}" class="lnk"/>')
    gui = m["one_d"].get("branch_gui", {})
    for name, bx, by in P.get("branches", []):
        kind = (gui.get(name, {}).get("type") or ("Pipe" if "pipe" in name.lower() else "Channel")).lower()
        cls = "pipe" if "pipe" in kind or "sewer" in kind else "chan"
        pts = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in zip(bx, by))
        o.append(f'<polyline points="{pts}" class="{cls}"><title>{_h.escape(name)} ({kind})</title></polyline>')
    if m1:
        for x, y in m1:
            o.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="1.6" class="m1"/>')
    stor = {n["node"] for n in m["one_d"].get("storage_nodes", [])} | {n["id"] for n in m["one_d"].get("storage_nodes", [])}
    for x, y, nid in P.get("net_nodes", []):
        cls = "stor" if nid in stor else "nn"
        o.append(f'<rect x="{X(x) - 4:.1f}" y="{Y(y) - 4:.1f}" width="8" height="8" class="{cls}"><title>{_h.escape(nid)}</title></rect>')
    for pl in bnd:
        pts = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b, _ in pl["_xy"])
        o.append(f'<polyline points="{pts}" class="bnd"><title>boundary {_h.escape(pl["name"])}</title></polyline>')
    o.append("</svg>")
    leg = (f'<div class="legend"><span><i class="sw" style="background:linear-gradient(90deg,'
           f'{",".join(ramp(t, ramp_stops) for t in (0, .25, .5, .75, 1))})"></i>{_h.escape(label)}: {g(lo)} .. {g(hi)}</span>'
           if vv else '<div class="legend"><span>2D cells (no values)</span>')
    leg += ('<span><i class="ln chan"></i>channel</span><span><i class="ln pipe"></i>pipe / sewer</span>'
            '<span><i class="ln lnk"></i>1D2D link</span><span><i class="ln bnd"></i>open boundary</span>'
            '<span><i class="bx stor"></i>storage node</span><span><i class="bx nn"></i>network node</span></div>')
    return "".join(o) + leg


# --------------------------------------------------------------------------- #
# HTML report (layout matches ras_config.py / hms_config.py)
# --------------------------------------------------------------------------- #
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<script>(function(){var t;try{t=localStorage.getItem('modelvue-theme')}catch(e){}if(t!=='dark')t='light';document.documentElement.dataset.theme=t;})();</script>
<title>D-Flow FM 1D2D · @@TITLE@@</title>
<style>
:root{--bg:#f6f7f9;--panel:#fff;--ink:#1d2433;--muted:#667085;--line:#e4e7ec;--acc:#3b4a9c;--accs:#eceffa;--warn:#b42318}
:root[data-theme=dark]{color-scheme:dark;--bg:#0f1318;--panel:#171c23;--ink:#e6e9ef;--muted:#98a2b3;--line:#2a313b;--acc:#9aa8f0;--accs:#1f2440;--warn:#f97066}
*{box-sizing:border-box} body{margin:0;font:14px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,"Noto Sans TC",sans-serif;background:var(--bg);color:var(--ink)}
header{position:sticky;top:0;z-index:5;background:var(--acc);color:#fff;padding:14px 24px;display:flex;gap:16px;align-items:baseline;flex-wrap:wrap}
header h1{margin:0;font-size:20px} header span{opacity:.85}
[data-theme=dark] header{background:#262f66}
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
<header><h1>D-Flow FM 1D2D · @@TITLE@@</h1><span>@@SUB@@</span><button class="theme" type="button" title="Switch between the light and dark theme">☾ Dark</button></header>
<div class="wrap"><nav><div class="tools"><button data-x="1">Expand all</button><button data-x="0">Collapse all</button></div>@@NAV@@</nav>
<main>@@BODY@@
<footer>Generated by modelvue.dfmvue on @@GEN@@ from <code>@@PATH@@</code></footer></main></div>
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


BLOCKS = [("summary", "Model summary"), ("network", "Network & mesh"), ("onedee", "1D: branches & cross sections"),
          ("nodes", "1D: storage nodes"), ("rough", "Roughness"), ("struct", "Structures"),
          ("forcing", "Boundaries & forcing"), ("init", "Initial fields"), ("time", "Time"),
          ("mdu", "MDU parameters"), ("files", "Referenced files"), ("results", "Results")]


def _model_section(m, mi) -> tuple[str, list]:
    pid = lambda k: f"m{mi}-{k}"
    present = []
    out = [f'<section id="m{mi}"><h2>{_e(m["name"])}']
    net, o, ext = m["network"], m["one_d"], m["external_forcing"]
    if net.get("branches"):
        out.append(' <span class="badge o">1D</span>')
    if net.get("mesh2d"):
        out.append(' <span class="badge o">2D</span>')
    if net.get("links_1d2d"):
        out.append(' <span class="badge o">1D2D</span>')
    r = m["results"]
    if r and r.get("dia"):
        ok = r["dia"]["finished"] and not r["dia"]["errors"]
        out.append(f' <span class="badge{"" if ok else " w"}">{"run finished" if ok else "run problem"}</span>')
    if m["warnings"]:
        out.append(f' <span class="badge w">{len(m["warnings"])} warning(s)</span>')
    out.append(f'</h2><p class="muted"><code>{_e(m["mdu"])}</code></p>')
    if m["warnings"]:
        out.append('<ul class="warns">' + "".join(f"<li>{_e(w)}</li>" for w in m["warnings"]) + "</ul>")

    def add(key, title, body, note="", open_=True):
        if body:
            present.append((pid(key), title))
            out.append(_blk(_e(title), body, note, open_, pid(key)))

    add("summary", "Model summary", _two(list(m["summary"].items())))

    # network
    P = net.get("_plot", {})
    vals, lab = None, "bed level (m)"
    if r and r.get("_plot", {}).get("bl_2d"):
        vals, lab = r["_plot"]["bl_2d"], "bed level used by model (m)"
    elif P.get("face_z") and any(v is not None for v in P["face_z"]):
        vals = P["face_z"]
    pairs = [("file", f"{net.get('file')} ({net.get('size_kB')} kB)"), ("projection", net.get("projection")),
             ("conventions", net.get("conventions"))]
    pairs += [(f"2D {k.replace('_', ' ')}", v) for k, v in (net.get("mesh2d") or {}).items() if k != "name"]
    pairs += [(f"1D mesh {k.replace('_', ' ')}", v) for k, v in (net.get("mesh1d") or {}).items() if k != "name"]
    if net.get("links_1d2d"):
        pairs.append(("1D2D links", f"{net['links_1d2d']['count']} - {_e(net['links_1d2d']['by_type'])}"))
    if net.get("legacy_net"):
        pairs.append(("legacy net", net["legacy_net"]))
    body = _two(pairs) + svg_planview(m, vals, RAMP_BED, lab)
    add("network", "Network & mesh", body, net.get("file", ""))

    # 1D
    if net.get("branches"):
        body = _table(["Branch", "Type", "Length (m)", "From", "To", "Material", "Cross sections"],
                      [[_e(c) for c in rr] for rr in branch_rows(m)])
        if o.get("cross_section_locations"):
            body += "<h4>Cross-section locations</h4>" + _table(
                ["Location", "Branch", "Chainage", "Shift", "Definition"],
                [[_e(L["id"]), _e(L["branch"]), _e(g(L["chainage"])), _e(L["shift"]), _e(L["definition"])]
                 for L in o["cross_section_locations"]])
            body += "<h4>Cross-section definitions</h4>" + _table(
                ["Definition", "Type", "Dimensions", "Friction", "Thalweg"],
                [[_e(d["id"]), _e(d["type"]), _e(d["dims"]), _e(d["friction"]), _e(d["thalweg"])]
                 for d in o["cross_section_definitions"]])
        add("onedee", "1D: branches & cross sections", body, f'{len(net["branches"])} branches')
    if o.get("storage_nodes"):
        add("nodes", "1D: storage nodes", _table(
            ["Id", "Network node", "Manhole", "Type", "Bed level", "Area", "Street level", "Street storage area", "Table"],
            [[_e(n[k]) for k in ("id", "node", "manhole", "type", "bedLevel", "area", "streetLevel", "streetStorageArea", "useTable")]
             for n in o["storage_nodes"]]), f'{len(o["storage_nodes"])} nodes'
            + (" · street storage on" if o.get("use_street_storage") == "1" else ""))
    if o.get("roughness"):
        add("rough", "Roughness", _table(["File", "Global", "Branch-specific"],
                                         [[_e(x["file"]), _e("; ".join(x["global"])),
                                           _e(f'{x["branch_specific"]} ({", ".join(x["branches"])})' if x["branch_specific"] else "-")]
                                          for x in o["roughness"]]))
    if o.get("structures"):
        add("struct", "Structures", _table(["Type", "Id", "Location", "Parameters"],
                                           [[_e(s["type"]), _e(s["id"]), _e(s["location"]),
                                             f'<span class="wrap">{"".join(f"<span class=chip>{_e(k)}={_e(v)}</span>" for k, v in s["params"].items())}</span>']
                                            for s in o["structures"]]), f'{len(o["structures"])}')
    elif m["files"] and any(f["key"].lower() == "structurefile" for f in m["files"]):
        add("struct", "Structures", '<p class="muted">Structure file present but it defines no structures.</p>', open_=False)

    # forcing
    body = ""
    if ext:
        rows, series = [], []
        for b in ext["boundaries"]:
            for f in b.get("forcings") or [{}]:
                rows.append([_e(b["quantity"]), _e(b["location"]), _e(f.get("name")), _e(f.get("function")),
                             _e(f.get("rows")), _e(f.get("time_span")), _e(f.get("value_range")), _e(b["forcing_file"])])
                for lab_, xs, ys in f.get("_series", []):
                    series.append((f"{f['name']}", xs, ys))
        body += _table(["Quantity", "Location", "Forcing", "Function", "Rows", "Time span", "Values", "File"], rows)
        if series:
            body += svg_lines(series[:8], "time since RefDate (h)", "boundary value")
        if ext["laterals"]:
            body += "<h4>Laterals</h4>" + _table(["Id", "Type", "Location", "Location type", "Discharge"],
                                                 [[_e(L["id"]), _e(L["type"]), _e(L["location"]), _e(L["locationType"]), _e(L["discharge"])]
                                                  for L in ext["laterals"]])
        if ext["meteo"]:
            body += "<h4>Meteo</h4>" + _table(["Quantity", "File", "File type", "Method", "Operand", "Exists"],
                                              [[_e(x["quantity"]), _e(x["forcingFile"]), _e(x["forcingFileType"]), _e(x["method"]),
                                                _e(x["operand"]), "yes" if x["exists"] else '<span class="warn">MISSING</span>']
                                               for x in ext["meteo"]])
    if m["external_forcing_old"]:
        body += "<h4>Old-style external forcing (ExtForceFile)</h4>" + _table(
            ["Quantity", "File", "Type", "Method", "Operand", "Exists"],
            [[_e(e.get("QUANTITY")), _e(e.get("FILENAME")), _e(e.get("FILETYPE")), _e(e.get("METHOD")),
              _e(e.get("OPERAND")), "yes" if e.get("exists") else '<span class="warn">MISSING</span>'] for e in m["external_forcing_old"]])
    add("forcing", "Boundaries & forcing", body or '<p class="muted">No external forcing file.</p>',
        ext["file"] if ext else "")

    if m["initial_fields"]:
        add("init", "Initial fields", _table(["Block", "Quantity", "Data file", "Type", "Method", "Location", "Content"],
                                             [[_e(f["block"]), _e(f["quantity"]), _e(f["dataFile"]), _e(f["dataFileType"]),
                                               _e(f["method"]), _e(f["locationType"]), (f'<span class="wrap">{_e(f.get("content"))}</span>')]
                                              for f in m["initial_fields"]]))
    t = m["time"]
    add("time", "Time", _two([(k, t.get(k)) for k in ("RefDate", "Tzone", "Tunit", "TStart", "TStop", "duration",
                                                         "period", "DtUser", "DtMax", "DtInit")]))

    # MDU parameters
    secs = []
    for sec, items in m["mdu_parameters"].items():
        nonempty = sum(1 for e in items if e["value"] != "")
        secs.append(_blk(_e(sec), _table(["Key", "Value", "Description"],
                                         [[f"<b>{_e(e['key'])}</b>", _e(e["value"]), (_e(e["description"]), "desc")]
                                          for e in items]), f"{len(items)} keys, {nonempty} set",
                         open_=sec.lower() in ("general", "geometry", "time", "numerics", "physics")))
    add("mdu", "MDU parameters", '<input class="pfilter" placeholder="filter keys / values / descriptions…" '
                                 'style="width:100%;max-width:420px;padding:6px 10px;border:1px solid var(--line);'
                                 'border-radius:6px;background:var(--panel);color:var(--ink);margin:4px 0 6px">' + "".join(secs),
        f'{sum(len(v) for v in m["mdu_parameters"].values())} keys', open_=False)

    add("files", "Referenced files", _table(["Section", "Key", "File", "Status"],
                                            [[_e(f["section"]), _e(f["key"]), f'<code>{_e(f["file"])}</code>',
                                              "ok" if f["exists"] else '<span class="warn">MISSING</span>'] for f in m["files"]]),
        f'{len(m["files"])} file(s)')

    # results
    if not r:
        body = '<p class="muted">No output folder found (model not run yet).</p>'
    else:
        body = _kv([("Output folder", r["output_dir"])])
        d = r.get("dia")
        if d:
            body += "<h4>Diagnostics (.dia)</h4>" + _two(
                [("file", d["file"]), ("finished", "yes" if d["finished"] else "NO"),
                 ("errors", len(d["errors"])), ("warnings", len(d["warnings"]))]
                + [(k, v) for k, v in d.items() if k not in ("file", "finished", "errors", "warnings", "counts")]
                + list(d["counts"].items()))
            if d["errors"] or d["warnings"]:
                body += _table(["Level", "Message"], [['<span class="warn">ERROR</span>', (_e(x), "msg")] for x in d["errors"]]
                               + [["warning", (_e(x), "msg")] for x in d["warnings"]], sortable=False)
        for kind in ("map", "his"):
            e = r.get(kind)
            if e:
                body += f"<h4>{kind} file</h4>" + _kv([(k, v) for k, v in e.items() if k != "variables"])
                if e.get("variables"):
                    body += f'<p>{"".join(f"<span class=chip>{_e(v)}</span>" for v in e["variables"])}</p>'
        md = r.get("_plot", {}).get("max_depth_2d")
        if md and net.get("_plot", {}).get("faces"):
            body += "<h4>Maximum 2D water depth over the run</h4>" + svg_planview(m, md, RAMP_DEPTH, "max water depth (m)")
    add("results", "Results", body)
    out.append("</section>")
    return "".join(out), present


def render_html(cfg: dict) -> str:
    nav, body = [], []
    for i, m in enumerate(cfg["models"]):
        sec, present = _model_section(m, i)
        body.append(sec)
        nav.append(f'<a href="#m{i}"><b>{_e(m["name"])}</b></a>')
        nav += [f'<a class="sub" href="#{bid}">{_e(t)}</a>' for bid, t in present]
    m0 = cfg["models"][0]
    sub = (f'{len(cfg["models"])} models' if len(cfg["models"]) > 1 else
           " · ".join(x for x in (m0["summary"].get("Simulation period"), m0["summary"].get("2D mesh", "").split(",")[0]) if x))
    rep = {"@@TITLE@@": _e(cfg["name"]), "@@SUB@@": _e(sub), "@@NAV@@": "".join(nav), "@@BODY@@": "".join(body),
           "@@GEN@@": _dt.datetime.now().strftime("%Y-%m-%d %H:%M"), "@@PATH@@": _e(cfg["root"])}
    out = HTML_TEMPLATE
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
    ap = argparse.ArgumentParser(prog=prog, description="Show the configuration of a Delft3D FM (D-Flow FM) 1D2D model.")
    ap.add_argument("model", help=".mdu file, .dsproj file, or a folder containing the model")
    ap.add_argument("--model", dest="only", metavar="NAME", help="when several .mdu files are found, show only this one")
    ap.add_argument("--all", action="store_true", help="terminal: list every MDU keyword, including empty ones")
    ap.add_argument("--json", metavar="FILE", help="also write the configuration to a JSON file")
    ap.add_argument("--html", nargs="?", const="", metavar="FILE",
                    help="write the configuration as a webpage and open it in the browser, instead of printing it "
                         "(default file: <name>_config.html next to the .mdu)")
    ap.add_argument("--no-browser", action="store_true", help="with --html, write the page but don't open it")
    a = ap.parse_args(argv)

    cfg = build_config(a.model, a.only)
    if a.html is None:
        show(cfg, a.all)
    else:
        m0 = cfg["models"][0]
        out = Path(a.html) if a.html else Path(m0["folder"]) / f"{m0['name'] if len(cfg['models']) == 1 else cfg['name']}_config.html"
        write_html(cfg, out, open_browser=not a.no_browser)
    if a.json:
        Path(a.json).write_text(json.dumps(strip_private(cfg), indent=2, default=str), encoding="utf-8")
        print(f"Configuration written to {a.json}")


if __name__ == "__main__":
    main()
