"""
modelvue.dfmvue23 - Read a Delft3D FM (D-Flow FM / D-HYDRO / Delft3D FM Suite) 2D or 3D
model and print its configuration.

Parses the MDU and the files it points to - UGRID mesh (*_net.nc), external
forcing (.ext new and old style), boundary conditions (.bc: time series,
harmonic, astronomic, 3D profiles), polylines (boundaries, thin dams, fixed
weirs, cross sections), observation points, structures, initial fields,
sediment / morphology (.sed / .mor), and a coupled D-Waves (SWAN) model (.mdw
with its curvilinear grids, bed levels, boundaries and obstacles). If the model
has been run it also summarises the .dia, map and his files.

Option codes (e.g. Turbulencemodel = 3) are decoded from the MDU's own comments,
so the labels match the D-Flow FM version that wrote the file.

No extra packages needed: classic NetCDF-3 files are read with a built-in
reader and the mesh picture is drawn as a PNG with the standard library.
NetCDF-4 files need netCDF4. numpy is optional (faster).

Usage
-----
    dfmvue23  path/to/FlowFM.mdu
    dfmvue23  path/to/project_folder          # finds the .mdu (incl. *.dsproj_data)
    dfmvue23  path/to/Project.dsproj
    dfmvue23  folder --model FlowFM           # pick one MDU when several exist
    dfmvue23  model.mdu --all                 # terminal: every MDU keyword
    dfmvue23  model.mdu --json config.json    # also dump to JSON
    dfmvue23  model.mdu --html                # webpage only, opens in browser
    dfmvue23  model.mdu --html report.html --no-browser
"""
from __future__ import annotations

import argparse
import base64
import datetime as _dt
import html as _h
import json
import math
import re
import struct
import sys
import zlib
from pathlib import Path


from ._common import _blk, _e, h1, h2, _kv, row, strip_private, _table, _two, wrap_pairs
from ._dflowfm import (
    bdict, bget, bgetall, decode, FILE_KEYS, find_ci, find_mdus, fmt_dur, g, MDU, mdu_groups,
    parse_ext_old, parse_ini, parse_initial_fields, parse_pli, ramp, RAMP_BED, RAMP_DEPTH,
    read_network, read_text, ref_datetime, results, sections, svg_lines, _time_unit, to_float,
    tunit_sec, vstats,
)


# --------------------------------------------------------------------------- #
# Forcing files
# --------------------------------------------------------------------------- #


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
# Model summary / MDU groups
# --------------------------------------------------------------------------- #


def referenced_files(mdu: MDU) -> list[dict]:
    out = []
    for b in mdu.blocks:
        for k, v, _ in b["items"]:
            if FILE_KEYS.search(k) and v and k.lower() not in ("hisfile", "mapfile", "flowgeomfile", "classmapfile",
                                                               "restartfile" if not v else ""):
                for f in [x for x in re.split(r"[;\s]+", v) if x and to_float(x) is None]:
                    out.append({"section": b["section"], "key": k, "file": f,
                                "exists": find_ci(mdu.dir, f) is not None})
    return out


# --------------------------------------------------------------------------- #
# HTML report (layout matches ras_config.py / hms_config.py)
# --------------------------------------------------------------------------- #
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<script>(function(){var t;try{t=localStorage.getItem('modelvue-theme')}catch(e){}if(t!=='dark')t='light';document.documentElement.dataset.theme=t;})();</script>
<title>D-Flow FM 2D/3D · @@TITLE@@</title>
<style>
:root{--bg:#f6f7f9;--panel:#fff;--ink:#1d2433;--muted:#667085;--line:#e4e7ec;--acc:#0b6e74;--accs:#e2f2f2;--warn:#b42318}
:root[data-theme=dark]{color-scheme:dark;--bg:#0f1318;--panel:#171c23;--ink:#e6e9ef;--muted:#98a2b3;--line:#2a313b;--acc:#7fd0d4;--accs:#163234;--warn:#f97066}
*{box-sizing:border-box} body{margin:0;font:14px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,"Noto Sans TC",sans-serif;background:var(--bg);color:var(--ink)}
header{position:sticky;top:0;z-index:5;background:var(--acc);color:#fff;padding:14px 24px;display:flex;gap:16px;align-items:baseline;flex-wrap:wrap}
header h1{margin:0;font-size:20px} header span{opacity:.85}
[data-theme=dark] header{background:#0f3e42}
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
<header><h1>D-Flow FM 2D/3D · @@TITLE@@</h1><span>@@SUB@@</span><button class="theme" type="button" title="Switch between the light and dark theme">☾ Dark</button></header>
<div class="wrap"><nav><div class="tools"><button data-x="1">Expand all</button><button data-x="0">Collapse all</button></div>@@NAV@@</nav>
<main>@@BODY@@
<footer>Generated by modelvue.dfmvue23 on @@GEN@@ from <code>@@PATH@@</code></footer></main></div>
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
# Boundary-condition files (.bc) incl. harmonic / astronomic / 3D profiles
# --------------------------------------------------------------------------- #
def parse_bc(path: Path) -> list[dict]:
    out = []
    for b in parse_ini(path):
        if b["section"].lower() != "forcing":
            continue
        qs, us = bgetall(b, "quantity"), bgetall(b, "unit")
        fn = (bget(b, "function") or "").lower()
        rows = [ln.split() for ln in b["data"] if ln.split()]
        f = {"name": bget(b, "name"), "function": bget(b, "function"),
             "interpolation": bget(b, "timeInterpolation"),
             "quantities": [f"{q} [{u}]" if u else q for q, u in zip(qs, us + [""] * len(qs))],
             "rows": len(rows)}
        if bget(b, "vertPositions"):
            f["vertical"] = f"{bget(b, 'vertPositionType', '')}: {bget(b, 'vertPositions')}"
        if fn.startswith("harmonic") or fn.startswith("astronomic"):
            comps = []
            for r in rows:
                if len(r) >= 3:
                    lab = r[0] if fn.startswith("astronomic") else f"{g(to_float(r[0]))} {us[0] if us else 'min'}"
                    comps.append({"component": lab, "amplitude": to_float(r[1]), "phase": to_float(r[2])})
            f["components"] = comps
            f["value_range"] = "; ".join(f"{c['component']}: A={g(c['amplitude'])}, phase={g(c['phase'])}"
                                         for c in comps[:5]) + (" ..." if len(comps) > 5 else "")
        elif fn.startswith("timeseries") or fn == "t3d":
            num = [[to_float(t) for t in r] for r in rows]
            num = [r for r in num if r and None not in r]
            if num and qs and qs[0].lower() == "time":
                fac, ref = _time_unit(us[0] if us else "")
                t = [r[0] * fac for r in num]
                f["reference"] = ref
                f["time_span"] = f"{g(t[0])} .. {g(t[-1])} s ({(t[-1] - t[0]) / 3600:.2f} h)"
                f["_t_end"] = t[-1]
                rng, f["_series"] = [], []
                for c in range(1, len(num[0])):
                    vals = [r[c] for r in num if len(r) > c]
                    rng.append(f"{g(min(vals))} .. {g(max(vals))}")
                    if c <= 4:
                        f["_series"].append((qs[c] if c < len(qs) else f"col{c + 1}", [x / 3600 for x in t], vals))
                f["value_range"] = "; ".join(dict.fromkeys(rng))
        elif fn == "constant" and rows:
            f["value_range"] = " ".join(rows[0])
        elif rows:
            flat = [to_float(t) for r in rows for t in r if to_float(t) is not None]
            if flat:
                f["value_range"] = f"{g(min(flat))} .. {g(max(flat))}"
        out.append(f)
    return out


def pli_length(pts):
    return sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(pts, pts[1:]))


def split_files(v):
    return [x for x in re.split(r"[;\s]+", v or "") if x]


def read_polyline_files(mdu, key, section="geometry", warnings=None) -> dict | None:
    files = split_files(mdu(section, key))
    if not files:
        return None
    out = {"key": key, "files": files, "polylines": [], "missing": []}
    for fn in files:
        fp = find_ci(mdu.dir, fn)
        if fp is None:
            out["missing"].append(fn)
            if warnings is not None:
                warnings.append(f"{key} '{fn}' not found")
            continue
        try:
            for p in parse_pli(fp):
                p["file"] = fp.name
                p["length_m"] = pli_length(p["_xy"])
                out["polylines"].append(p)
        except OSError as e:
            if warnings is not None:
                warnings.append(str(e))
    out["count"] = len(out["polylines"])
    out["total_length_m"] = sum(p["length_m"] for p in out["polylines"])
    return out


def parse_pliz(fp: Path) -> list[dict]:
    """Fixed weirs: x y crest [left right ...]"""
    res = []
    for p in parse_pli(fp):
        res.append(p)
    lines = [ln.split() for ln in read_text(fp).splitlines() if ln.strip() and not ln.strip().startswith("*")]
    crest = [to_float(r[2]) for r in lines if len(r) >= 3 and to_float(r[0]) is not None and to_float(r[2]) is not None]
    for p in res:
        p["_crest"] = crest
    return res


def read_obs(mdu, warnings) -> list[dict]:
    out = []
    for fn in split_files(mdu("output", "ObsFile")):
        fp = find_ci(mdu.dir, fn)
        if fp is None:
            warnings.append(f"ObsFile '{fn}' not found")
            continue
        if fp.suffix.lower() == ".ini":
            for b in parse_ini(fp):
                if b["section"].lower() == "observationpoint":
                    out.append({"name": bget(b, "name"), "x": to_float(bget(b, "x")), "y": to_float(bget(b, "y")),
                                "branch": bget(b, "branchId"), "chainage": bget(b, "chainage"), "file": fp.name})
        else:
            for ln in read_text(fp).splitlines():
                t = ln.split(None, 2)
                if len(t) >= 2 and to_float(t[0]) is not None:
                    out.append({"name": (t[2].strip().strip("'\"") if len(t) > 2 else ""),
                                "x": to_float(t[0]), "y": to_float(t[1]), "file": fp.name})
    return out


def read_crs(mdu, warnings) -> list[dict]:
    out = []
    for fn in split_files(mdu("output", "CrsFile")):
        fp = find_ci(mdu.dir, fn)
        if fp is None:
            warnings.append(f"CrsFile '{fn}' not found")
            continue
        if fp.suffix.lower() == ".ini":
            for b in parse_ini(fp):
                if b["section"].lower() == "observationcrosssection":
                    xs = [to_float(x) for x in (bget(b, "xCoordinates") or "").split()]
                    ys = [to_float(x) for x in (bget(b, "yCoordinates") or "").split()]
                    out.append({"name": bget(b, "name"), "_xy": list(zip(xs, ys, [""] * len(xs))), "file": fp.name})
        else:
            for p in parse_pli(fp):
                p["file"] = fp.name
                out.append(p)
    return out


def read_structures(mdu, warnings) -> list[dict]:
    out = []
    for fn in split_files(mdu("geometry", "StructureFile")):
        fp = find_ci(mdu.dir, fn)
        if fp is None:
            warnings.append(f"StructureFile '{fn}' not found")
            continue
        for b in parse_ini(fp):
            if b["section"].lower() != "structure":
                continue
            d = bdict(b)
            loc = d.get("polylinefile") or (f"{d.get('branchId')} @ {d.get('chainage')}" if d.get("branchId") else
                                            f"{d.get('numCoordinates', '')} coords" if d.get("numCoordinates") else "")
            s = {"id": d.get("id"), "type": d.get("type"), "location": loc,
                 "params": {k: v for k, v in d.items() if k.lower() not in ("id", "name", "type", "polylinefile",
                                                                              "branchid", "chainage")}}
            if d.get("polylinefile"):
                pf = find_ci(fp.parent, d["polylinefile"])
                if pf:
                    try:
                        s["_xy"] = [pt for p in parse_pli(pf) for pt in p["_xy"]]
                    except OSError:
                        pass
                else:
                    warnings.append(f"structure {d.get('id')}: polylinefile '{d['polylinefile']}' not found")
            elif d.get("xCoordinates"):
                xs = [to_float(x) for x in d["xCoordinates"].split()]
                ys = [to_float(x) for x in d.get("yCoordinates", "").split()]
                s["_xy"] = list(zip(xs, ys, [""] * len(xs)))
            out.append(s)
    return out


# --------------------------------------------------------------------------- #
# 3D / physics / sediment
# --------------------------------------------------------------------------- #
def vertical_info(mdu) -> dict:
    kmx = int(to_float(mdu("geometry", "Kmx", "0"), 0))
    out = {"Kmx": kmx, "mode": "2D depth-averaged" if kmx == 0 else f"3D, {kmx} layers"}
    if kmx == 0:
        return out
    lt = (mdu("geometry", "Layertype", "1") or "1").strip()
    out["Layertype"] = decode(mdu, "geometry", "Layertype", "1")
    if lt == "1":
        gf = to_float(mdu("geometry", "SigmaGrowthFactor", "1"), 1.0)
        th = [gf ** k for k in range(kmx)]
        s = sum(th)
        out["SigmaGrowthFactor"] = g(gf)
        out["relative_thickness"] = " : ".join(f"{100 * t / s:.1f}%" for t in th) + \
                                    ("  (uniform)" if gf == 1 else "  (ratio set by SigmaGrowthFactor)")
        zkeys = [k for k in ("DzTop", "FloorLevTopLay", "DzTopUniAboveZ", "NumTopSig", "NumTopSigUniform")
                 if mdu("geometry", k) not in (None, "", "0")]
        if zkeys:
            out["note"] = f"{', '.join(zkeys)} set, but these only apply to z-layers (Layertype 2)"
    else:
        for k in ("DzTop", "FloorLevTopLay", "DzTopUniAboveZ", "NumTopSig", "NumTopSigUniform",
                  "SigmaGrowthFactor", "ZlayBot", "ZlayTop", "VertplizFile"):
            if mdu("geometry", k) not in (None, ""):
                out[k] = mdu("geometry", k)
    for k in ("StretchType", "StretchCoef"):
        if mdu("geometry", k):
            out[k] = mdu("geometry", k)
    out["Turbulencemodel"] = decode(mdu, "numerics", "Turbulencemodel")
    out["Turbulenceadvection"] = mdu("numerics", "Turbulenceadvection")
    out["Vicoww / Dicoww"] = f"{mdu('physics', 'Vicoww')} / {mdu('physics', 'Dicoww')}"
    return out


def physics_info(mdu) -> dict:
    sal = mdu("physics", "Salinity", "0")
    tem = mdu("physics", "Temperature", "0")
    out = {
        "Salinity": decode(mdu, "physics", "Salinity", "0") + (f"; initial {mdu('physics', 'InitialSalinity')} ppt" if sal == "1" else ""),
        "Temperature": decode(mdu, "physics", "Temperature", "0") + (
            f"; initial {mdu('physics', 'InitialTemperature')} degC, Secchi {mdu('physics', 'Secchidepth')} m"
            if tem not in ("0", None) else ""),
        "Density": (decode(mdu, "physics", "Idensform") or "not set") + ("" if sal == "1" or tem not in ("0", None) else " (no salinity/temperature, so density is uniform)"),
        "Background salinity / temperature": f"{mdu('physics', 'Backgroundsalinity')} ppt / {mdu('physics', 'Backgroundwatertemperature')} degC",
        "Friction (uniform)": f"{decode(mdu, 'physics', 'UnifFrictType', '1')}, coef {mdu('physics', 'UnifFrictCoef')}",
        "Horizontal viscosity / diffusivity": f"{mdu('physics', 'Vicouv')} / {mdu('physics', 'Dicouv')} m2/s"
                                              + (f", Smagorinsky {mdu('physics', 'Smagorinsky')}" if to_float(mdu("physics", "Smagorinsky"), 0) else ""),
        "Wall roughness": decode(mdu, "physics", "Irov"),
        "Secondary flow": decode(mdu, "physics", "SecondaryFlow"),
        "Tidal forcing (astronomic, in-domain)": decode(mdu, "physics", "TidalForcing"),
        "Wind drag": (decode(mdu, "wind", "ICdtyp") or "not set") + (f"; Cd {mdu('wind', 'Cdbreakpoints')} at {mdu('wind', 'Windspeedbreakpoints')} m/s"
                                                       if mdu("wind", "Cdbreakpoints") else ""),
        "Coriolis (AngLat)": mdu("geometry", "AngLat"),
    }
    return {k: v for k, v in out.items() if v not in (None, "", "None")}


def numerics_info(mdu) -> dict:
    out = {
        "TimeStepType": decode(mdu, "numerics", "TimeStepType"),
        "CFLMax / DtMax / DtInit": f"{mdu('numerics', 'CFLMax')} / {mdu('time', 'DtMax')} s / {mdu('time', 'DtInit')} s",
        "AdvecType": decode(mdu, "numerics", "AdvecType"),
        "Icgsolver": decode(mdu, "numerics", "Icgsolver"),
        "Limtypmom / Limtypsa": f"{mdu('numerics', 'Limtypmom')} / {mdu('numerics', 'Limtypsa')}",
        "Teta0": mdu("numerics", "Teta0"),
        "Epshu (dry threshold)": mdu("numerics", "Epshu"),
        "BedlevType": decode(mdu, "geometry", "BedlevType"),
        "Conveyance2D": decode(mdu, "geometry", "Conveyance2D"),
        "FixedWeirScheme": decode(mdu, "numerics", "FixedWeirScheme"),
        "Tlfsmo (boundary spin-up)": fmt_dur(mdu("numerics", "Tlfsmo")),
    }
    return {k: v for k, v in out.items() if v not in (None, "")}


def parse_d3d_ini(path: Path) -> list[dict]:
    """Delft3D-4 style .mor / .sed: 'Key = value [unit] description', strings in #...#."""
    blocks, cur = [], None
    for raw in read_text(path).splitlines():
        s = raw.strip()
        if not s or s.startswith("*"):
            continue
        if s.startswith("[") and "]" in s:
            cur = {"section": s[1:s.index("]")], "items": []}
            blocks.append(cur)
            continue
        if "=" not in s or cur is None:
            continue
        k, rest = s.split("=", 1)
        rest = rest.strip()
        unit, desc = "", ""
        if rest.startswith("#"):
            j = rest.find("#", 1)
            val, tail = rest[1:j], rest[j + 1:]
        else:
            m = re.match(r"(\S+)(.*)", rest)
            val, tail = (m.group(1), m.group(2)) if m else ("", "")
        tail = tail.strip()
        m = re.match(r"\[([^\]]*)\]\s*(.*)", tail)
        if m:
            unit, desc = m.group(1), m.group(2)
        else:
            desc = tail
        cur["items"].append((k.strip(), val.strip(), f"[{unit}] {desc}".strip() if unit else desc))
    return blocks


def sediment_info(mdu, warnings) -> dict | None:
    if not mdu("sediment", "MorFile") and not mdu("sediment", "SedFile"):
        return None
    out = {"Sedimentmodelnr": decode(mdu, "sediment", "Sedimentmodelnr")}
    mf = find_ci(mdu.dir, mdu("sediment", "MorFile") or "")
    if mdu("sediment", "MorFile") and mf is None:
        warnings.append(f"MorFile '{mdu('sediment', 'MorFile')}' not found")
    if mf:
        bl = parse_d3d_ini(mf)
        out["morphology"] = [{"section": b["section"], "items": b["items"]} for b in bl
                             if "fileinformation" not in b["section"].lower()]
        mor = {k.lower(): v for b in bl if b["section"].lower() == "morphology" for k, v, _ in b["items"]}
        out["MorFac"] = mor.get("morfac")
        out["MorStt (spin-up)"] = mor.get("morstt")
        out["BedUpd"] = mor.get("bedupd")
        out["Bed / Sus transport factors"] = f"{mor.get('bed', '-')} / {mor.get('sus', '-')} (wave-related {mor.get('bedw', '-')} / {mor.get('susw', '-')})"
    sf = find_ci(mdu.dir, mdu("sediment", "SedFile") or "")
    if mdu("sediment", "SedFile") and sf is None:
        warnings.append(f"SedFile '{mdu('sediment', 'SedFile')}' not found")
    if sf:
        fr = []
        for b in parse_d3d_ini(sf):
            if b["section"].lower() != "sediment":
                continue
            d = {k: v for k, v, _ in b["items"]}
            e = {"name": d.get("Name"), "type": d.get("SedTyp"), "D50 (m)": d.get("SedDia"),
                 "RhoSol": d.get("RhoSol"), "CDryB": d.get("CDryB"), "TraFrm": d.get("TraFrm"),
                 "IniSedThick": d.get("IniSedThick"), "_all": b["items"]}
            if d.get("IniSedThick") and to_float(d["IniSedThick"]) is None:
                tf = find_ci(sf.parent, d["IniSedThick"])
                if tf is None:
                    warnings.append(f"sediment {d.get('Name')}: IniSedThick file '{d['IniSedThick']}' not found")
                else:
                    vals = [to_float(r.split()[2]) for r in read_text(tf).splitlines()
                            if len(r.split()) >= 3 and to_float(r.split()[2]) is not None]
                    s = vstats(vals)
                    if s:
                        e["IniSedThick"] += f" ({s['n']} samples, {g(s['min'])} .. {g(s['max'])} m, mean {g(s['mean'])})"
            fr.append(e)
        out["fractions"] = fr
    return out


# --------------------------------------------------------------------------- #
# D-Waves (SWAN) model: .mdw, curvilinear .grd/.dep, obstacles
# --------------------------------------------------------------------------- #
def read_rgf_grid(fp: Path) -> dict:
    lines = read_text(fp).splitlines()
    i, missing, cs = 0, 0.0, None
    while i < len(lines) and (lines[i].startswith("*") or "=" in lines[i] and "ETA" not in lines[i]):
        if "Missing Value" in lines[i]:
            missing = to_float(lines[i].split("=")[1], 0.0)
        if "Coordinate System" in lines[i]:
            cs = lines[i].split("=")[1].strip()
        i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    m, n = (int(x) for x in lines[i].split()[:2])
    i += 2                                          # skip "0 0 0"
    vals = []
    for ln in lines[i:]:
        t = ln.split()
        if t and t[0].upper().startswith("ETA"):
            t = t[2:] if t[0].upper() == "ETA=" else t[1:]
        vals += [to_float(x) for x in t]
    X, Y = vals[:m * n], vals[m * n:2 * m * n]
    ok = lambda k: X[k] is not None and not (abs(X[k] - missing) < 1e-6 and abs(Y[k] - missing) < 1e-6) and not (X[k] == 0 and Y[k] == 0)
    good = [k for k in range(len(X)) if ok(k)]
    out = {"file": fp.name, "M x N": f"{m} x {n}", "coordinate_system": cs, "_m": m, "_n": n}
    if good:
        gx, gy = [X[k] for k in good], [Y[k] for k in good]
        out["extent_x"] = f"{g(min(gx), 9)} .. {g(max(gx), 9)} ({g(max(gx) - min(gx))} m)"
        out["extent_y"] = f"{g(min(gy), 9)} .. {g(max(gy), 9)} ({g(max(gy) - min(gy))} m)"
        dm = [math.hypot(X[k + 1] - X[k], Y[k + 1] - Y[k]) for k in good[:: max(1, len(good) // 4000)]
              if (k + 1) % m and k + 1 < len(X) and ok(k + 1)]
        dn = [math.hypot(X[k + m] - X[k], Y[k + m] - Y[k]) for k in good[:: max(1, len(good) // 4000)]
              if k + m < len(X) and ok(k + m)]
        if dm and dn:
            out["cell_size_M_N"] = f"~{g(sorted(dm)[len(dm) // 2])} x {g(sorted(dn)[len(dn) // 2])} m (median)"
        ring = [(r, 0) for r in range(m)] + [(m - 1, c) for c in range(n)] + \
               [(r, n - 1) for r in range(m - 1, -1, -1)] + [(0, c) for c in range(n - 1, -1, -1)]
        out["_outline"] = [(X[c * m + r], Y[c * m + r]) for r, c in ring if ok(c * m + r)]
        out["_bbox"] = (min(gx), max(gx), min(gy), max(gy))
    return out


def read_dep(fp: Path, m: int, n: int) -> str | None:
    vals = [v for v in (to_float(t) for t in read_text(fp).split()) if v is not None]
    good = [v for v in vals if v > -998]
    if not good:
        return None
    note = "" if len(vals) in ((m + 1) * (n + 1), m * n) else f" - {len(vals)} values, expected {(m + 1) * (n + 1)}"
    return f"{g(min(good))} .. {g(max(good))} (depth, positive down; {len(good)} valid points){note}"


def find_mdw(mdu) -> list[Path]:
    seen, out = set(), []
    for base in (mdu.dir, mdu.dir.parent, mdu.dir.parent.parent):
        for f in sorted(base.rglob("*.mdw")) if base.exists() else []:
            if f.resolve() in seen or any(x.lower() in ("output",) or x.lower().startswith("dfm_output")
                                          for x in f.parts[:-1]):
                continue
            seen.add(f.resolve())
            out.append(f)
        if out:
            break
    return out


def _num(v):
    f = to_float(v)
    return g(f) if f is not None else v


def _nums(vals):
    return ", ".join(_num(v) for v in vals)


def parse_wave(fp: Path, warnings) -> dict:
    bl = parse_ini(fp)
    d = fp.parent
    gen = {k: v for b in sections(bl, "General") for k, v, _ in b["items"]}
    proc = {k: v for b in sections(bl, "Processes") for k, v, _ in b["items"]}
    out = {"file": str(fp), "general": gen, "processes": proc,
           "numerics": {k: v for b in sections(bl, "Numerics") for k, v, _ in b["items"]},
           "output": {k: v for b in sections(bl, "Output") for k, v, _ in b["items"]},
           "constants": {k: v for b in sections(bl, "Constants") for k, v, _ in b["items"]},
           "domains": [], "boundaries": [], "obstacles": [], "time_points": []}
    for i, b in enumerate(sections(bl, "Domain")):
        dd = bdict(b)
        dom = {"#": i + 1, "grid": dd.get("Grid"), "bed": dd.get("BedLevel"),
               "nested_in": dd.get("NestedInDomain"), "output": dd.get("Output")}
        gf = find_ci(d, dd.get("Grid", ""))
        if gf is None:
            warnings.append(f"wave grid '{dd.get('Grid')}' not found")
        else:
            try:
                gr = read_rgf_grid(gf)
                dom.update({k: v for k, v in gr.items() if k != "file"})
                bf = find_ci(d, dd.get("BedLevel", ""))
                if bf is None:
                    warnings.append(f"wave bed level '{dd.get('BedLevel')}' not found")
                else:
                    dom["depth"] = read_dep(bf, gr["_m"], gr["_n"])
            except (OSError, ValueError, IndexError) as e:
                warnings.append(f"wave grid {gf.name}: {e}")
        out["domains"].append(dom)
    for b in sections(bl, "Boundary"):
        dd = bdict(b)
        e = {"name": dd.get("Name"), "definition": dd.get("Definition"),
             "spectrum": dd.get("SpectrumSpec"), "shape": dd.get("SpShapeType"),
             "Hs": _nums(bgetall(b, "WaveHeight")), "period": _nums(bgetall(b, "Period")),
             "period_type": dd.get("PeriodType"), "direction": _nums(bgetall(b, "Direction")),
             "spreading": _nums(bgetall(b, "DirSpreading")), "gamma": _num(dd.get("PeakEnhanceFac")),
             "orientation": dd.get("Orientation"), "spectrum_file": dd.get("Spectrum")}
        if dd.get("StartCoordX"):
            e["_seg"] = ((to_float(dd["StartCoordX"]), to_float(dd["StartCoordY"])),
                         (to_float(dd["EndCoordX"]), to_float(dd["EndCoordY"])))
            e["segment"] = f"({g(e['_seg'][0][0], 9)}, {g(e['_seg'][0][1], 9)}) -> ({g(e['_seg'][1][0], 9)}, {g(e['_seg'][1][1], 9)})"
        out["boundaries"].append(e)
    for b in sections(bl, "TimePoint"):
        out["time_points"].append(bdict(b))
    of = gen.get("ObstacleFile")
    if of:
        ofp = find_ci(d, of)
        if ofp is None:
            warnings.append(f"wave ObstacleFile '{of}' not found")
        else:
            obl = parse_ini(ofp)
            pol = next((v for b in obl for k, v, _ in b["items"] if k.lower() == "polylinefile"), None)
            lines = {}
            if pol and find_ci(d, pol):
                lines = {p["name"]: p["_xy"] for p in parse_pli(find_ci(d, pol))}
            elif pol:
                warnings.append(f"obstacle polyline file '{pol}' not found")
            for b in sections(obl, "Obstacle"):
                dd = bdict(b)
                out["obstacles"].append({"name": dd.get("Name"), "type": dd.get("Type"),
                                         "TransmCoef": _num(dd.get("TransmCoef")), "Reflections": dd.get("Reflections"),
                                         "_xy": lines.get(dd.get("Name"), [])})
    return out


# --------------------------------------------------------------------------- #
# Time
# --------------------------------------------------------------------------- #
def _dt14(s):
    try:
        s = re.sub(r"\D", "", s or "")
        return _dt.datetime.strptime((s + "000000")[:14], "%Y%m%d%H%M%S")
    except ValueError:
        return None


def time_info(mdu) -> dict:
    ref = ref_datetime(mdu)
    f = tunit_sec(mdu)
    sd, ed = _dt14(mdu("time", "StartDateTime")), _dt14(mdu("time", "StopDateTime"))
    if sd and ed and ref:
        t0, t1, src = (sd - ref).total_seconds(), (ed - ref).total_seconds(), "StartDateTime / StopDateTime"
    else:
        t0, t1, src = to_float(mdu("time", "TStart"), 0) * f, to_float(mdu("time", "TStop"), 0) * f, "TStart / TStop"
    out = {"RefDate": mdu("time", "RefDate"), "Tzone": mdu("time", "Tzone"), "Tunit": mdu("time", "Tunit", "S"),
           "defined by": src, "start": fmt_dur(t0), "stop": fmt_dur(t1), "duration": fmt_dur(t1 - t0),
           "DtUser": mdu("time", "DtUser"), "DtMax": mdu("time", "DtMax"), "DtInit": mdu("time", "DtInit"),
           "DtNodal": mdu("time", "DtNodal"), "_t0": t0, "_t1": t1}
    if ref:
        out["period"] = f"{ref + _dt.timedelta(seconds=t0):%Y-%m-%d %H:%M} -> {ref + _dt.timedelta(seconds=t1):%Y-%m-%d %H:%M}"
    if sd and (mdu("time", "TStart") or mdu("time", "TStop")):
        out["note"] = "both StartDateTime/StopDateTime and TStart/TStop are set"
    return out


# --------------------------------------------------------------------------- #
# Build
# --------------------------------------------------------------------------- #
BND_COLOR = {"waterlevelbnd": "#d32f2f", "dischargebnd": "#2e7d32", "neumannbnd": "#ef6c00",
             "velocitybnd": "#6a1b9a", "riemannbnd": "#00838f", "qhbnd": "#5d4037"}


def model_summary(mdu, net, ti, vert, ext, feats, sed, waves) -> dict:
    m2 = net.get("mesh2d", {})
    nb = sum(len(b.get("polylines", [])) or 1 for b in ext["boundaries"]) if ext else 0
    q = {}
    for b in (ext or {}).get("boundaries", []):
        q[b["quantity"]] = q.get(b["quantity"], 0) + 1
    out = {
        "Program": re.sub(r"^Deltares,\s*", "", mdu("General", "Version") or mdu("General", "Program") or ""),
        "Delta Shell GUI": mdu("General", "GuiVersion"),
        "Simulation period": ti.get("period") or f"{ti['start']} -> {ti['stop']}",
        "Duration": ti["duration"],
        "2D mesh": (f"{m2.get('faces')} cells ({m2.get('cell_types')}), {m2.get('approx_cell_size')}, {m2.get('total_area')}"
                    if m2 else None),
        "Vertical": vert["mode"] + (f", {vert.get('Layertype')}" if vert["Kmx"] else ""),
        "Projection": net.get("projection"),
        "Bed level": (m2.get("bed_level_face_z") if "not set" not in str(m2.get("bed_level_face_z", "not set")) else m2.get("bed_level_node_z")),
        "Time step": (decode(mdu, "numerics", "TimeStepType") or "not set") + f"; CFLMax {mdu('numerics', 'CFLMax', '0.7')}, DtMax {mdu('time', 'DtMax')} s",
        "Open boundaries": f"{nb} polylines: " + ", ".join(f"{v} {k}" for k, v in q.items()) if q else None,
        "Laterals / sources": (len(ext["laterals"]) + len(ext["sources"])) if ext and (ext["laterals"] or ext["sources"]) else None,
        "Meteo forcing": ", ".join(x["quantity"] for x in ext["meteo"]) if ext and ext["meteo"] else None,
        "Salinity / temperature": f"{decode(mdu, 'physics', 'Salinity', '0')} / {decode(mdu, 'physics', 'Temperature', '0')}",
        "Friction": f"{decode(mdu, 'physics', 'UnifFrictType', '1')} {mdu('physics', 'UnifFrictCoef')} (uniform)",
        "Waves": decode(mdu, "waves", "Wavemodelnr", "0") + (f" - D-Waves model {Path(waves[0]['file']).name} found" if waves else ""),
        "Sediment / morphology": (f"{sed['Sedimentmodelnr']}; {len(sed.get('fractions', []))} fraction(s)"
                                  + (f", MorFac {sed.get('MorFac')}" if sed.get("MorFac") else "")) if sed else None,
        "Thin dams": f"{feats['thin_dams']['count']} ({feats['thin_dams']['total_length_m'] / 1000:.2f} km)" if feats.get("thin_dams") else None,
        "Fixed weirs": f"{feats['fixed_weirs']['count']}" if feats.get("fixed_weirs") else None,
        "Structures": len(feats["structures"]) or None,
        "Observation points / cross sections": f"{len(feats['obs'])} / {len(feats['crs'])}" if feats["obs"] or feats["crs"] else None,
        "Output": f"his {fmt_dur(mdu('output', 'HisInterval')) or '-'}, map {fmt_dur((mdu('output', 'MapInterval') or '').split()[0] if mdu('output', 'MapInterval') else None) or '-'}"
                  + (f", restart {fmt_dur(mdu('output', 'RstInterval'))}" if to_float(mdu("output", "RstInterval"), 0) else ""),
        "Restart from": mdu("restart", "RestartFile"),
    }
    return {k: v for k, v in out.items() if v not in (None, "", 0)}


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
    ext = None
    nf = mdu("external forcing", "ExtForceFileNew")
    if nf:
        f = find_ci(mdu.dir, nf)
        if f is None:
            warnings.append(f"ExtForceFileNew '{nf}' not found")
        else:
            ext = parse_ext_new(f, f.parent, warnings)
            for b in ext["boundaries"]:
                npts = sum(p["points"] for p in b.get("polylines", []))
                if npts:
                    nf = len(b.get("forcings", []))
                    b["support_points"] = (f"all {npts} points forced" if nf >= npts else
                                           f"{nf} of {npts} points forced, rest interpolated along the polyline")
    old = None
    of = mdu("external forcing", "ExtForceFile")
    if of:
        f = find_ci(mdu.dir, of)
        if f is None:
            warnings.append(f"ExtForceFile '{of}' not found")
        else:
            old = parse_ext_old(f, f.parent)
    feats = {"thin_dams": read_polyline_files(mdu, "ThinDamFile", warnings=warnings),
             "obs": read_obs(mdu, warnings), "crs": read_crs(mdu, warnings),
             "structures": read_structures(mdu, warnings)}
    fw = split_files(mdu("geometry", "FixedWeirFile"))
    if fw:
        pl = []
        for fn in fw:
            p = find_ci(mdu.dir, fn)
            if p is None:
                warnings.append(f"FixedWeirFile '{fn}' not found")
            else:
                pl += parse_pliz(p)
        crest = [c for p in pl for c in p.get("_crest", [])]
        feats["fixed_weirs"] = {"files": fw, "count": len(pl), "polylines": pl,
                                "crest_range": f"{g(min(crest))} .. {g(max(crest))} m" if crest else None,
                                "total_length_m": sum(pli_length(p["_xy"]) for p in pl)}
    for key in ("DryPointsFile", "GridEnclosureFile", "LandBoundaryFile", "PartitionFile"):
        for fn in split_files(mdu("geometry", key)):
            if find_ci(mdu.dir, fn) is None:
                warnings.append(f"{key} '{fn}' not found")
    m["time"] = ti
    m["network"] = net
    m["vertical"] = vert = vertical_info(mdu)
    m["physics"] = physics_info(mdu)
    m["numerics"] = numerics_info(mdu)
    m["external_forcing"] = ext
    m["external_forcing_old"] = old
    m["features"] = feats
    m["initial_fields"] = parse_initial_fields(mdu, warnings)
    m["sediment"] = sed = sediment_info(mdu, warnings)
    waves = []
    if (mdu("waves", "Wavemodelnr", "0") or "0") != "0":
        for w in find_mdw(mdu):
            try:
                waves.append(parse_wave(w, warnings))
            except OSError as e:
                warnings.append(str(e))
        if mdu("waves", "Wavemodelnr") == "3" and not waves:
            warnings.append("Wavemodelnr = 3 (SWAN) but no .mdw wave model found next to the flow model")
    m["waves"] = waves
    m["files"] = referenced_files(mdu)
    for f in m["files"]:
        if not f["exists"] and not any(f["file"] in w for w in warnings):
            warnings.append(f"{f['key']} '{f['file']}' not found")
    m["summary"] = model_summary(mdu, net, ti, vert, ext, feats, sed, waves)
    m["mdu_parameters"] = mdu_groups(mdu)
    m["results"] = results(mdu, ti, warnings)
    # checks
    if ext:
        for b in ext["boundaries"]:
            if not b.get("forcings"):
                warnings.append(f"boundary {b['location']}: no matching forcing found")
            for fo in b.get("forcings", []):
                if fo.get("_t_end") is not None and fo["_t_end"] < ti["_t1"] - 1e-6:
                    warnings.append(f"forcing {fo['name']} ends at {g(fo['_t_end'])} s, before the stop time ({g(ti['_t1'])} s)")
    if vert.get("note"):
        warnings.append(vert["note"])
    if vert["Kmx"] and mdu("physics", "Salinity", "0") == "0" and mdu("physics", "Temperature", "0") == "0":
        warnings.append("3D model without salinity or temperature: no density stratification (check this is intended)")
    if waves and mdu("waves", "Wavemodelnr") == "3":
        com = waves[0]["output"].get("COMFile")
        if com and Path(com).stem.replace("_com", "") != path.stem:
            warnings.append(f"D-Waves COMFile '{com}' does not match the flow model name '{path.stem}'")
    m["warnings"] = warnings
    return m


def build_config(arg: str, only=None) -> dict:
    mdus = find_mdus(arg)
    sel = [p for p in mdus if not only or p.stem.lower() == only.lower()]
    if not sel:
        sys.exit(f"Model '{only}' not found. MDUs: {', '.join(p.stem for p in mdus)}")
    root = Path(arg).expanduser().resolve()
    root = root if root.is_dir() else root.parent
    proj = next(iter(sorted(root.glob("*.dsproj"))), None)
    if proj is None:
        for q in sel[0].parents:
            if q.name.endswith(".dsproj_data"):
                proj = q.with_name(q.name[:-5])
                break
    return {"name": Path(str(proj)).stem if proj else root.name, "root": str(root),
            "dsproj": str(proj) if proj else None, "models": [build_model(p) for p in sel]}


# --------------------------------------------------------------------------- #
# Terminal output
# --------------------------------------------------------------------------- #
def show_model(m, all_params=False):
    h1(f"D-FLOW FM 2D/3D MODEL: {m['name']}")
    row("MDU file", m["mdu"])
    h2("Model summary")
    for k, v in m["summary"].items():
        row(k, v)
    net = m["network"]
    h2(f"Mesh ({net.get('file', '-')})")
    if net.get("error"):
        print(f"  ERROR: {net['error']}")
    row("Projection", net.get("projection"))
    for k, v in (net.get("mesh2d") or {}).items():
        if k != "name":
            row(k.replace("_", " "), v)
    if net.get("branches"):
        row("1D branches", f"{len(net['branches'])} (use dfmvue for 1D detail)")
    h2("Vertical schematisation")
    for k, v in m["vertical"].items():
        row(k, v)
    h2("Physics")
    for k, v in m["physics"].items():
        row(k, v)
    h2("Numerics")
    for k, v in m["numerics"].items():
        row(k, v)
    ext = m["external_forcing"]
    h2("Boundaries & forcing")
    if ext:
        for b in ext["boundaries"]:
            pl = ", ".join(f"{p['name']} ({p['points']} pts, {pli_length(p['_xy']) / 1000:.2f} km)" for p in b.get("_pli", []))
            print(f"  BND {b['quantity']:<15} {pl}  <- {b['forcing_file']}   [{b.get('support_points', '')}]")
            for f in b.get("forcings", []):
                print(f"      {f['name']:<18} {f['function']:<10} {f.get('time_span', '')} {f.get('value_range', '')}")
        for L in ext["laterals"]:
            print(f"  LAT {L['id']} ({L['type'] or ''}) @ {L['location']}  discharge = {L['discharge']}")
        for mt in ext["meteo"]:
            print(f"  MET {mt['quantity']} <- {mt['forcingFile']} ({mt['forcingFileType']}){'' if mt['exists'] else '  MISSING'}")
    for e in m["external_forcing_old"] or []:
        print(f"  OLD {e.get('QUANTITY')} <- {e.get('FILENAME')} ({e.get('FILETYPE')}, method {e.get('METHOD')}, operand {e.get('OPERAND')})")
    if not ext and not m["external_forcing_old"]:
        print("  no external forcing file")
    F = m["features"]
    h2("Thin dams, weirs, structures, observations")
    if F.get("thin_dams"):
        t = F["thin_dams"]
        row("Thin dams", f"{t['count']} polylines in {len(t['files'])} file(s), {t['total_length_m'] / 1000:.2f} km")
    if F.get("fixed_weirs"):
        row("Fixed weirs", f"{F['fixed_weirs']['count']} polylines, crest {F['fixed_weirs']['crest_range']}")
    for s in F["structures"]:
        print(f"  STR {s['type'] or '?':<12} {s['id']:<20} @ {s['location']}")
    if F["obs"]:
        row("Observation points", ", ".join(o["name"] for o in F["obs"]))
    if F["crs"]:
        row("Cross sections", ", ".join(c["name"] for c in F["crs"]))
    if m["initial_fields"]:
        h2("Initial fields & parameters")
        for f in m["initial_fields"]:
            print(f"  [{f['block']}] {f['quantity']:<22} <- {f['dataFile']} ({f['dataFileType']}{', ' + f['method'] if f['method'] else ''})")
            if f.get("content"):
                print(f"      {f['content']}")
    s = m["sediment"]
    if s:
        h2("Sediment & morphology")
        for k in ("Sedimentmodelnr", "MorFac", "MorStt (spin-up)", "BedUpd", "Bed / Sus transport factors"):
            row(k, s.get(k))
        for fr in s.get("fractions", []):
            print(f"  fraction {fr['name']} ({fr['type']}): D50 {fr['D50 (m)']} m, RhoSol {fr['RhoSol']}, "
                  f"CDryB {fr['CDryB']}, TraFrm {fr['TraFrm']}")
            row("IniSedThick", fr["IniSedThick"], indent=4)
    for w in m["waves"]:
        h2(f"Waves: D-Waves {Path(w['file']).name}")
        gnl = w["general"]
        row("Mode / time step", f"{gnl.get('SimMode')}, TimeStep {gnl.get('TimeStep')}, TimeInterval {gnl.get('TimeInterval')}")
        row("Flow coupling (bed/wl/vel/wind)", "/".join(gnl.get(k, "-") for k in ("FlowBedLevel", "FlowWaterLevel", "FlowVelocity", "FlowWind")))
        row("Spectral resolution", f"{gnl.get('NDir')} directions ({gnl.get('DirSpace')}), {gnl.get('NFreq')} freq {gnl.get('FreqMin')}-{gnl.get('FreqMax')} Hz")
        row("Wind (uniform)", f"{gnl.get('WindSpeed')} m/s from {gnl.get('WindDir')} deg")
        p = w["processes"]
        row("Physics", f"GEN{p.get('GenModePhys')}, breaking {p.get('Breaking')} (gamma {p.get('BreakGamma')}), "
                       f"friction {p.get('BedFriction')} {p.get('BedFricCoef')}, triads {p.get('Triads')}, "
                       f"whitecapping {p.get('WhiteCapping')}, diffraction {p.get('Diffraction')}, setup {p.get('WaveSetup')}")
        row("Wave forces", p.get("WaveForces"))
        row("COM file", w["output"].get("COMFile"))
        for dm in w["domains"]:
            print(f"  DOMAIN {dm['#']}: {dm['grid']} {dm.get('M x N', '')}  {dm.get('cell_size_M_N', '')}"
                  + (f"  nested in {dm['nested_in']}" if dm.get("nested_in") else ""))
            row("depth", dm.get("depth"), indent=6)
        for b in w["boundaries"]:
            print(f"  BND {b['name']:<14} {b['spectrum']}/{b['shape']}: Hs {b['Hs']} m, T{(b['period_type'] or 'p')[0]} {b['period']} s, "
                  f"dir {b['direction']}, spreading {b['spreading']}")
        if w["obstacles"]:
            row("Obstacles", f"{len(w['obstacles'])} ({', '.join(sorted({o['type'] or '?' for o in w['obstacles']}))}, "
                             f"transmission {', '.join(sorted({o['TransmCoef'] or '?' for o in w['obstacles']}))})")
    t = m["time"]
    h2("Time")
    for k in ("RefDate", "defined by", "start", "stop", "period", "duration", "DtUser", "DtMax", "DtNodal", "note"):
        row(k, t.get(k))
    h2("MDU parameters" + ("" if all_params else " (non-empty; Wri* output switches summarised; --all for every key)"))
    for sec, items in m["mdu_parameters"].items():
        it = [e for e in items if all_params or e["value"] != ""]
        on = []
        if not all_params:
            on = [e["key"] for e in it if re.match(r"wri(his|map|shp)_", e["key"], re.I) and e["value"] == "1"]
            it = [e for e in it if not re.match(r"wri(his|map|shp)_", e["key"], re.I)]
        if not it and not on:
            continue
        print(f"  [{sec}]")
        wrap_pairs([(e["key"], e["value"]) for e in it])
        if on:
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
            for k in ("Computation started at", "Computation finished at", "total computation time (s)"):
                row(k, d.get(k))
            for x in d["errors"]:
                print(f"    ERROR   {x}")
            for x in d["warnings"][:20]:
                print(f"    warning {x}")
        for kind in ("map", "his"):
            e = r.get(kind)
            if e:
                print(f"  {kind.upper()} {e['file']} ({e['size_MB']} MB)")
                for k, v in e.items():
                    if k not in ("file", "size_MB", "variables"):
                        row(k, v, indent=4)
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
# Plan view: mesh rasterised to an embedded PNG (stdlib only) + SVG vectors
# --------------------------------------------------------------------------- #
RAMP_SEA = [(0, (8, 48, 107)), (0.5, (66, 146, 198)), (1, (198, 219, 239))]
RAMP_LAND = [(0, (229, 245, 224)), (0.4, (161, 217, 155)), (0.75, (199, 162, 104)), (1, (120, 84, 50))]


def _rgb(hexc):
    return tuple(int(hexc[i:i + 2], 16) for i in (1, 3, 5))


def bed_color(v, lo, hi):
    if lo < 0 < hi:
        return ramp((v - lo) / (0 - lo), RAMP_SEA) if v < 0 else ramp(v / hi, RAMP_LAND)
    return ramp((v - lo) / ((hi - lo) or 1), RAMP_BED)


def png_bytes(w, h, buf) -> bytes:
    raw = b"".join(b"\x00" + bytes(buf[r * w * 4:(r + 1) * w * 4]) for r in range(h))

    def chunk(t, d):
        return struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d) & 0xFFFFFFFF)
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 6)) + chunk(b"IEND", b""))


def rasterize(polys, colors, X, Y, W, H) -> str:
    buf = bytearray(W * H * 4)
    fast = len(polys) > 400000

    def put(i, j, c):
        if 0 <= i < W and 0 <= j < H:
            k = 4 * (j * W + i)
            buf[k:k + 4] = bytes((c[0], c[1], c[2], 255))
    for (px, py), col in zip(polys, colors):
        if col is None or not px:
            continue
        c = _rgb(col)
        xs, ys = [X(a) for a in px], [Y(b) for b in py]
        x0, x1, y0, y1 = int(min(xs)), int(max(xs)), int(min(ys)), int(max(ys))
        cx, cy = sum(xs) / len(xs), sum(ys) / len(ys)
        if fast or (x1 - x0 <= 1 and y1 - y0 <= 1):
            put(int(cx), int(cy), c)
            continue
        n = len(xs)
        for j in range(max(0, y0), min(H - 1, y1) + 1):
            yc = j + 0.5
            for i in range(max(0, x0), min(W - 1, x1) + 1):
                xc, inside = i + 0.5, False
                for a in range(n):
                    b = a - 1
                    if (ys[a] > yc) != (ys[b] > yc) and xc < (xs[b] - xs[a]) * (yc - ys[a]) / (ys[b] - ys[a]) + xs[a]:
                        inside = not inside
                if inside:
                    put(i, j, c)
        put(int(cx), int(cy), c)
    return base64.b64encode(png_bytes(W, H, buf)).decode()


def planview(m, values=None, kind="bed", label="bed level (m)", extent="flow", w=900) -> str:
    P = m["network"].get("_plot", {})
    faces = P.get("faces") or []
    ext = m["external_forcing"] or {}
    bnds = [(b["quantity"], pl) for b in ext.get("boundaries", []) for pl in b.get("_pli", [])]
    F = m["features"]
    wave = m["waves"][0] if m["waves"] else None
    xs, ys = [], []
    if extent == "wave" and wave:
        for dm in wave["domains"]:
            if dm.get("_bbox"):
                xs += dm["_bbox"][:2]
                ys += dm["_bbox"][2:]
    if not xs:
        for px, py in faces[::max(1, len(faces) // 20000)]:
            xs += px
            ys += py
        for _, pl in bnds:
            xs += [p[0] for p in pl["_xy"]]
            ys += [p[1] for p in pl["_xy"]]
    if not xs:
        return ""
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    pad = max(x1 - x0, y1 - y0) * 0.02 or 1
    x0, x1, y0, y1 = x0 - pad, x1 + pad, y0 - pad, y1 + pad
    sc = w / (x1 - x0)
    h = (y1 - y0) * sc
    if h > 900:
        sc = 900 / (y1 - y0)
        h, w = 900, (x1 - x0) * sc
    W, H = int(w), int(h)
    X = lambda v: (v - x0) * sc
    Y = lambda v: (y1 - v) * sc
    o = [f'<svg class="plan" viewBox="0 0 {W} {H}">']
    vv = [v for v in (values or []) if v is not None]
    lo, hi = (min(vv), max(vv)) if vv else (0, 1)
    if faces:
        if values:
            cols = [(bed_color(v, lo, hi) if kind == "bed" else ramp((v - lo) / ((hi - lo) or 1), RAMP_DEPTH))
                    if v is not None and not (kind == "depth" and v <= 0.01) else None for v in values]
        else:
            cols = ["#c9ced6"] * len(faces)
        png = rasterize(faces, cols, X, Y, W, H)
        o.append(f'<image href="data:image/png;base64,{png}" x="0" y="0" width="{W}" height="{H}"/>')
    pts = lambda L: " ".join(f"{X(p[0]):.1f},{Y(p[1]):.1f}" for p in L if p[0] is not None)
    if wave:
        for dm in wave["domains"]:
            if dm.get("_outline"):
                o.append(f'<polygon points="{pts(dm["_outline"])}" class="wdom"><title>wave domain {dm["#"]}: {_h.escape(str(dm["grid"]))}</title></polygon>')
        for b in wave["boundaries"]:
            if b.get("_seg"):
                o.append(f'<polyline points="{pts(b["_seg"])}" class="wbnd"><title>wave boundary {_h.escape(str(b["name"]))}: Hs {b["Hs"]} m</title></polyline>')
        for ob in wave["obstacles"]:
            if ob["_xy"]:
                o.append(f'<polyline points="{pts(ob["_xy"])}" class="wobs"><title>wave obstacle {_h.escape(str(ob["name"]))}</title></polyline>')
    if F.get("thin_dams"):
        for p in F["thin_dams"]["polylines"]:
            o.append(f'<polyline points="{pts(p["_xy"])}" class="thd"><title>thin dam {_h.escape(p["name"])} ({p["file"]})</title></polyline>')
    if F.get("fixed_weirs"):
        for p in F["fixed_weirs"]["polylines"]:
            o.append(f'<polyline points="{pts(p["_xy"])}" class="fxw"><title>fixed weir {_h.escape(p["name"])}</title></polyline>')
    for s in F["structures"]:
        if s.get("_xy"):
            o.append(f'<polyline points="{pts(s["_xy"])}" class="str"><title>{_h.escape(str(s["type"]))} {_h.escape(str(s["id"]))}</title></polyline>')
    for c in F["crs"]:
        o.append(f'<polyline points="{pts(c["_xy"])}" class="crs"><title>cross section {_h.escape(str(c["name"]))}</title></polyline>')
    for q, pl in bnds:
        col = BND_COLOR.get(q, "#455a64")
        o.append(f'<polyline points="{pts(pl["_xy"])}" fill="none" stroke="{col}" stroke-width="5" stroke-linecap="round">'
                 f'<title>{_h.escape(pl["name"])}: {_h.escape(q)}</title></polyline>')
        a = pl["_xy"][0]
        o.append(f'<text x="{X(a[0]) + 6:.1f}" y="{Y(a[1]) - 6:.1f}" class="blab" fill="{col}">{_h.escape(pl["name"])}</text>')
    for ob in F["obs"]:
        if ob["x"] is None:
            continue
        o.append(f'<circle cx="{X(ob["x"]):.1f}" cy="{Y(ob["y"]):.1f}" r="3.5" class="obs"><title>{_h.escape(ob["name"])}</title></circle>')
        if len(F["obs"]) <= 40:
            o.append(f'<text x="{X(ob["x"]) + 5:.1f}" y="{Y(ob["y"]) + 3:.1f}" class="olab">{_h.escape(ob["name"])}</text>')
    o.append("</svg>")
    if vv and kind == "bed":
        grad = (",".join(ramp(t, RAMP_SEA) for t in (0, .5, 1)) + "," + ",".join(ramp(t, RAMP_LAND) for t in (0, .4, .75, 1))
                if lo < 0 < hi else ",".join(ramp(t, RAMP_BED) for t in (0, .25, .5, .75, 1)))
    else:
        grad = ",".join(ramp(t, RAMP_DEPTH) for t in (0, .3, .65, 1))
    leg = [f'<span><i class="sw" style="background:linear-gradient(90deg,{grad})"></i>{_h.escape(label)}: {g(lo)} .. {g(hi)}</span>'
           if vv else '<span>mesh (no values)</span>']
    for q in dict.fromkeys(q for q, _ in bnds):
        leg.append(f'<span><i class="ln" style="border-color:{BND_COLOR.get(q, "#455a64")}"></i>{_h.escape(q)}</span>')
    if F.get("thin_dams"):
        leg.append('<span><i class="ln" style="border-color:#111"></i>thin dam</span>')
    if F.get("fixed_weirs"):
        leg.append('<span><i class="ln" style="border-color:#8d6e63"></i>fixed weir</span>')
    if F["obs"]:
        leg.append('<span><i class="dot"></i>observation point</span>')
    if wave:
        leg.append('<span><i class="ln" style="border-color:#7b1fa2;border-top-style:dashed"></i>wave domain / boundary</span>')
        if wave["obstacles"]:
            leg.append('<span><i class="ln" style="border-color:#616161;border-top-style:dotted"></i>wave obstacle</span>')
    return "".join(o) + f'<div class="legend">{"".join(leg)}</div>'


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #
EXTRA_CSS = """
svg.plan{width:100%;max-width:980px;height:auto;display:block;margin:8px 0;background:var(--bg);border:1px solid var(--line);border-radius:6px}
svg.plan .thd{fill:none;stroke:#111;stroke-width:2.2} svg.plan .fxw{fill:none;stroke:#8d6e63;stroke-width:2.4}
svg.plan .str{fill:none;stroke:#c2185b;stroke-width:3} svg.plan .crs{fill:none;stroke:#00897b;stroke-width:1.6;stroke-dasharray:4 3}
svg.plan .obs{fill:#fff;stroke:#c62828;stroke-width:1.6} svg.plan .olab{font-size:10px;fill:#222;paint-order:stroke;stroke:#fff;stroke-width:2.5px}
svg.plan .blab{font-size:11px;font-weight:600;paint-order:stroke;stroke:#fff;stroke-width:3px}
svg.plan .wdom{fill:none;stroke:#7b1fa2;stroke-width:1.6;stroke-dasharray:7 4} svg.plan .wbnd{fill:none;stroke:#7b1fa2;stroke-width:4;opacity:.7}
svg.plan .wobs{fill:none;stroke:#616161;stroke-width:2;stroke-dasharray:2 2}
.legend{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:12px;color:var(--muted);margin-bottom:8px}
.legend i{display:inline-block;vertical-align:middle;margin-right:5px}
.legend .sw{width:110px;height:10px;border-radius:2px} .legend .ln{width:22px;height:0;border-top:3px solid}
.legend .dot{width:8px;height:8px;border-radius:50%;border:1.6px solid #c62828;background:#fff}
"""


def _model_section(m, mi):
    pid = lambda k: f"m{mi}-{k}"
    present, out = [], [f'<section id="m{mi}"><h2>{_e(m["name"])}']
    net, ext, F, r = m["network"], m["external_forcing"], m["features"], m["results"]
    out.append(f' <span class="badge o">{"3D · " + str(m["vertical"]["Kmx"]) + " layers" if m["vertical"]["Kmx"] else "2D"}</span>')
    if m["waves"]:
        out.append(' <span class="badge o">+ D-Waves</span>')
    if m["sediment"]:
        out.append(' <span class="badge o">+ morphology</span>')
    if r and r.get("dia"):
        ok = r["dia"]["finished"] and not r["dia"]["errors"]
        out.append(f' <span class="badge{"" if ok else " w"}">{"run finished" if ok else "run problem"}</span>')
    elif not r:
        out.append(' <span class="badge o">not run</span>')
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

    P = net.get("_plot", {})
    vals, lab = None, "bed level (m)"
    if r and r.get("_plot", {}).get("bl_2d"):
        vals, lab = r["_plot"]["bl_2d"], "bed level used by model (m)"
    elif P.get("face_z") and any(v is not None for v in P["face_z"]):
        vals = P["face_z"]
    pairs = [("file", f"{net.get('file')} ({net.get('size_kB')} kB)"), ("created by", net.get("created_by")),
             ("projection", net.get("projection")), ("conventions", net.get("conventions"))]
    pairs += [(k.replace("_", " "), v) for k, v in (net.get("mesh2d") or {}).items() if k != "name"]
    add("mesh", "Mesh & bathymetry", _two(pairs) + planview(m, vals, "bed", lab), net.get("file", ""))

    add("vert", "Vertical & physics", '<div class="grid2"><div><h4>Vertical</h4>' + _kv(list(m["vertical"].items()))
        + "</div><div><h4>Physics</h4>" + _kv(list(m["physics"].items())) + "</div></div><h4>Numerics</h4>"
        + _two(list(m["numerics"].items())))

    body = ""
    if ext:
        rows, series = [], []
        for b in ext["boundaries"]:
            for f in b.get("forcings") or [{}]:
                rows.append([f'<b style="color:{BND_COLOR.get(b["quantity"], "inherit")}">{_e(b["quantity"])}</b>',
                             _e(b["location"]), _e(f.get("name")), _e(f.get("function")), _e(f.get("time_span", "")),
                             (_e(f.get("value_range")), "msg"), _e(b.get("support_points", "")), _e(b["forcing_file"])])
                for lab_, xs, ys in f.get("_series", []):
                    series.append((f"{f['name']}", xs, ys))
        body += _table(["Quantity", "Location", "Forcing point", "Function", "Time span", "Values / components",
                        "Coverage", "File"], rows)
        if series:
            body += "<h4>Time series</h4>" + svg_lines(series[:8], "time since RefDate (h)", "boundary value")
        if ext["laterals"]:
            body += "<h4>Laterals</h4>" + _table(["Id", "Type", "Location", "Location type", "Discharge"],
                                                 [[_e(L["id"]), _e(L["type"]), _e(L["location"]), _e(L["locationType"]),
                                                   _e(L["discharge"])] for L in ext["laterals"]])
        if ext["meteo"]:
            body += "<h4>Meteo</h4>" + _table(["Quantity", "File", "File type", "Method", "Operand", "Exists"],
                                              [[_e(x["quantity"]), _e(x["forcingFile"]), _e(x["forcingFileType"]),
                                                _e(x["method"]), _e(x["operand"]),
                                                "yes" if x["exists"] else '<span class="warn">MISSING</span>'] for x in ext["meteo"]])
        if ext["sources"]:
            body += "<h4>Sources / sinks</h4>" + _table(["Settings"], [[_e(s)] for s in ext["sources"]])
    if m["external_forcing_old"]:
        body += "<h4>Old-style external forcing (ExtForceFile)</h4>" + _table(
            ["Quantity", "File", "Type", "Method", "Operand", "Exists"],
            [[_e(e.get("QUANTITY")), _e(e.get("FILENAME")), _e(e.get("FILETYPE")), _e(e.get("METHOD")),
              _e(e.get("OPERAND")), "yes" if e.get("exists") else '<span class="warn">MISSING</span>']
             for e in m["external_forcing_old"]])
    add("forcing", "Boundaries & forcing", body or '<p class="muted">No external forcing.</p>', ext["file"] if ext else "")

    body = ""
    if F.get("thin_dams"):
        t = F["thin_dams"]
        body += "<h4>Thin dams</h4>" + _table(["Name", "File", "Points", "Length (m)"],
                                              [[_e(p["name"]), _e(p["file"]), _e(p["points"]), _e(g(p["length_m"]))]
                                               for p in t["polylines"]])
    if F.get("fixed_weirs"):
        fw = F["fixed_weirs"]
        body += f"<h4>Fixed weirs</h4><p>{fw['count']} polylines, crest {_e(fw['crest_range'])}, " \
                f"{fw['total_length_m'] / 1000:.2f} km total</p>"
    if F["structures"]:
        body += "<h4>Structures</h4>" + _table(["Type", "Id", "Location", "Parameters"],
                                               [[_e(s["type"]), _e(s["id"]), _e(s["location"]),
                                                 (" ".join(f"<span class=chip>{_e(k)}={_e(v)}</span>" for k, v in s["params"].items()), "msg")]
                                                for s in F["structures"]])
    if F["obs"]:
        body += "<h4>Observation points</h4>" + _table(["Name", "x", "y", "File"],
                                                       [[_e(o["name"]), _e(g(o["x"], 10)), _e(g(o["y"], 10)), _e(o["file"])]
                                                        for o in F["obs"]])
    if F["crs"]:
        body += "<h4>Observation cross sections</h4>" + _table(["Name", "Points", "File"],
                                                               [[_e(c["name"]), _e(len(c["_xy"])), _e(c["file"])] for c in F["crs"]])
    add("feat", "Thin dams, structures & observations", body,
        f"{F['thin_dams']['count'] if F.get('thin_dams') else 0} thin dams · {len(F['obs'])} obs points")

    if m["initial_fields"]:
        add("init", "Initial fields", _table(["Block", "Quantity", "Data file", "Type", "Method", "Location", "Content"],
                                             [[_e(f["block"]), _e(f["quantity"]), _e(f["dataFile"]), _e(f["dataFileType"]),
                                               _e(f["method"]), _e(f["locationType"]), (_e(f.get("content")), "msg")]
                                              for f in m["initial_fields"]]))
    s = m["sediment"]
    if s:
        body = _two([(k, s.get(k)) for k in ("Sedimentmodelnr", "MorFac", "MorStt (spin-up)", "BedUpd",
                                              "Bed / Sus transport factors")])
        if s.get("fractions"):
            body += "<h4>Sediment fractions</h4>" + _table(
                ["Name", "Type", "D50 (m)", "RhoSol", "CDryB", "TraFrm", "IniSedThick"],
                [[_e(f["name"]), _e(f["type"]), _e(f["D50 (m)"]), _e(f["RhoSol"]), _e(f["CDryB"]), _e(f["TraFrm"]),
                  (_e(f["IniSedThick"]), "msg")] for f in s["fractions"]])
            for f in s["fractions"]:
                body += _blk(f"All settings · {_e(f['name'])}", _table(["Key", "Value", "Unit / description"],
                             [[f"<b>{_e(k)}</b>", _e(v), (_e(c), "desc")] for k, v, c in f["_all"]]), open_=False)
        for b in s.get("morphology", []):
            body += _blk(f"Morphology file · [{_e(b['section'])}]", _table(["Key", "Value", "Unit / description"],
                         [[f"<b>{_e(k)}</b>", _e(v), (_e(c), "desc")] for k, v, c in b["items"]]), open_=False)
        add("sed", "Sediment & morphology", body)

    for wi, w in enumerate(m["waves"]):
        gnl, p = w["general"], w["processes"]
        body = '<div class="grid2"><div><h4>General & coupling</h4>' + _kv(
            [("file", w["file"]), ("SimMode", gnl.get("SimMode")), ("TimeStep / TimeInterval", f"{gnl.get('TimeStep')} / {gnl.get('TimeInterval')}"),
             ("Flow -> wave coupling (bed / wl / vel / wind)", "/".join(gnl.get(k, "-") for k in ("FlowBedLevel", "FlowWaterLevel", "FlowVelocity", "FlowWind"))),
             ("FlowVelocityType", gnl.get("FlowVelocityType")), ("Directions", f"{gnl.get('NDir')} ({gnl.get('DirSpace')}, {gnl.get('StartDir')}-{gnl.get('EndDir')})"),
             ("Frequencies", f"{gnl.get('NFreq')}, {gnl.get('FreqMin')}-{gnl.get('FreqMax')} Hz"),
             ("Uniform wind", f"{gnl.get('WindSpeed')} m/s, {gnl.get('WindDir')} deg ({gnl.get('DirConvention')})"),
             ("COM file", w["output"].get("COMFile")), ("COM write interval", w["output"].get("COMWriteInterval")),
             ("Hot file", w["output"].get("UseHotFile"))]) + "</div><div><h4>Processes</h4>" + _kv(list(p.items())) + "</div></div>"
        body += "<h4>Computational domains</h4>" + _table(
            ["#", "Grid", "M x N", "Cell size", "Bed level file", "Depth range", "Nested in", "Extent x", "Extent y"],
            [[_e(d["#"]), _e(d["grid"]), _e(d.get("M x N")), _e(d.get("cell_size_M_N")), _e(d["bed"]),
              (_e(d.get("depth")), "msg"), _e(d.get("nested_in") or "-"), _e(d.get("extent_x")), _e(d.get("extent_y"))]
             for d in w["domains"]])
        body += "<h4>Wave boundaries</h4>" + _table(
            ["Name", "Definition", "Spectrum", "Hs (m)", "Period (s)", "Dir (deg)", "Spreading", "Gamma", "Segment / file"],
            [[_e(b["name"]), _e(b["definition"]), _e(f"{b['spectrum']} / {b['shape']} / {b['period_type']}"),
              _e(b["Hs"]), _e(b["period"]), _e(b["direction"]), _e(b["spreading"]), _e(b["gamma"]),
              _e(b.get("segment") or b.get("spectrum_file") or b.get("orientation"))] for b in w["boundaries"]])
        if w["obstacles"]:
            body += "<h4>Obstacles</h4>" + _table(["Name", "Type", "Transmission", "Reflections", "Points"],
                                                  [[_e(o["name"]), _e(o["type"]), _e(o["TransmCoef"]), _e(o["Reflections"]),
                                                    _e(len(o["_xy"]))] for o in w["obstacles"]])
        if w["time_points"]:
            body += "<h4>Time points</h4>" + _table(list(w["time_points"][0].keys()),
                                                    [[_e(v) for v in t.values()] for t in w["time_points"]])
        if wi == 0:
            body += "<h4>Wave domains over the flow model</h4>" + planview(m, vals, "bed", lab, extent="wave")
        body += _blk("Numerics & constants", _two(list(w["numerics"].items()) + list(w["constants"].items())), open_=False)
        add(f"wave{wi}", f"Waves: D-Waves {Path(w['file']).name}", body, f"{len(w['domains'])} domain(s)")

    t = m["time"]
    add("time", "Time", _two([(k, t.get(k)) for k in ("RefDate", "Tzone", "Tunit", "defined by", "start", "stop",
                                                         "duration", "period", "DtUser", "DtMax", "DtInit", "DtNodal", "note")]))
    secs = []
    for sec, items in m["mdu_parameters"].items():
        nonempty = sum(1 for e in items if e["value"] != "")
        secs.append(_blk(_e(sec), _table(["Key", "Value", "Description"],
                                         [[f"<b>{_e(e['key'])}</b>", _e(e["value"]), (_e(e["description"]), "desc")]
                                          for e in items]), f"{len(items)} keys, {nonempty} set",
                         open_=sec.lower() in ("geometry", "physics", "time")))
    add("mdu", "MDU parameters", '<input class="pfilter" placeholder="filter keys / values / descriptions…" '
                                 'style="width:100%;max-width:420px;padding:6px 10px;border:1px solid var(--line);'
                                 'border-radius:6px;background:var(--panel);color:var(--ink);margin:4px 0 6px">' + "".join(secs),
        f'{sum(len(v) for v in m["mdu_parameters"].values())} keys', open_=False)
    add("files", "Referenced files", _table(["Section", "Key", "File", "Status"],
                                            [[_e(f["section"]), _e(f["key"]), f'<code>{_e(f["file"])}</code>',
                                              "ok" if f["exists"] else '<span class="warn">MISSING</span>'] for f in m["files"]]),
        f'{len(m["files"])} file(s)')
    if not r:
        body = '<p class="muted">No output folder found (model not run yet).</p>'
    else:
        body = _kv([("Output folder", r["output_dir"])])
        d = r.get("dia")
        if d:
            body += "<h4>Diagnostics (.dia)</h4>" + _two(
                [("file", d["file"]), ("finished", "yes" if d["finished"] else "NO"), ("errors", len(d["errors"])),
                 ("warnings", len(d["warnings"]))] + [(k, v) for k, v in d.items()
                                                      if k not in ("file", "finished", "errors", "warnings", "counts")]
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
            body += "<h4>Maximum water depth over the run</h4>" + planview(m, md, "depth", "max water depth (m)")
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
           " · ".join(x for x in (m0["summary"].get("Simulation period"), m0["summary"].get("Vertical"),
                                  (m0["summary"].get("2D mesh") or "").split(",")[0]) if x))
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
    ap = argparse.ArgumentParser(prog=prog, description="Show the configuration of a Delft3D FM (D-Flow FM) 2D/3D model.")
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
