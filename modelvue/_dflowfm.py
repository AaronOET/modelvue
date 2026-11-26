"""Delft3D FM (D-Flow FM) parsing shared by the 1D2D reader (modelvue.dfmvue) and the 2D/3D reader (modelvue.dfmvue23)."""
from __future__ import annotations
import datetime as _dt
import html as _h
import math
import re
import struct
import sys
from pathlib import Path
from ._common import PALETTE, _ticks
from ._netcdf import open_nc


def read_text(path: Path) -> str:
    data = path.read_bytes()
    if data and not data.strip(b"\0"):
        raise OSError(f"{path.name} holds only NUL bytes (cloud placeholder not downloaded?)")
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("latin-1", "replace")


def to_float(v, default=None):
    try:
        return float(str(v).strip().replace("D", "E").replace("d", "e"))
    except (TypeError, ValueError):
        return default


def g(v, nd=4):
    if v is None:
        return None
    if isinstance(v, int) or (isinstance(v, float) and v.is_integer() and abs(v) < 1e12):
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


def find_ci(folder: Path, name: str) -> Path | None:
    if not name:
        return None
    p = (folder / name.strip().replace("\\", "/"))
    if p.exists():
        return p
    if p.parent.is_dir():
        low = p.name.lower()
        for q in p.parent.iterdir():
            if q.name.lower() == low:
                return q
    return None


def vstats(vals, fill=None) -> dict:
    v = [x for x in vals if x is not None and not (isinstance(x, float) and math.isnan(x))
         and (fill is None or abs(x - fill) > 1e-9)]
    if not v:
        return {}
    return {"min": min(v), "max": max(v), "mean": sum(v) / len(v), "n": len(v)}


def parse_ini(path: Path) -> list[dict]:
    """List of blocks: {'section', 'items': [(key, value, comment)], 'data': [lines]}."""
    blocks, cur = [], {"section": "", "items": [], "data": []}
    for raw in read_text(path).splitlines():
        s = raw.strip()
        if not s or s.startswith("#") or s.startswith("*"):
            continue
        if s.startswith("[") and "]" in s:
            cur = {"section": s[1:s.index("]")].strip(), "items": [], "data": []}
            blocks.append(cur)
            continue
        if "=" in s:
            k, rest = s.split("=", 1)
            val, com = rest, ""
            if "#" in rest:
                val, com = rest.split("#", 1)
            cur["items"].append((k.strip(), val.strip(), com.strip()))
        else:
            cur["data"].append(s)
    return blocks


def bget(block, key, default=None):
    kl = key.lower()
    for k, v, _ in block["items"]:
        if k.lower() == kl:
            return v
    return default


def bgetall(block, key):
    kl = key.lower()
    return [v for k, v, _ in block["items"] if k.lower() == kl]


def bdict(block):
    return {k: v for k, v, _ in block["items"]}


def sections(blocks, name):
    return [b for b in blocks if b["section"].lower() == name.lower()]


class MDU:
    def __init__(self, path: Path):
        self.path = path
        self.dir = path.parent
        self.blocks = parse_ini(path)

    def __call__(self, section, key, default=None):
        for b in sections(self.blocks, section):
            v = bget(b, key)
            if v not in (None, ""):
                return v
        return default

    def any(self, key, default=None):
        for b in self.blocks:
            v = bget(b, key)
            if v not in (None, ""):
                return v
        return default


BRANCH_GUI_TYPE = {"0": "Channel", "1": "Sewer connection", "2": "Pipe"}


MATERIAL = {"0": "Unknown", "1": "Concrete", "2": "Cast iron", "3": "Stoneware", "4": "HDPE",
            "5": "Masonry", "6": "Sheet metal", "7": "Polyester", "8": "PVC", "9": "Steel"}


LINK_TYPE = {3: "embedded / internal", 4: "lateral", 5: "street inlet", 7: "roof gutter pipe"}


def _poly_area(xs, ys):
    a = 0.0
    for i in range(len(xs)):
        j = (i + 1) % len(xs)
        a += xs[i] * ys[j] - xs[j] * ys[i]
    return abs(a) / 2


def read_network(path: Path, warnings: list) -> dict:
    out = {"file": path.name, "size_kB": round(path.stat().st_size / 1024, 1), "_plot": {}}
    try:
        nc = open_nc(path)
    except (ValueError, OSError) as e:
        warnings.append(f"{path.name}: {e}")
        out["error"] = str(e)
        return out
    with nc:
        V = nc.vars
        out["conventions"] = nc.attrs.get("Conventions")
        out["created_by"] = nc.attrs.get("source")
        # projection
        for cand in ("projected_coordinate_system", "wgs84", "crs"):
            if cand in V:
                a = V[cand]["attrs"]
                epsg = a.get("epsg") or a.get("EPSG_code")
                out["projection"] = " ".join(str(x) for x in (a.get("name"), f"(EPSG {epsg})" if epsg else "") if x)
                break
        topo = {k: v for k, v in V.items() if v["attrs"].get("cf_role") == "mesh_topology"}
        net = m1 = m2 = None
        for k, v in topo.items():
            td = v["attrs"].get("topology_dimension")
            nco = str(v["attrs"].get("node_coordinates", ""))
            if td == 2:
                m2 = k
            elif td == 1 and ("branch" in nco or "offset" in nco):
                m1 = k
            elif td == 1:
                net = k

        def first(*names):
            return next((n for n in names if n in V), None)

        # ---- 1D network (branches)
        if net:
            bid = first(f"{net}_branch_id")
            ids = nc.strings(bid) if bid else []
            lname = first(f"{net}_branch_long_name")
            lnames = nc.strings(lname) if lname else [""] * len(ids)
            lens = nc.read(f"{net}_edge_length") if f"{net}_edge_length" in V else [None] * len(ids)
            btyp = nc.read(f"{net}_branch_type") if f"{net}_branch_type" in V else [None] * len(ids)
            nodes = nc.strings(f"{net}_node_id") if f"{net}_node_id" in V else []
            en = nc.read(f"{net}_edge_nodes") if f"{net}_edge_nodes" in V else []
            out["branches"] = [{"id": ids[i], "long_name": lnames[i] if i < len(lnames) else "",
                                "length_m": lens[i], "type_code": btyp[i],
                                "from": nodes[en[2 * i]] if en and nodes and en[2 * i] < len(nodes) else None,
                                "to": nodes[en[2 * i + 1]] if en and nodes and en[2 * i + 1] < len(nodes) else None}
                               for i in range(len(ids))]
            out["network_nodes"] = nodes
            if f"{net}_node_x" in V:
                out["_plot"]["net_nodes"] = list(zip(nc.read(f"{net}_node_x"), nc.read(f"{net}_node_y"), nodes))
            gx, gy = first(f"{net}_geom_x"), first(f"{net}_geom_y")
            cnt = first(f"{net}_geom_node_count")
            if gx and gy and cnt:
                xs, ys, cs = nc.read(gx), nc.read(gy), nc.read(cnt)
                lines, k = [], 0
                for i, c in enumerate(cs):
                    lines.append((ids[i] if i < len(ids) else str(i), xs[k:k + c], ys[k:k + c]))
                    k += c
                out["_plot"]["branches"] = lines
        # ---- 1D computational mesh
        if m1:
            out["mesh1d"] = {"name": m1, "nodes": nc.dim(V[f"{m1}_node_x"]["dims"][0]) if f"{m1}_node_x" in V else None,
                             "edges": nc.dim(V[f"{m1}_edge_nodes"]["dims"][0]) if f"{m1}_edge_nodes" in V else None}
            if f"{m1}_node_x" in V and f"{m1}_edge_nodes" in V:
                x, y = nc.read(f"{m1}_node_x"), nc.read(f"{m1}_node_y")
                en = nc.read(f"{m1}_edge_nodes")
                si = nc.att(f"{m1}_edge_nodes", "start_index", 0) or 0
                dx = [math.hypot(x[en[i + 1] - si] - x[en[i] - si], y[en[i + 1] - si] - y[en[i] - si])
                      for i in range(0, len(en), 2)]
                s = vstats(dx)
                if s:
                    out["mesh1d"]["dx_min_mean_max"] = f"{g(s['min'])} / {g(s['mean'])} / {g(s['max'])} m (straight-line)"
                out["_plot"]["mesh1d_xy"] = list(zip(x, y))
        # ---- 2D mesh
        if m2:
            nx_, ny_ = first(f"{m2}_node_x"), first(f"{m2}_node_y")
            fn = next((k for k, v in V.items() if v["attrs"].get("cf_role") == "face_node_connectivity"
                       and k.lower().startswith(m2.lower())), None)
            en = first(f"{m2}_edge_nodes")
            info = {"name": m2,
                    "nodes": nc.dim(V[nx_]["dims"][0]) if nx_ else None,
                    "edges": nc.dim(V[en]["dims"][0]) if en else None,
                    "faces": nc.dim(V[fn]["dims"][0]) if fn else None}
            if nx_ and fn:
                X, Y = nc.read(nx_), nc.read(ny_)
                F = nc.read(fn)
                mx = V[fn]["shape"][-1]
                si = nc.att(fn, "start_index", 0) or 0
                fill = nc.att(fn, "_FillValue", -999)
                polys, areas, nverts = [], [], {}
                for i in range(0, len(F), mx):
                    idx = [j - si for j in F[i:i + mx] if j != fill and j >= si]
                    px, py = [X[j] for j in idx], [Y[j] for j in idx]
                    polys.append((px, py))
                    areas.append(_poly_area(px, py))
                    nverts[len(idx)] = nverts.get(len(idx), 0) + 1
                sa = vstats(areas)
                info["cell_types"] = ", ".join(f"{n} {'triangles' if k == 3 else 'quads' if k == 4 else f'{k}-gons'}"
                                               for k, n in sorted(nverts.items()))
                if sa:
                    info["cell_area_min_mean_max"] = f"{g(sa['min'])} / {g(sa['mean'])} / {g(sa['max'])} m2"
                    info["approx_cell_size"] = f"~{g(math.sqrt(sa['mean']))} m"
                    info["total_area"] = f"{sum(areas) / 1e6:.4g} km2"
                info["extent_x"] = f"{g(min(X), 9)} .. {g(max(X), 9)}  ({g(max(X) - min(X))} m)"
                info["extent_y"] = f"{g(min(Y), 9)} .. {g(max(Y), 9)}  ({g(max(Y) - min(Y))} m)"
                out["_plot"]["faces"] = polys
                for zname, loc in ((f"{m2}_face_z", "face"), (f"{m2}_node_z", "node")):
                    if zname in V:
                        z = nc.read(zname)
                        s = vstats(z, fill=-999.0)
                        if s:
                            info[f"bed_level_{loc}_z"] = f"{g(s['min'])} .. {g(s['max'])} m (mean {g(s['mean'])}, {s['n']} values)"
                            if loc == "face":
                                out["_plot"]["face_z"] = [None if abs(v + 999) < 1e-9 else v for v in z]
                        else:
                            info[f"bed_level_{loc}_z"] = "not set (all fill values)"
                fx, fy = first(f"{m2}_face_x"), first(f"{m2}_face_y")
                if fx:
                    out["_plot"]["face_xy"] = list(zip(nc.read(fx), nc.read(fy)))
            out["mesh2d"] = info
        # ---- 1D2D links
        ck = next((k for k, v in V.items() if v["attrs"].get("cf_role") == "mesh_topology_contact"), None)
        if ck:
            n = V[ck]["shape"][0]
            types = nc.read(f"{ck}_contact_type") if f"{ck}_contact_type" in V else []
            cnt = {}
            for t in types:
                cnt[t] = cnt.get(t, 0) + 1
            out["links_1d2d"] = {"name": ck, "count": n, "contact": V[ck]["attrs"].get("contact"),
                                 "by_type": {f"type {t} {LINK_TYPE.get(t, '')}".strip(): c for t, c in sorted(cnt.items())}}
            out["_plot"]["links"] = nc.read(ck)
            out["_plot"]["links_start"] = nc.att(ck, "start_index", 0) or 0
        # ---- legacy (pre-UGRID) net files
        if not (net or m1 or m2) and "NetNode_x" in V:
            out["legacy_net"] = {"net_nodes": V["NetNode_x"]["shape"][0],
                                 "net_links": V["NetLink"]["shape"][0] if "NetLink" in V else None,
                                 "net_cells": V["NetElemNode"]["shape"][0] if "NetElemNode" in V else None}
    return out


_UNIT_SEC = {"seconds": 1, "second": 1, "s": 1, "minutes": 60, "minute": 60, "hours": 3600, "hour": 3600, "days": 86400}


def _time_unit(u):
    m = re.match(r"\s*(\w+)\s+since\s+(.+)", u or "", re.I)
    if not m:
        return 1, None
    return _UNIT_SEC.get(m.group(1).lower(), 1), m.group(2).strip()


def parse_pli(path: Path) -> list[dict]:
    lines = [ln.strip() for ln in read_text(path).splitlines() if ln.strip() and not ln.strip().startswith("*")]
    out, i = [], 0
    while i < len(lines):
        name = lines[i]
        i += 1
        if i >= len(lines):
            break
        hdr = lines[i].split()
        i += 1
        try:
            n = int(float(hdr[0]))
        except (ValueError, IndexError):
            break
        pts = []
        for ln in lines[i:i + n]:
            f = ln.split()
            pts.append((to_float(f[0]), to_float(f[1]), f[3] if len(f) > 3 else (f[2] if len(f) > 2 and to_float(f[2]) is None else "")))
        i += n
        out.append({"name": name, "points": n, "_xy": pts})
    return out


def parse_ext_old(path: Path, folder: Path) -> list[dict]:
    out, cur = [], None
    for ln in read_text(path).splitlines():
        s = ln.split("#", 1)[0].strip()
        if not s or s.startswith("*") or "=" not in s:
            continue
        k, v = (x.strip() for x in s.split("=", 1))
        if k.upper() == "QUANTITY":
            cur = {"QUANTITY": v}
            out.append(cur)
        elif cur is not None:
            cur[k.upper()] = v
    for e in out:
        if e.get("FILENAME"):
            e["exists"] = find_ci(folder, e["FILENAME"]) is not None
    return out


def parse_initial_fields(mdu: MDU, warnings: list) -> list[dict]:
    fn = mdu("geometry", "IniFieldFile")
    if not fn:
        return []
    fp = find_ci(mdu.dir, fn)
    if fp is None:
        warnings.append(f"IniFieldFile '{fn}' not found")
        return []
    out = []
    for b in parse_ini(fp):
        if b["section"].lower() not in ("initial", "parameter"):
            continue
        e = {"block": b["section"], "quantity": bget(b, "quantity"), "dataFile": bget(b, "dataFile"),
             "dataFileType": bget(b, "dataFileType"), "method": bget(b, "interpolationMethod"),
             "operand": bget(b, "operand"), "locationType": bget(b, "locationType"),
             "averaging": bget(b, "averagingType")}
        df = find_ci(fp.parent, e["dataFile"] or "")
        if e["dataFile"] and df is None:
            warnings.append(f"initial field dataFile '{e['dataFile']}' not found")
            e["content"] = "MISSING"
        elif df:
            t = (e["dataFileType"] or "").lower()
            try:
                if t == "1dfield":
                    bl = parse_ini(df)
                    gl = sections(bl, "Global")
                    e["content"] = "; ".join(f"global {bget(x, 'quantity')} = {bget(x, 'value')} {bget(x, 'unit', '')}" for x in gl)
                    nb = len(sections(bl, "Branch"))
                    if nb:
                        e["content"] += f"; {nb} branch-specific block(s)"
                elif t == "sample":
                    vals = []
                    for ln in read_text(df).splitlines():
                        f = ln.split()
                        if len(f) >= 3 and to_float(f[2]) is not None:
                            vals.append(to_float(f[2]))
                    s = vstats(vals)
                    e["content"] = f"{len(vals)} samples, value {g(s['min'])} .. {g(s['max'])} (mean {g(s['mean'])})" if s else "no samples"
                elif t == "polygon":
                    e["content"] = f"polygon, value = {bget(b, 'value')}"
                else:
                    e["content"] = f"{df.suffix} file, {df.stat().st_size / 1024:.1f} kB"
            except OSError as ex:
                e["content"] = str(ex)
        out.append(e)
    if (mdu("grw", "Infiltrationmodel") or "0") != "0" and any("infiltration" in (x["quantity"] or "").lower() and x["block"].lower() == "initial" for x in out):
        warnings.append("infiltrationcapacity is given as an [Initial] block; newer D-Flow FM versions expect [Parameter]")
    return out


_CODE_RE = re.compile(r"(?<![\d.\w-])(-?\d+)\s*[:=]\s*")


FALLBACK_CODES = {
    "layertype": {"1": "all sigma", "2": "all z"},
    "turbulencemodel": {"0": "none", "1": "constant", "2": "algebraic", "3": "k-epsilon", "4": "k-tau"},
    "temperature": {"0": "no", "1": "only transport", "3": "excess model of D3D", "5": "composite (ocean) model"},
    "salinity": {"0": "no", "1": "yes"},
    "wavemodelnr": {"0": "none", "1": "fetch/depth limited Hurdle-Stive", "2": "Young-Verhagen", "3": "SWAN"},
    "sedimentmodelnr": {"0": "no", "1": "Krone", "2": "SvR2007", "3": "E-H", "4": "MorphologyModule"},
    "uniffricttype": {"0": "Chezy", "1": "Manning", "2": "White-Colebrook", "3": "z0"},
    "timesteptype": {"0": "only transport", "1": "transport + velocity update", "2": "full implicit step-reduce",
                     "3": "step-Jacobi", "4": "explicit"},
    "bedlevtype": {"1": "at cell centers", "2": "at cell interfaces", "3": "at nodes, face levels mean of node values",
                   "4": "at nodes, face levels min. of node values", "5": "at nodes, face levels max. of node values",
                   "6": "at cell centers, cell levels mean of node values"},
    "idensform": {"0": "uniform", "1": "Eckart", "2": "UNESCO"},
}


def code_table(comment: str) -> dict:
    out, ms = {}, list(_CODE_RE.finditer(comment or ""))
    for i, m in enumerate(ms):
        end = ms[i + 1].start() if i + 1 < len(ms) else len(comment)
        txt = comment[m.end():end].strip().rstrip(",;").strip()
        while txt.endswith(")") and txt.count(")") > txt.count("("):
            txt = txt[:-1].rstrip(",; ")
        if txt:
            out.setdefault(m.group(1), txt)
    return out


def mdu_comment(mdu, key):
    kl = key.lower()
    for b in mdu.blocks:
        for k, v, c in b["items"]:
            if k.lower() == kl:
                return c
    return ""


def decode(mdu, section, key, default=None):
    """'3 = k-epsilon' using the MDU comment (or a fallback table)."""
    v = mdu(section, key, default)
    if v is None:
        return None
    tab = code_table(mdu_comment(mdu, key)) or FALLBACK_CODES.get(key.lower(), {})
    t = tab.get(str(v).strip())
    return f"{v} = {t}" if t else str(v)


def ref_datetime(mdu: MDU):
    rd = mdu("time", "RefDate")
    try:
        return _dt.datetime.strptime(rd.strip()[:8], "%Y%m%d")
    except (AttributeError, ValueError):
        return None


def tunit_sec(mdu: MDU) -> float:
    return {"S": 1, "M": 60, "H": 3600, "D": 86400}.get((mdu("time", "Tunit") or "S").strip().upper()[:1], 1)


def mdu_groups(mdu: MDU) -> dict:
    out = {}
    for b in mdu.blocks:
        items = [{"key": k, "value": v, "description": c} for k, v, c in b["items"]]
        out.setdefault(b["section"] or "(none)", []).extend(items)
    return out


FILE_KEYS = re.compile(r"(file|filenew)$", re.I)


def results(mdu: MDU, ti: dict, warnings: list) -> dict | None:
    stem = mdu.path.stem
    od = mdu("output", "OutputDir") or f"DFM_OUTPUT_{stem}"
    folder = None
    for base in (mdu.dir, mdu.dir.parent):          # Delta Shell keeps output/ next to input/
        f = find_ci(base, od)
        if f is not None and f.is_dir():
            folder = f
            break
    if folder is None:
        return None
    out = {"output_dir": str(folder), "_plot": {}}
    dia = next((p for p in sorted(folder.glob("*.dia")) if p.stem.lower().startswith(stem.lower())), None)
    if dia:
        txt = read_text(dia).splitlines()
        d = {"file": dia.name}
        warn = [ln.split(":", 1)[1].strip() for ln in txt if ln.startswith("** WARNING")]
        err = [ln.split(":", 1)[1].strip() for ln in txt if ln.startswith("** ERROR")]
        d["warnings"] = list(dict.fromkeys(warn))
        d["errors"] = list(dict.fromkeys(err))
        counts = {}
        for ln in txt:
            m = re.match(r"\*\* INFO\s*:\s*(nr of [^(]+?)\s*\(.*?\)\s*:\s*(\S+)", ln)
            if m:
                counts[m.group(1).strip()] = m.group(2)
            m = re.match(r"\*\* INFO\s*:\s*(Computation (started|finished)\s+at|total computation time\s*\(s\)|"
                         r"simulation period\s*\(h\)|my model area\s*\(m2\)|my model volume\s*\(m3\)|"
                         r"time steps \(\+ plots\)\s*\(s\))\s*:\s*(.+)", ln, re.I)
            if m:
                d[re.sub(r"\s+", " ", m.group(1)).strip()] = m.group(3).strip()
        d["counts"] = counts
        d["finished"] = any("Computation finished" in ln for ln in txt)
        out["dia"] = d
        if err:
            warnings.append(f"{dia.name}: {len(err)} ERROR message(s)")
    for kind in ("map", "his"):
        f = next((p for p in sorted(folder.glob(f"*_{kind}.nc")) if p.name.lower().startswith(stem.lower())), None) \
            or next(iter(sorted(folder.glob(f"*_{kind}.nc"))), None)
        if f is None:
            continue
        e = {"file": f.name, "size_MB": round(f.stat().st_size / 1e6, 2)}
        try:
            with open_nc(f) as nc:
                tv = "time" if "time" in nc.vars else None
                if tv:
                    t = nc.read(tv)
                    fac, ref = _time_unit(nc.att(tv, "units", ""))
                    e["time_steps_written"] = len(t)
                    if t:
                        e["time_written"] = f"{fmt_dur(t[0] * fac)} .. {fmt_dur(t[-1] * fac)}"
                        e["completed"] = "yes" if t[-1] * fac >= ti["_t1"] - 1e-6 else \
                            f"NO - last output at {100 * (t[-1] * fac - ti['_t0']) / max(1e-9, ti['_t1'] - ti['_t0']):.0f} % of the period"
                e["variables"] = [k for k, v in nc.vars.items() if v["record"] and k not in ("time", "timestep", "time_bds")]
                if kind == "map":
                    _map_stats(nc, e, out["_plot"], mdu)
                else:
                    _his_stats(nc, e)
        except (ValueError, OSError, KeyError, struct.error) as ex:
            e["error"] = str(ex)
        out[kind] = e
    return out


def _map_stats(nc, e, plot, mdu):
    thr = to_float(mdu("output", "Wrimap_wet_waterdepth_threshold"), 0.01) or 0.01
    for var in [k for k in nc.vars if k.lower().endswith("_waterdepth") and "on_ground" not in k.lower()]:
        n = nc.vars[var]["count"] if "count" in nc.vars[var] else nc.shape(var)[-1]
        nt = nc.numrecs
        if nt * n > 3e7:
            e[f"{var} note"] = "too large - max over last time step only"
            mx = nc.read(var, rec=-1)
        else:
            mx = None
            for r in range(nt):
                a = nc.read(var, rec=r)
                mx = a if mx is None else [max(x, y) for x, y in zip(mx, a)]
        if not mx:
            continue
        s = vstats(mx)
        wet = sum(1 for v in mx if v > thr)
        e[f"max {var}"] = f"{g(s['max'])} m; {wet} of {len(mx)} elements wet (> {g(thr)} m) at some time"
        if "mesh2d" in var.lower() and not var.lower().startswith("mesh1d"):
            plot["max_depth_2d"] = mx
            ba = next((k for k in nc.vars if k.lower() == var.lower().replace("_waterdepth", "_flowelem_ba")), None)
            if ba:
                area = nc.read(ba)
                e["max flooded 2D area"] = f"{sum(a for a, d in zip(area, mx) if d > thr) / 1e6:.4g} km2 (depth > {g(thr)} m)"
    for var in [k for k in nc.vars if k.lower().endswith("flowelem_bl")]:
        s = vstats(nc.read(var))
        if s:
            e[f"{var} (bed level used)"] = f"{g(s['min'])} .. {g(s['max'])} m"
            if var.lower().startswith("mesh2d"):
                plot["bl_2d"] = nc.read(var)


def _his_stats(nc, e):
    for nm in ("station_name", "station_id", "cross_section_name", "weirgen_id", "pump_id", "culvert_id"):
        if nm in nc.vars:
            try:
                e[nm] = nc.strings(nm)
            except Exception:
                pass
    for nm in ("water_balance_volume_error", "water_balance_total_volume", "water_balance_boundaries_total",
               "water_balance_exchange_with_1D_total", "water_balance_precipitation_total"):
        if nm in nc.vars:
            v = nc.read(nm, rec=-1)
            if v:
                e[nm.replace("water_balance_", "balance: ") + " (end)"] = f"{g(v[0])} m3"


def find_mdus(arg: str) -> list[Path]:
    p = Path(arg).expanduser().resolve()
    if not p.exists():
        sys.exit(f"Path not found: {p}")
    if p.is_file() and p.suffix.lower() == ".mdu":
        return [p]
    if p.is_file() and p.suffix.lower() == ".dsproj":
        p = p.with_name(p.name + "_data")
        if not p.is_dir():
            sys.exit(f"{p} not found next to the .dsproj")
    root = p if p.is_dir() else p.parent
    out = []
    for f in sorted(root.rglob("*.mdu")):
        parts = [x.lower() for x in f.relative_to(root).parts]
        if any(x.startswith("dfm_output") or x == "output" for x in parts[:-1]) or "_interpreted" in f.name.lower():
            continue
        out.append(f)
    if not out:
        sys.exit(f"No .mdu file found under {root}")
    return out


def text_table(headers, rows, indent=2):
    if not rows:
        return
    w = [max(len(str(h)), *(len(str(r[i])) for r in rows)) for i, h in enumerate(headers)]
    w = [min(x, 40) for x in w]
    pad = " " * indent
    print(pad + "  ".join(str(h)[:w[i]].ljust(w[i]) for i, h in enumerate(headers)))
    print(pad + "  ".join("-" * x for x in w))
    for r in rows:
        print(pad + "  ".join(str(c if c is not None else "")[:w[i]].ljust(w[i]) for i, c in enumerate(r)))


RAMP_BED = [(0, (44, 123, 182)), (0.25, (171, 217, 233)), (0.5, (255, 255, 191)), (0.75, (253, 174, 97)), (1, (166, 97, 26))]


RAMP_DEPTH = [(0, (237, 248, 251)), (0.3, (158, 202, 225)), (0.65, (49, 130, 189)), (1, (8, 48, 107))]


def ramp(t, stops):
    t = max(0.0, min(1.0, t))
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t <= t1:
            f = (t - t0) / ((t1 - t0) or 1)
            return "#%02x%02x%02x" % tuple(round(a + (b - a) * f) for a, b in zip(c0, c1))
    return "#%02x%02x%02x" % stops[-1][1]


def svg_lines(series, xlabel="", ylabel="", w=820, h=220) -> str:
    allx = [v for s in series for v in s[1]]
    ally = [v for s in series for v in s[2]]
    if not allx:
        return ""
    x0, x1, y0, y1 = min(allx), max(allx), min(ally), max(ally)
    pad = (y1 - y0) * 0.08 or max(0.5, abs(y1) * 0.1)
    y0, y1 = y0 - pad, y1 + pad
    x1 = x1 if x1 > x0 else x0 + 1
    L, R, T, B = 58, 14, 14, 38
    X = lambda v: L + (v - x0) / (x1 - x0) * (w - L - R)
    Y = lambda v: T + (y1 - v) / (y1 - y0) * (h - T - B)
    o = [f'<svg class="plot" viewBox="0 0 {w} {h}">']
    for t in _ticks(y0, y1):
        o.append(f'<line x1="{L}" x2="{w - R}" y1="{Y(t):.1f}" y2="{Y(t):.1f}" class="gl"/>'
                 f'<text x="{L - 6}" y="{Y(t) + 4:.1f}" text-anchor="end">{g(t)}</text>')
    for t in _ticks(x0, x1, 8):
        o.append(f'<text x="{X(t):.1f}" y="{h - B + 16}" text-anchor="middle">{g(t)}</text>')
    for i, (lab, xs, ys) in enumerate(series):
        col = PALETTE[i % len(PALETTE)]
        pts = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in zip(xs, ys))
        o.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="1.8"><title>{_h.escape(lab)}</title></polyline>')
    o.append(f'<rect x="{L}" y="{T}" width="{w - L - R}" height="{h - T - B}" class="frame"/>')
    o.append(f'<text x="{(L + w - R) / 2}" y="{h - 4}" text-anchor="middle" class="lab">{_h.escape(xlabel)}</text>')
    o.append(f'<text transform="translate(13,{(T + h - B) / 2}) rotate(-90)" text-anchor="middle" class="lab">{_h.escape(ylabel)}</text>')
    if 1 < len(series) <= 8:
        lx = L + 10
        for i, (lab, _, _) in enumerate(series):
            o.append(f'<rect x="{lx}" y="{T + 6}" width="12" height="3" fill="{PALETTE[i % len(PALETTE)]}"/>'
                     f'<text x="{lx + 16}" y="{T + 11}">{_h.escape(lab)}</text>')
            lx += 24 + 7 * len(lab)
    o.append("</svg>")
    return "".join(o)
