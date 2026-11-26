"""HTML page, SVG plots, map scaling and a small shapefile reader shared by rasvue and hmsvue.

Standard library only. The page has the same layout as the other modelvue reports.
"""
from __future__ import annotations

import datetime as _dt
import html as _h
import math
import struct
from pathlib import Path

from ._common import _e, _ticks

PALETTE = ["#1f5fa8", "#d9622b", "#2e7d32", "#8e5bb5", "#c62828", "#b8860b", "#16a2a2", "#7f8c8d"]


def g(v, nd=4):
    """Compact number text: integers without decimals, others with ``nd`` significant digits."""
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    try:
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return ""
        if float(v).is_integer() and abs(v) < 1e12:
            return str(int(v))
        return f"{v:.{nd}g}"
    except (TypeError, ValueError):
        return str(v)


def num(v, default=None):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return default


def read_text(path: Path) -> str:
    data = Path(path).read_bytes()
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("latin-1")


def find_ci(folder: Path, name: str) -> Path | None:
    """File ``name`` in ``folder``, ignoring letter case (HEC files mix EX1.F01 / ex1.f01)."""
    if not name:
        return None
    p = folder / name
    if p.exists():
        return p
    low = name.lower()
    try:
        for q in folder.iterdir():
            if q.name.lower() == low:
                return q
    except OSError:
        pass
    return None


# --------------------------------------------------------------------------- #
# Terminal helpers
# --------------------------------------------------------------------------- #
def text_table(headers, rows, indent=2, maxw=22, limit=None):
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


# --------------------------------------------------------------------------- #
# Shapefile (.shp) reader - points, polylines and polygons, enough to draw a map
# --------------------------------------------------------------------------- #
SHP_KIND = {1: "point", 11: "point", 21: "point", 8: "point", 18: "point", 28: "point",
            3: "line", 13: "line", 23: "line", 5: "polygon", 15: "polygon", 25: "polygon"}


def shp_kind(path: Path) -> str | None:
    """'point', 'line' or 'polygon' from the shapefile header."""
    with open(path, "rb") as f:
        hdr = f.read(36)
    return SHP_KIND.get(struct.unpack("<i", hdr[32:36])[0]) if len(hdr) == 36 else None


def read_shp(path: Path, max_shapes=20000) -> list[list[list[tuple[float, float]]]]:
    """Return a list of shapes; each shape is a list of parts; each part a list of (x, y)."""
    out = []
    with open(path, "rb") as f:
        hdr = f.read(100)
        if len(hdr) < 100 or struct.unpack(">i", hdr[:4])[0] != 9994:
            raise ValueError("not a shapefile")
        while len(out) < max_shapes:
            rh = f.read(8)
            if len(rh) < 8:
                break
            _, clen = struct.unpack(">2i", rh)
            rec = f.read(clen * 2)
            if len(rec) < 4:
                break
            st = struct.unpack("<i", rec[:4])[0]
            if st in (1, 11, 21):
                x, y = struct.unpack("<2d", rec[4:20])
                out.append([[(x, y)]])
            elif st in (3, 5, 13, 15, 23, 25):
                nparts, npts = struct.unpack("<2i", rec[36:44])
                parts = list(struct.unpack(f"<{nparts}i", rec[44:44 + 4 * nparts])) + [npts]
                base = 44 + 4 * nparts
                xy = struct.unpack(f"<{2 * npts}d", rec[base:base + 16 * npts])
                pts = list(zip(xy[0::2], xy[1::2]))
                out.append([pts[parts[i]:parts[i + 1]] for i in range(nparts)])
            elif st == 8 or st in (18, 28):
                npts = struct.unpack("<i", rec[36:40])[0]
                xy = struct.unpack(f"<{2 * npts}d", rec[40:40 + 16 * npts])
                out.append([[p] for p in zip(xy[0::2], xy[1::2])])
    return out


# --------------------------------------------------------------------------- #
# SVG
# --------------------------------------------------------------------------- #
class MapFrame:
    """Scales world coordinates (y up) into an SVG of at most ``w`` x ``hmax`` pixels."""

    def __init__(self, pts, w=920, hmax=720, pad=0.04):
        pts = [p for p in pts if p and p[0] is not None and p[1] is not None and math.isfinite(p[0]) and math.isfinite(p[1])]
        self.ok = bool(pts)
        if not pts:
            self.w = self.h = 0
            return
        x0, x1 = min(p[0] for p in pts), max(p[0] for p in pts)
        y0, y1 = min(p[1] for p in pts), max(p[1] for p in pts)
        span = max(x1 - x0, y1 - y0) or 1.0
        x0, x1, y0, y1 = x0 - span * pad, x1 + span * pad, y0 - span * pad, y1 + span * pad
        sc = w / ((x1 - x0) or 1)
        if (y1 - y0) * sc > hmax:
            sc = hmax / ((y1 - y0) or 1)
        self.x0, self.y1, self.sc = x0, y1, sc
        self.w, self.h = (x1 - x0) * sc, (y1 - y0) * sc

    def X(self, v):
        return (v - self.x0) * self.sc

    def Y(self, v):
        return (self.y1 - v) * self.sc

    def pts(self, P):
        return " ".join(f"{self.X(x):.1f},{self.Y(y):.1f}" for x, y in P)

    def svg(self, parts, cls="map"):
        # max-width keeps tall, narrow maps from being stretched to the full page width
        return (f'<svg class="{cls}" viewBox="0 0 {self.w:.0f} {self.h:.0f}" style="max-width:{max(self.w, 320):.0f}px">'
                + "".join(parts) + "</svg>")


def frame(w, h, L, R, T, B, x0, x1, y0, y1, xlabel, ylabel, xfmt=None, flip=False):
    """Axes for an x/y plot; returns (X, Y, svg parts)."""
    if x1 == x0:
        x1 = x0 + 1
    if y1 == y0:
        y1 = y0 + 1
    X = lambda v: L + (v - x0) / (x1 - x0) * (w - L - R)
    Y = (lambda v: T + (v - y0) / (y1 - y0) * (h - T - B)) if flip else (lambda v: T + (y1 - v) / (y1 - y0) * (h - T - B))
    o = []
    for t in _ticks(y0, y1):
        o.append(f'<line x1="{L}" x2="{w - R}" y1="{Y(t):.1f}" y2="{Y(t):.1f}" class="gl"/>'
                 f'<text x="{L - 6}" y="{Y(t) + 4:.1f}" text-anchor="end">{g(t)}</text>')
    for t in _ticks(x0, x1, 7):
        o.append(f'<text x="{X(t):.1f}" y="{h - B + 15}" text-anchor="middle">{_e(xfmt(t) if xfmt else g(t))}</text>')
    o.append(f'<rect x="{L}" y="{T}" width="{w - L - R}" height="{h - T - B}" class="frame"/>')
    o.append(f'<text x="{(L + w - R) / 2}" y="{h - 3}" text-anchor="middle" class="lab">{_e(xlabel)}</text>')
    o.append(f'<text transform="translate(13,{(T + h - B) / 2}) rotate(-90)" text-anchor="middle" class="lab">{_e(ylabel)}</text>')
    return X, Y, o


def svg_lines(series, xlabel="", ylabel="", w=820, h=230, xfmt=None, step=False, ymin=None) -> str:
    """series: [(label, colour, xs, ys)] -> line chart; ``step`` draws a stepped (bar-like) line."""
    series = [s for s in series if s[2] and s[3]]
    if not series:
        return ""
    xs = [x for s in series for x in s[2]]
    ys = [y for s in series for y in s[3] if y is not None]
    if not ys:
        return ""
    y0 = min(ys) if ymin is None else min(ymin, min(ys))
    y1 = max(ys)
    pad = (y1 - y0) * 0.06 or abs(y1) * 0.1 or 1
    X, Y, o = frame(w, h, 62, 12, 14, 34, min(xs), max(xs), y0 if ymin is not None else y0 - pad, y1 + pad, xlabel, ylabel, xfmt)
    for lab, col, sx, sy in series:
        P = []
        for i, (a, b) in enumerate(zip(sx, sy)):
            if b is None:
                continue
            if step and i:
                P.append(f"{X(a):.1f},{Y(sy[i - 1] if sy[i - 1] is not None else b):.1f}")
            P.append(f"{X(a):.1f},{Y(b):.1f}")
        o.append(f'<polyline points="{" ".join(P)}" fill="none" stroke="{col}" stroke-width="1.8"><title>{_e(lab)}</title></polyline>')
    if len(series) > 1 or series[0][0]:
        lx = 72
        for lab, col, _, _ in series[:8]:
            o.append(f'<rect x="{lx}" y="20" width="12" height="3" fill="{col}"/><text x="{lx + 16}" y="25">{_e(lab)}</text>')
            lx += 30 + 6.2 * len(str(lab))
    return f'<svg class="plot" viewBox="0 0 {w} {h}">' + "".join(o) + "</svg>"


def svg_bars(labels, values, unit="", colour="#1565c0", w=820, maxbars=40) -> str:
    """Horizontal bar chart, one bar per label."""
    rows = [(lab, v) for lab, v in zip(labels, values) if v is not None][:maxbars]
    if not rows:
        return ""
    vmax = max(v for _, v in rows) or 1
    L, bh = min(240, 12 + 7 * max(len(str(lab)) for lab, _ in rows)), 18
    h = 16 + bh * len(rows)
    o = []
    for k, (lab, v) in enumerate(rows):
        y = 8 + k * bh
        bw = max(1.0, (w - L - 90) * max(v, 0) / vmax)
        o.append(f'<text x="{L - 6}" y="{y + 12}" text-anchor="end" class="lab">{_e(lab)}</text>'
                 f'<rect x="{L}" y="{y + 2}" width="{bw:.1f}" height="{bh - 5}" fill="{colour}"><title>{_e(lab)}: {g(v)} {_e(unit)}</title></rect>'
                 f'<text x="{L + bw + 5:.1f}" y="{y + 12}">{g(v)} {_e(unit)}</text>')
    return f'<svg class="plot" viewBox="0 0 {w} {h}">' + "".join(o) + "</svg>"


def legend(items) -> str:
    """items: [(css style for the swatch, label, kind)] with kind 'ln' (line), 'bx' (box) or 'dot'."""
    out = ['<div class="legend">']
    for style, label, kind in items:
        out.append(f'<span><i class="{kind}" style="{style}"></i>{_e(label)}</span>')
    return "".join(out) + "</div>"


# --------------------------------------------------------------------------- #
# HTML page
# --------------------------------------------------------------------------- #
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<script>(function(){var t;try{t=localStorage.getItem('modelvue-theme')}catch(e){}if(t!=='dark')t='light';document.documentElement.dataset.theme=t;})();</script>
<title>@@APP@@ · @@TITLE@@</title>
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
.tw{overflow-x:auto} table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
td{white-space:nowrap} td.wrap,td .wrap{white-space:normal}
thead th{background:var(--accs);color:var(--acc);cursor:pointer;user-select:none;white-space:nowrap}
thead th:hover{text-decoration:underline} tbody tr:hover{background:var(--accs)}
table.kv{width:auto;min-width:50%} table.kv th{color:var(--muted);font-weight:500;background:none;cursor:auto;padding-right:24px;white-space:nowrap}
table.kv td{white-space:normal;overflow-wrap:anywhere}
td.msg{white-space:normal}
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
svg.map{width:100%;max-width:1000px;height:auto;display:block;margin:8px 0;background:var(--bg);border:1px solid var(--line);border-radius:6px}
svg.map text{font-size:11px;fill:var(--ink)}
svg.plot{width:100%;max-width:980px;height:auto;display:block;margin:8px 0;font-size:11px}
svg.plot text{fill:var(--muted)} svg.plot .lab{fill:var(--ink)} svg.plot .gl{stroke:var(--line)} svg.plot .frame{fill:none;stroke:var(--line)}
.legend{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:12px;color:var(--muted);margin-bottom:8px}
.legend i{display:inline-block;vertical-align:middle;margin-right:5px}
.legend .ln{width:22px;height:0;border-top:3px solid} .legend .bx{width:11px;height:11px} .legend .dot{width:9px;height:9px;border-radius:50%}
.desc{white-space:pre-wrap;color:var(--muted);font-size:13px;max-height:260px;overflow:auto;border:1px solid var(--line);border-radius:6px;padding:8px 12px;margin:6px 0}
pre.log{white-space:pre-wrap;font-size:12px;max-height:320px;overflow:auto;border:1px solid var(--line);border-radius:6px;padding:8px 12px;margin:6px 0}
ul.warns{margin:6px 0;padding-left:20px;color:var(--warn)}
@media (max-width:1000px){.wrap{grid-template-columns:1fr} nav{display:none}}
footer{color:var(--muted);font-size:12px;text-align:center;padding:10px}
@media print{nav,header{display:none} .wrap{display:block} details>*{display:block}}
</style></head>
<body>
<header><h1>@@APP@@ · @@TITLE@@</h1><span>@@SUB@@</span><button class="theme" type="button" title="Switch between the light and dark theme">☾ Dark</button></header>
<div class="wrap"><nav><div class="tools"><button data-x="1">Expand all</button><button data-x="0">Collapse all</button></div>@@NAV@@</nav>
<main>@@BODY@@
<footer>Generated by modelvue.@@TOOL@@ on @@GEN@@ from <code>@@PATH@@</code></footer></main></div>
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
document.querySelectorAll('input.pfilter').forEach(inp=>inp.addEventListener('input',()=>{
  const q=inp.value.toLowerCase(), box=inp.closest('.cb');
  box.querySelectorAll('tbody tr').forEach(r=>r.style.display=r.innerText.toLowerCase().includes(q)?'':'none');}));
const links=[...document.querySelectorAll('nav a')];
const obs=new IntersectionObserver(es=>es.forEach(e=>{if(e.isIntersecting){links.forEach(l=>l.classList.toggle('on',l.getAttribute('href')==='#'+e.target.id));}}),{rootMargin:'-30% 0px -60% 0px'});
document.querySelectorAll('main section, main details.blk[id]').forEach(s=>obs.observe(s));
const themeBtn=document.querySelector('header .theme');
const themeLabel=()=>{themeBtn.textContent=document.documentElement.dataset.theme==='dark'?'☀ Light':'☾ Dark';};
themeLabel();
themeBtn.addEventListener('click',()=>{const t=document.documentElement.dataset.theme==='dark'?'light':'dark';
  document.documentElement.dataset.theme=t; themeLabel(); try{localStorage.setItem('modelvue-theme',t)}catch(e){}});
</script></body></html>"""


def page(app, tool, title, sub, nav, body, path, extra_css="") -> str:
    rep = {"@@APP@@": _e(app), "@@TOOL@@": tool, "@@TITLE@@": _e(title), "@@SUB@@": _e(sub), "@@NAV@@": "".join(nav),
           "@@BODY@@": "".join(body), "@@GEN@@": _dt.datetime.now().strftime("%Y-%m-%d %H:%M"), "@@PATH@@": _e(path)}
    out = HTML_TEMPLATE.replace("</style>", extra_css + "</style>")
    for k, v in rep.items():
        out = out.replace(k, v)
    return out


def write_page(text, out: Path, open_browser=True) -> Path:
    out = Path(out)
    out.write_text(text, encoding="utf-8")
    print(f"Webpage written to {out.resolve()}")
    if open_browser:
        import webbrowser
        webbrowser.open(out.resolve().as_uri())
    return out


def esc(x) -> str:
    return _h.escape("" if x is None else str(x))
