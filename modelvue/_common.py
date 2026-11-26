"""Helpers shared by several modelvue readers: number/text formatting, terminal tables, HTML building blocks."""
from __future__ import annotations
import html as _h
import math


def h1(t):
    print("\n" + "=" * 78 + f"\n {t}\n" + "=" * 78)


def h2(t):
    print(f"\n--- {t} " + "-" * max(0, 72 - len(t)))


def row(label, value, indent=2):
    if value in (None, "", [], {}):
        return
    if isinstance(value, dict):
        value = ", ".join(f"{k}: {v}" for k, v in value.items())
    elif isinstance(value, list):
        value = ", ".join(map(str, value))
    print(f"{' ' * indent}{label:<30}: {value}")


def wrap_pairs(pairs, indent=4, width=78):
    line = " " * indent
    for k, v in pairs:
        item = f"{k} = {v}"
        if len(line) + len(item) + 3 > width and line.strip():
            print(line.rstrip())
            line = " " * indent
        line += item + "   "
    if line.strip():
        print(line.rstrip())


PALETTE = ["#2f80c0", "#d9622b", "#3a9d5d", "#8e5bb5", "#c0392b", "#b8860b", "#16a2a2", "#7f8c8d"]


def _ticks(lo, hi, n=5):
    if hi == lo:
        hi = lo + 1
    span = hi - lo
    step = 10 ** math.floor(math.log10(span / n))
    for k in (1, 2, 2.5, 5, 10):
        if span / (step * k) <= n:
            step *= k
            break
    t, out = math.ceil(lo / step) * step, []
    while t <= hi + 1e-9 * span:
        out.append(round(t, 10))
        t += step
    return out


def _e(x) -> str:
    if isinstance(x, dict):
        return ", ".join(f"{_h.escape(str(k))}: {_h.escape(str(v))}" for k, v in x.items())
    if isinstance(x, (list, tuple)):
        return ", ".join(_h.escape(str(i)) for i in x)
    return _h.escape("" if x is None else str(x))


def _kv(pairs) -> str:
    rows = "".join(f"<tr><th>{_e(k)}</th><td>{_e(v)}</td></tr>" for k, v in pairs if v not in (None, "", [], {}))
    return f'<div class="tw"><table class="kv">{rows}</table></div>' if rows else ""


def _table(headers, rows, sortable=True) -> str:
    if not rows:
        return ""
    th = "".join(f"<th>{_e(h)}</th>" for h in headers)
    tr = "".join("<tr>" + "".join(f'<td class="{c[1]}">{c[0]}</td>' if isinstance(c, tuple) else f"<td>{c}</td>" for c in r)
                 + "</tr>" for r in rows)
    return f'<div class="tw"><table class="{"sortable" if sortable else ""}"><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></div>'


def _blk(title, body, note="", open_=True, bid=None) -> str:
    if not body:
        return ""
    n = f"<small>{_e(note)}</small>" if note else ""
    i = f' id="{bid}"' if bid else ""
    return f'<details class="blk"{i}{" open" if open_ else ""}><summary>{title}{n}</summary><div class="cb">{body}</div></details>'


def _two(pairs):
    h = (len(pairs) + 1) // 2
    return f'<div class="grid2"><div>{_kv(pairs[:h])}</div><div>{_kv(pairs[h:])}</div></div>'


def strip_private(o):
    if isinstance(o, dict):
        return {k: strip_private(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, list):
        return [strip_private(v) for v in o]
    return o
