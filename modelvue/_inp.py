"""Code shared by the EPANET and SWMM readers (both use the EPA .inp text format)."""
from __future__ import annotations


def num(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def table(inp, sec):
    return inp["sections"].get(sec, [])
