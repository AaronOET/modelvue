"""Work out which reader (engine) handles a given file or folder."""
from __future__ import annotations

import re
import struct
from pathlib import Path

#: tool (command / module) name -> one-line description
ENGINES = {
    "dfmvue": "Delft3D FM / D-HYDRO 1D2D model (.mdu, .dsproj)",
    "dfmvue23": "Delft3D FM 2D/3D model, incl. coupled D-Waves (.mdu, .dsproj)",
    "epanetvue": "EPANET 2.x water distribution model (.inp, .out)",
    "swmmvue": "EPA SWMM 5 drainage model (.inp, .out)",
    "flovue": "MODFLOW 6 / 2005 / NWT / USG, MT3D, SEAWAT (mfsim.nam, .nam) - needs flopy",
    "xbvue": "XBeach model or case set (params.txt, .config)",
    "rasvue": "HEC-RAS project (.prj, plus .g/.f/.u/.p files)",
    "hmsvue": "HEC-HMS project (.hms, plus .basin/.met/.control/.run files)",
}

EPANET_OUT_MAGIC = 516114521
SWMM_OUT_MAGIC = 516114522

_SWMM_SECTIONS = {"SUBCATCHMENTS", "CONDUITS", "OUTFALLS", "RAINGAGES", "XSECTIONS", "SUBAREAS", "INFILTRATION",
                  "STORAGE", "ORIFICES", "WEIRS", "OUTLETS", "DIVIDERS", "TRANSECTS", "LID_CONTROLS"}
_EPANET_SECTIONS = {"PIPES", "RESERVOIRS", "TANKS", "VALVES", "DEMANDS", "EMITTERS", "ENERGY", "REACTIONS", "MIXING"}
_ONED_KEYS = ("crossdeffile", "crosslocfile", "storagenodefile", "branchfile")
_SKIP_DIRS = {".git", "__pycache__", "output", "dflowfmoutput"}


def _inp_kind(path: Path) -> str | None:
    try:
        text = path.read_text(encoding="latin-1", errors="replace")
    except OSError:
        return None
    secs = {s.upper() for s in re.findall(r"^\s*\[([A-Za-z_]+)\]", text, re.M)}
    s, e = len(secs & _SWMM_SECTIONS), len(secs & _EPANET_SECTIONS)
    if s > e:
        return "swmmvue"
    if e > s:
        return "epanetvue"
    if "JUNCTIONS" in secs:  # both have JUNCTIONS; SWMM junctions have an elevation + max depth, EPANET has demand
        return "epanetvue"
    return None


def _out_kind(path: Path) -> str | None:
    try:
        with open(path, "rb") as f:
            head = f.read(4)
    except OSError:
        return None
    if len(head) < 4:
        return None
    magic = struct.unpack("<i", head)[0]
    return {EPANET_OUT_MAGIC: "epanetvue", SWMM_OUT_MAGIC: "swmmvue"}.get(magic)


def _mdu_kind(path: Path) -> str:
    """1D2D (dfmvue) when the MDU names 1D input files or its net file holds a 1D network, else 2D/3D (dfmvue23)."""
    try:
        text = path.read_text(encoding="latin-1", errors="replace")
    except OSError:
        return "dfmvue23"
    netfile = None
    for line in text.splitlines():
        k, _, v = line.partition("=")
        k, v = k.strip().lower(), v.split("#")[0].strip()
        if k in _ONED_KEYS and v:
            return "dfmvue"
        if k == "netfile" and v:
            netfile = v
    if netfile:
        nf = path.parent / netfile
        if nf.is_file():
            try:
                with open(nf, "rb") as f:
                    head = f.read(1 << 20)
                if re.search(rb"network1d|mesh1d|Network_branch|network_edge_nodes", head, re.I):
                    return "dfmvue"
            except OSError:
                pass
    return "dfmvue23"


def _is_ras_prj(p: Path) -> bool:
    """HEC-RAS project files start with 'Proj Title='; ESRI projection files (.prj too) do not."""
    try:
        with open(p, "rb") as f:
            return f.read(200).decode("latin-1", "replace").lstrip().startswith("Proj Title=")
    except OSError:
        return False


def _file_kind(p: Path) -> str | None:
    name, suf = p.name.lower(), p.suffix.lower()
    if suf == ".prj":
        return "rasvue" if _is_ras_prj(p) else None
    if suf == ".rasmap":
        return "rasvue"
    if suf in (".hms", ".basin", ".met", ".control"):
        return "hmsvue"
    if suf == ".mdu":
        return _mdu_kind(p)
    if suf == ".dsproj":
        mdus = sorted(p.parent.glob(p.stem + ".dsproj_data/**/*.mdu"))
        return _mdu_kind(mdus[0]) if mdus else "dfmvue23"
    if suf == ".inp":
        return _inp_kind(p)
    if suf in (".out", ".bin"):
        return _out_kind(p)
    if name == "mfsim.nam" or suf in (".nam", ".mfn"):
        return "flovue"
    if name == "params.txt" or suf == ".config":
        return "xbvue"
    return None


def _walk(folder: Path, depth: int):
    try:
        entries = sorted(folder.iterdir())
    except OSError:
        return
    for e in entries:
        if e.is_file():
            yield e
    if depth > 0:
        for e in entries:
            if e.is_dir() and e.name.lower() not in _SKIP_DIRS:
                yield from _walk(e, depth - 1)


def candidates(path, depth: int = 5) -> dict[str, list[Path]]:
    """Every engine that recognises something at ``path``, with the files that triggered it."""
    p = Path(path)
    found: dict[str, list[Path]] = {}
    files = [p] if p.is_file() else list(_walk(p, depth)) if p.is_dir() else []
    for f in files:
        k = _file_kind(f)
        if k:
            found.setdefault(k, []).append(f)
    return found


def detect(path, depth: int = 5) -> str:
    """Return the tool (engine) name for a model file or folder.

    Raises ``ValueError`` when nothing is recognised, or when a folder holds models of
    more than one kind (pass the engine explicitly then).
    """
    p = Path(path)
    if not p.exists():
        raise ValueError(f"{p} does not exist")
    found = candidates(p, depth)
    if not found:
        raise ValueError(f"no supported model found at {p}")
    if len(found) > 1:
        detail = "; ".join(f"{k}: {v[0]}" + (f" (+{len(v) - 1})" if len(v) > 1 else "") for k, v in found.items())
        raise ValueError(f"{p} holds more than one kind of model ({detail}); name the engine, e.g. "
                         f"'{next(iter(found))} {p}'")
    return next(iter(found))
