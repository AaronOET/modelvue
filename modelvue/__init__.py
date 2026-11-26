"""
MODELVUE - Read hydraulic / hydrological model setups and show their configuration.

Tools (each is a module and a command of the same name):

    dfmvue      Delft3D FM / D-HYDRO 1D2D model
    dfmvue23    Delft3D FM 2D / 3D model, incl. D-Waves
    epanetvue   EPANET 2.x model
    swmmvue     EPA SWMM 5 model
    flovue      MODFLOW model, read with FloPy (needs flopy)
    xbvue       XBeach model or case set
    rasvue      HEC-RAS project (h5py optional, for HDF files)
    hmsvue      HEC-HMS project

Every tool module offers ``build_config(path, ...)``, ``show(cfg)`` and
``write_html(cfg, out)``. ``load`` / ``show`` / ``report`` below pick the tool
from the path.
"""

import importlib
from pathlib import Path

__version__ = '0.2.1'

__all__ = [
    'dfmvue',
    'dfmvue23',
    'epanetvue',
    'swmmvue',
    'flovue',
    'xbvue',
    'rasvue',
    'hmsvue',
    'describe',
    'cli',
    'detect',
    'load',
    'show',
    'report',
]

from . import dfmvue
from . import dfmvue23
from . import epanetvue
from . import swmmvue
from . import xbvue
from . import rasvue
from . import hmsvue
from . import describe
from . import cli
from .detection import ENGINES, detect

# flovue needs flopy, so it is imported on first use: `modelvue.flovue` or `from modelvue import flovue`


def __getattr__(name):
    if name == 'flovue':
        return importlib.import_module('.flovue', __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def get_tool(name):
    """Return the module of a tool (see ``ENGINES``)."""
    if name not in ENGINES:
        raise ValueError(f"unknown tool {name!r}; choose from {', '.join(ENGINES)}")
    return importlib.import_module(f".{name}", __name__)


def load(path, tool=None, **kwargs):
    """Read a model and return its configuration dict.

    ``tool`` is detected from ``path`` when not given. Extra keyword arguments go to
    the tool's ``build_config`` (e.g. ``only=`` for dfmvue/dfmvue23/flovue/rasvue/hmsvue,
    ``rpt=``/``out=`` for epanetvue/swmmvue, ``params_name=``/``only_case=`` for xbvue).
    """
    tool = tool or detect(path)
    cfg = get_tool(tool).build_config(str(path), **kwargs)
    cfg.setdefault("_tool", tool)
    return cfg


def _tool_of(cfg, tool):
    tool = tool or cfg.get("_tool")
    if not tool:
        raise ValueError("cannot tell which tool built this configuration; pass tool=")
    return get_tool(tool)


def show(cfg, tool=None, full=False):
    """Print a configuration (from :func:`load`) in the terminal."""
    _tool_of(cfg, tool).show(cfg, full)


def report(cfg, out, tool=None, open_browser=False):
    """Write a configuration (from :func:`load`) as a self-contained HTML page."""
    out = Path(out)
    _tool_of(cfg, tool).write_html(cfg, out, open_browser=open_browser)
    return out
