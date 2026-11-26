"""Small synthetic models written on the fly, so the tests need no external data."""
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"

SWMM_INP = """[TITLE]
tiny test network

[OPTIONS]
FLOW_UNITS           CMS
INFILTRATION         HORTON
FLOW_ROUTING         DYNWAVE
START_DATE           01/01/2024
START_TIME           00:00:00
END_DATE             01/01/2024
END_TIME             06:00:00
REPORT_STEP          00:05:00
WET_STEP             00:05:00
DRY_STEP             01:00:00
ROUTING_STEP         0:00:10

[RAINGAGES]
RG1   INTENSITY 0:15   1.0   TIMESERIES  TS1

[SUBCATCHMENTS]
S1    RG1   J1    2.0   50   200   0.5   0

[SUBAREAS]
S1    0.01  0.1   1.0   5.0   25   OUTLET

[INFILTRATION]
S1    75    5     4     7     0

[JUNCTIONS]
J1    10.0  3.0   0     0     0
J2     9.0  3.0   0     0     0

[OUTFALLS]
O1     8.0  FREE  NO

[CONDUITS]
C1    J1    J2    100   0.013  0  0  0  0
C2    J2    O1    100   0.013  0  0  0  0

[XSECTIONS]
C1    CIRCULAR  1.0  0  0  0  1
C2    CIRCULAR  1.0  0  0  0  1

[TIMESERIES]
TS1   0:00   0
TS1   1:00   20
TS1   2:00   5
TS1   3:00   0

[COORDINATES]
J1    0     0
J2    100   0
O1    200   0
"""

MDU = """[General]
Program      = D-Flow FM
fileVersion  = 1.09
fileType     = modelDef

[geometry]
NetFile      = tiny_net.nc
BedlevType   = 3
Kmx          = 0

[numerics]
CFLMax       = 0.7

[physics]
UnifFrictCoef = 0.023

[time]
RefDate      = 20240101
Tunit        = S
TStart       = 0
TStop        = 3600

[output]
OutputDir    = output
"""


def _xbeach(folder: Path) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    nx = 40
    xs = [i * 10.0 for i in range(nx + 1)]
    (folder / "x.grd").write_text(" ".join(f"{x:.1f}" for x in xs) + "\n")
    (folder / "bed.dep").write_text(" ".join(f"{-8 + 10 * (x / 400) ** 1.5:.3f}" for x in xs) + "\n")
    (folder / "jonswap.txt").write_text("Hm0 = 1.5\nTp = 8\nmainang = 270\ngammajsp = 3.3\ns = 10\n")
    (folder / "params.txt").write_text(
        f"nx = {nx}\nny = 0\nxfile = x.grd\ndepfile = bed.dep\nposdwn = -1\nvardx = 1\n"
        "wbctype = jonswap\nbcfile = jonswap.txt\ntstop = 3600\nmorfac = 5\n")
    return folder


@pytest.fixture
def epanet_inp(tmp_path):
    p = tmp_path / "Net1.inp"
    p.write_bytes((DATA / "Net1.inp").read_bytes())
    return p


@pytest.fixture
def swmm_inp(tmp_path):
    p = tmp_path / "tiny.inp"
    p.write_text(SWMM_INP)
    return p


@pytest.fixture
def xbeach_case(tmp_path):
    return _xbeach(tmp_path / "case1")


@pytest.fixture
def mdu_2d(tmp_path):
    p = tmp_path / "tiny.mdu"
    p.write_text(MDU)
    return p


@pytest.fixture
def mdu_1d(tmp_path):
    p = tmp_path / "tiny1d.mdu"
    p.write_text(MDU.replace("Kmx          = 0", "Kmx          = 0\nCrossDefFile = crsdef.ini\nCrossLocFile = crsloc.ini"))
    return p


# --------------------------------------------------------------------------- HEC-RAS
def _fw(vals, width, per_line):
    out = []
    for k in range(0, len(vals), per_line):
        out.append("".join(f"{v:>{width}}" for v in vals[k:k + per_line]))
    return "\n".join(out)


def _ras_xs(rs, x0):
    se = [0, 105, 20, 100, 40, 95, 60, 100, 80, 106]
    cut = [x0, 0.0, x0, 80.0]
    return (f"Type RM Length L Ch R = 1 ,{rs},100,100,100\n"
            f"XS GIS Cut Line=2\n{_fw([f'{v:.1f}' for v in cut], 16, 4)}\n"
            f"#Sta/Elev= 5 \n{_fw(se, 8, 10)}\n"
            "#Mann= 3 , 0 , 0 \n" + _fw([0, .06, 0, 20, .035, 0, 60, .06, 0], 8, 9) + "\n"
            "Bank Sta=20,60\nExp/Cntr=0.3,0.1\n\n")


@pytest.fixture
def ras_project(tmp_path):
    d = tmp_path / "ras"
    d.mkdir()
    (d / "Tiny.prj").write_text("Proj Title=Tiny test\nCurrent Plan=p01\nSI Units\nGeom File=g01\nUnsteady File=u01\n"
                                "Plan File=p01\nBEGIN DESCRIPTION:\nA test project.\nEND DESCRIPTION:\n")
    reach_xy = _fw([f"{v:.1f}" for v in (10.0, 40.0, 200.0, 40.0)], 16, 4)
    (d / "Tiny.g01").write_text("Geom Title=Tiny geometry\nProgram Version=6.30\n\nRiver Reach=Creek           ,Main            \n"
                                f"Reach XY= 2 \n{reach_xy}\n\n" + _ras_xs("200", 10.0) + _ras_xs("100", 110.0) + _ras_xs("0", 200.0))
    flows = [10, 20, 40, 30, 20, 10]
    (d / "Tiny.u01").write_text("Flow Title=Tiny flow\nProgram Version=6.30\n"
                                "Boundary Location=Creek           ,Main            ,200     ,        ,                ,                \n"
                                f"Interval=1HOUR\nFlow Hydrograph= {len(flows)} \n{_fw(flows, 8, 10)}\nUse DSS=False\n"
                                "Boundary Location=Creek           ,Main            ,0       ,        ,                ,                \n"
                                "Friction Slope=0.001,0\nInitial Flow Loc=Creek           ,Main            ,200     ,10\n")
    (d / "Tiny.p01").write_text("Plan Title=Tiny plan\nProgram Version=6.30\nShort Identifier=tiny\n"
                                "Simulation Date=01JAN2024,0000,01JAN2024,0500\nGeom File=g01\nFlow File=u01\nSubcritical Flow\n"
                                "Computation Interval=1MIN\nOutput Interval=10MIN\nRun HTab= 1 \nRun UNet= 1 \n")
    (d / "Tiny.p01.comp_msgs.txt").write_text("Simulation started at: 01Jan2024 10:00\nFinished Unsteady Flow Simulation\n"
                                              "Complete Process  12.3 sec\n")
    (d / "gis.prj").write_text('PROJCS["not a ras project"]')
    return d / "Tiny.prj"


# --------------------------------------------------------------------------- HEC-HMS
@pytest.fixture
def hms_project(tmp_path):
    d = tmp_path / "hms"
    (d / "results").mkdir(parents=True)
    (d / "tiny.hms").write_text("Project: tiny\n     Description: Tiny HMS\n     Version: 4.11\nEnd:\n\n"
                                "Precipitation: Met 1\n     Filename: Met_1.met\nEnd:\n\n"
                                "Basin: Basin 1\n     Filename: Basin_1.basin\nEnd:\n\n"
                                "Control: Ctl\n     FileName: Ctl.control\nEnd:\n")
    (d / "Basin_1.basin").write_text(
        "Basin: Basin 1\n     Unit System: Metric\nEnd:\n\n"
        "Subbasin: S1\n     Canvas X: 0\n     Canvas Y: 100\n     Area: 2.5\n     Downstream: J1\n"
        "     LossRate: SCS\n     Transform: SCS\n     Baseflow: None\nEnd:\n\n"
        "Subbasin: S2\n     Canvas X: 100\n     Canvas Y: 100\n     Area: 1.5\n     Downstream: J1\n"
        "     LossRate: SCS\n     Transform: Clark\nEnd:\n\n"
        "Junction: J1\n     Canvas X: 50\n     Canvas Y: 0\nEnd:\n\n"
        "Reach: R-lost\n     Canvas X: 60\n     Canvas Y: 10\n     Downstream: Nowhere\n     Route: Lag\nEnd:\n")
    (d / "Met_1.met").write_text("Meteorology: Met 1\n     Precipitation Method: Specified Average\n     Use Basin Model: Basin 1\nEnd:\n\n"
                                 "Subbasin: S1\n     Gage: G1\nEnd:\n")
    (d / "Ctl.control").write_text("Control: Ctl\n     Start Date: 1 January 2024\n     Start Time: 00:00\n"
                                   "     End Date: 2 January 2024\n     End Time: 00:00\n     Time Interval: 15\nEnd:\n")
    (d / "tiny.run").write_text("Run: Run 1\n     Log File: Run_1.log\n     Basin: Basin 1\n     Precip: Met 1\n     Control: Ctl\nEnd:\n")
    (d / "Run_1.log").write_text("NOTE 10000: Began computing.\nWARNING 20000: Something odd.\nNOTE 10001: Finished computing simulation run.\n")
    (d / "results" / "RUN_Run_1.results").write_text(
        '<?xml version="1.0" encoding="UTF-8"?><RunResults><RunName>Run 1</RunName>'
        '<BasinElement name="J1" type="Junction"><Statistics>'
        '<StatisticMeasure type="Outflow Maximum" value="12.5" units="M3/S" />'
        '<StatisticMeasure type="Outflow Maximum Time" value="01Jan2024, 06:00" />'
        '</Statistics></BasinElement></RunResults>')
    return d / "tiny.hms"
