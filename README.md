# MODELVUE

A collection of Python tools that read hydraulic and hydrological model setups and show
their configuration, in the terminal or as a self-contained HTML report.

## Installation

```bash
pip install -e .
pip install -e ".[modflow]"    # + flopy, numpy, matplotlib for flovue
pip install -e ".[all]"        # + netCDF4 for NetCDF-4 map / his / xboutput files
```

`dfmvue`, `dfmvue23`, `epanetvue`, `swmmvue` and `xbvue` need only the Python standard
library (classic NetCDF-3 files are read by a built-in reader); numpy makes them faster.

## Features

- **dfmvue**: Configuration of a Delft3D FM / D-HYDRO 1D2D model (`.mdu`, `.dsproj`)
- **dfmvue23**: Configuration of a Delft3D FM 2D/3D model, incl. sediment/morphology and D-Waves
- **epanetvue**: Configuration of an EPANET 2.x model, with network checks and results (`.inp`, `.rpt`, `.out`)
- **swmmvue**: Configuration of an EPA SWMM 5 model, with network checks and results (`.inp`, `.rpt`, `.out`)
- **flovue**: Configuration of a MODFLOW 6 / 2005 / NWT / USG, MT3D or SEAWAT model, read with FloPy
- **xbvue**: Configuration of an XBeach model or OET `.config` case set
- **rasvue**: Configuration of a HEC-RAS project: geometry (1D and 2D), steady / unsteady flow, plans, run messages
- **hmsvue**: Configuration of a HEC-HMS project: basin and meteorologic models, controls, runs and their results
- **modelvue**: Detect the model type of a file or folder and run the matching tool

Every tool prints to the terminal by default. `--html [FILE]` writes an HTML report instead
and opens it in the browser (`--no-browser` to skip that); `--json FILE` also dumps the
configuration to JSON.

## Usage

### Check version

```bash
modelvue --version
```

### List all available commands

```bash
modelvue-info
modelvue-info swmmvue
```

### modelvue

Detect the model type and run the matching tool. `.inp` files are told apart by their
sections (SWMM or EPANET), binary `.out` files by their header, and an `.mdu` goes to
`dfmvue` when it names 1D input files or its net file holds a 1D network, otherwise to
`dfmvue23`. A folder holding several kinds of model is reported as ambiguous.

```bash
modelvue model.inp
modelvue FlowFM.mdu --html
modelvue detect folder
modelvue dfmvue23 FlowFM.mdu --all     # same as: dfmvue23 FlowFM.mdu --all
```

### dfmvue

Delft3D FM 1D2D model: MDU, UGRID network/mesh, branches, cross sections, roughness,
storage nodes, structures, external forcing (`.ext` / `.bc`), polylines, initial fields,
and the `.dia`, map and his output when the model has been run.

```bash
dfmvue FlowFM.mdu
dfmvue project_folder                  # finds the .mdu (incl. *.dsproj_data)
dfmvue Project.dsproj --model FlowFM   # pick one MDU when several exist
dfmvue FlowFM.mdu --all                # every MDU keyword
dfmvue FlowFM.mdu --html report.html --no-browser
```

### dfmvue23

Delft3D FM 2D/3D model: mesh, boundary conditions (time series, harmonic, astronomic,
3D profiles), thin dams, fixed weirs, observation points, cross sections, structures,
vertical layers, physics and numerics, sediment / morphology (`.sed` / `.mor`), and a
coupled D-Waves model (`.mdw`). Option codes are decoded from the MDU's own comments.

```bash
dfmvue23 FlowFM.mdu
dfmvue23 FlowFM.mdu --html
dfmvue23 FlowFM.mdu --json config.json
```

### epanetvue

EPANET 2.x model: every section of the `.inp`, network checks (references, connectivity,
supply sources, valves, curves, patterns), and the `.rpt` and binary `.out` results.

```bash
epanetvue Net1.inp
epanetvue folder                       # every .inp in the folder (and below)
epanetvue Net1.inp --rpt run.rpt --out run.out
epanetvue Net1.out                     # only a binary output file
epanetvue Net1.inp --all --html
```

### swmmvue

EPA SWMM 5 model: every section of the `.inp`, network checks, and the `.rpt`
(continuity errors, flooding, instability) and binary `.out` results.

```bash
swmmvue model.inp
swmmvue model.inp --rpt run.rpt --out run.out
swmmvue model.inp --html
```

### flovue

MODFLOW model read with FloPy: MODFLOW 6 simulations and classic name-file models
(MODFLOW-2005, NWT, USG, MODFLOW-2000, MT3DMS / MT3D-USGS, SEAWAT). Grid, time
discretisation, solver, every package, FloPy's model check, and the listing, head and
budget files. Needs `flopy` and `numpy` (`matplotlib` for the figures).

```bash
flovue mfsim.nam
flovue model.nam --version mfnwt
flovue folder --model sagehen
flovue mfsim.nam --no-check --html
```

### xbvue

XBeach model: `params.txt`, grid and bathymetry, wave boundary, water-level and
vegetation files, OET test-bed `.config` case sets, and `xboutput.nc` when present.

```bash
xbvue params.txt
xbvue case_folder
xbvue Holland1976 --case Raai1085      # one case of a set
xbvue Holland1976 --details            # full detail for every case
```

### rasvue

HEC-RAS project: the `.prj` and every file it lists.

- **Geometry** (`.gNN`): river reaches, cross sections, bridges, culverts, inline and
  lateral structures, junctions, storage areas, 2D flow areas, SA/2D connections, BC
  lines, break lines and pump stations. The report adds a map, longitudinal profiles,
  sample cross sections and data checks (bank stations, Manning's n, station order).
- **Steady flow** (`.fNN`) and **unsteady flow** (`.uNN`): boundary and initial
  conditions, plus hydrographs and rating curves stored as tables.
- **Plans** (`.pNN`): simulation window, intervals, solver settings, and the computation
  messages and volume accounting when the plan has been run.
- **RAS Mapper** (`.rasmap`): layers and projection.

With `h5py` installed the geometry HDF (`.gNN.hdf`) adds 2D mesh statistics, and the
plan HDF (`.pNN.hdf`) adds a results summary. A file name or folder is enough: GIS
`.prj` projection files are skipped.

```bash
rasvue Muncie.prj
rasvue Muncie.g04                      # any file of the project
rasvue "7.0.2 Examples"                # every project below
rasvue "7.0.2 Examples" --project Muncie
rasvue Muncie.prj --all --html
```

### hmsvue

HEC-HMS project: the `.hms` and its basin models (element network, loss / transform /
baseflow / routing methods, map with the GIS layers), meteorologic models (methods,
gage weights), control specifications, simulation runs, time-series gages, paired data and
grids. For computed runs it reads `results/RUN_<name>.results` (peak flow and time,
volumes, Nash-Sutcliffe and bias against observed flow) and the run log. DSS files are
listed, not read.

```bash
hmsvue castro.hms
hmsvue samples                         # every .hms below
hmsvue castro.hms --all                # every element parameter
hmsvue castro.hms --html
```

## Python API

```python
import modelvue

cfg = modelvue.load("Net1.inp")        # tool picked from the path
modelvue.show(cfg)
modelvue.report(cfg, "net1.html")

from modelvue import swmmvue           # or use one tool directly
cfg = swmmvue.build_config("model.inp", rpt="run.rpt", out="run.out")
swmmvue.write_html(cfg, "model.html", open_browser=False)
```

## Development

```bash
pip install -e ".[dev]"
pytest
ruff check modelvue
```
