# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

---

## [0.2.1] - 2026-11-26

### Added

- All HTML reports have a light / dark theme switch in the header. Reports open in the
  light theme; a switch to dark is remembered by the browser for every modelvue report

---

## [0.2.0] - 2026-11-26

### Added

- `rasvue`: configuration of a HEC-RAS project - geometry (reaches, cross sections,
  bridges, culverts, inline / lateral structures, junctions, storage areas, 2D flow
  areas, SA/2D connections, BC lines, break lines, pump stations), steady and unsteady
  flow, plans, RAS Mapper, computation messages; map, longitudinal profiles and cross
  section plots in the HTML report; checks of bank stations, Manning's n and station
  order. Standard library only; with h5py also the geometry HDF (2D mesh statistics)
  and the plan HDF (results summary)
- `hmsvue`: configuration of a HEC-HMS project - basin models (element network and
  methods, schematic with the basin's GIS layers), meteorologic models, control
  specifications, runs with their results (peak flow, volumes, observed-flow
  statistics) and logs, gages, paired data and grids; network and reference checks
- `modelvue` detects HEC-RAS (`.prj` starting with `Proj Title=`, `.rasmap`) and
  HEC-HMS (`.hms`, `.basin`, `.met`, `.control`) projects; folder search now goes
  5 levels deep
- `hdf` extra (`pip install modelvue[hdf]`) for h5py

---

## [0.1.0] - 2026-11-26

### Added

- Initial release of `modelvue`, built from the prototype scripts
- `dfmvue`: configuration of a Delft3D FM / D-HYDRO 1D2D model (MDU, UGRID network,
  cross sections, roughness, structures, forcing, .dia / map / his output)
- `dfmvue23`: configuration of a Delft3D FM 2D/3D model, incl. sediment /
  morphology and a coupled D-Waves (SWAN) model
- `epanetvue`: configuration of an EPANET 2.x model, with network checks and the
  .rpt / binary .out results
- `swmmvue`: configuration of an EPA SWMM 5 model, with network checks and the
  .rpt / binary .out results
- `flovue`: configuration of a MODFLOW 6 / 2005 / NWT / USG, MT3D or SEAWAT model,
  read with FloPy
- `xbvue`: configuration of an XBeach model or OET `.config` case set
- Every tool prints to the terminal, or writes a self-contained HTML report
  (`--html`) and a JSON dump (`--json`)
- `modelvue` command: detects the model type from a file or folder and runs the
  matching tool; `modelvue detect PATH` only reports the type
- `modelvue-info` command: lists the tools and their usage
- Python API: `modelvue.load()`, `modelvue.show()`, `modelvue.report()`
- GitHub Actions workflow for publishing to PyPI and TestPyPI

### Fixed

- `dfmvue23`: no longer crashes on an MDU without `Idensform`, `ICdtyp` or
  `TimeStepType` (these were unhandled in the prototype script)
