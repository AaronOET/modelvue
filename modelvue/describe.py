#!/usr/bin/env python
"""
Command-line utility to display descriptions of modelvue functionality.
"""

import argparse
from textwrap import dedent

TOOL_DESCRIPTIONS = {
    'dfmvue': """
        Show the configuration of a Delft3D FM (D-Flow FM / D-HYDRO) 1D2D model.

        Reads the MDU and the files it points to: UGRID network/mesh (*_net.nc),
        branches, cross sections, roughness, storage nodes, structures, external
        forcing (.ext / .bc), polylines and initial fields. If the model has been
        run, also the .dia, map and his output. Standard library only.

        Examples:
            dfmvue FlowFM.mdu
            dfmvue project_folder              # finds the .mdu (incl. *.dsproj_data)
            dfmvue Project.dsproj --model FlowFM
            dfmvue FlowFM.mdu --all            # every MDU keyword
            dfmvue FlowFM.mdu --html           # HTML report, opens in the browser
            dfmvue FlowFM.mdu --json config.json
    """,
    'dfmvue23': """
        Show the configuration of a Delft3D FM 2D or 3D model.

        Reads the MDU, UGRID mesh, boundary conditions (time series, harmonic,
        astronomic, 3D profiles), thin dams, fixed weirs, observation points and
        cross sections, structures, initial fields, sediment / morphology, and a
        coupled D-Waves (SWAN) model. Option codes are decoded from the MDU's own
        comments. Summarises the .dia, map and his files when present.

        Examples:
            dfmvue23 FlowFM.mdu
            dfmvue23 project_folder --all
            dfmvue23 FlowFM.mdu --html report.html --no-browser
    """,
    'epanetvue': """
        Show the configuration of an EPANET 2.x water distribution model.

        Parses the .inp file, checks the network (references, connectivity,
        supply sources, valves, curves, patterns) and, if the model has been run,
        reads the .rpt and the binary .out file (pressures, velocities, tank
        levels, pump energy, water quality). Standard library only.

        Examples:
            epanetvue Net1.inp
            epanetvue folder                   # every .inp in the folder
            epanetvue Net1.inp --rpt run.rpt --out run.out
            epanetvue Net1.out                 # only a binary output file
            epanetvue Net1.inp --html
    """,
    'swmmvue': """
        Show the configuration of an EPA SWMM 5 model.

        Parses the .inp file section by section, checks the network and, if the
        model has been run, reads the .rpt (continuity errors, flooding,
        instability) and the binary .out file (peak values, time series,
        profiles). Standard library only.

        Examples:
            swmmvue model.inp
            swmmvue model.inp --all            # every row of every section
            swmmvue model.inp --rpt run.rpt --out run.out
            swmmvue model.inp --html
    """,
    'flovue': """
        Show the configuration of a MODFLOW model, read with FloPy.

        MODFLOW 6 simulations (every model in mfsim.nam) and classic name-file
        models: MODFLOW-2005, NWT, USG, MODFLOW-2000, MT3DMS / MT3D-USGS and
        SEAWAT. Shows grid, time discretisation, solver, every package, FloPy's
        model check, and the listing / head / budget files when present.
        Needs flopy and numpy (matplotlib for figures).

        Examples:
            flovue mfsim.nam
            flovue model.nam --version mfnwt
            flovue folder --model sagehen
            flovue mfsim.nam --no-check --html
    """,
    'xbvue': """
        Show the configuration of an XBeach model or case set.

        Reads params.txt, grid and bathymetry files, wave boundary files,
        water-level and vegetation files, and the Deltares OET test-bed .config
        file that lists several runs. Summarises xboutput.nc when present.

        Examples:
            xbvue params.txt
            xbvue case_folder
            xbvue Holland1976 --case Raai1085  # one case of a set
            xbvue Holland1976 --details        # full detail for every case
            xbvue case_folder --html
    """,
    'rasvue': """
        Show the configuration of a HEC-RAS project.

        Reads the project (.prj) and every file it lists: geometries (reaches,
        cross sections, bridges, culverts, inline / lateral structures,
        junctions, storage areas, 2D flow areas, connections, BC lines),
        steady flow, unsteady flow (boundary and initial conditions,
        hydrographs), plans (simulation window, intervals, solver settings)
        and RAS Mapper. Reports computation messages and volume accounting.
        Text files need only the standard library; with h5py also the 2D mesh
        (.gNN.hdf) and a results summary (.pNN.hdf).

        Examples:
            rasvue Muncie.prj
            rasvue project_folder              # every HEC-RAS project below
            rasvue examples --project Muncie   # one project of many
            rasvue Muncie.prj --all            # every plan, XS and boundary
            rasvue Muncie.prj --html
    """,
    'hmsvue': """
        Show the configuration of a HEC-HMS project.

        Reads the project (.hms) and its basin models (element network, loss /
        transform / baseflow / routing methods), meteorologic models, control
        specifications, simulation runs, gages, paired data and grids, and -
        for computed runs - the results (peak flow, volume, observed-flow
        statistics such as Nash-Sutcliffe) and the run logs.
        Standard library only (DSS files are listed, not read).

        Examples:
            hmsvue castro.hms
            hmsvue samples_folder              # every .hms below
            hmsvue castro.hms --all            # every element parameter
            hmsvue castro.hms --html
    """,
}

GENERAL = """
    modelvue picks the right tool for you:

        modelvue PATH [options]            detect the model type and run its tool
        modelvue detect PATH               only print which tool would read PATH
        modelvue TOOL PATH [options]       run one tool (same as calling TOOL directly)

    Common options: --html [FILE], --no-browser, --json FILE.
"""


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='modelvue-info',
        description='Display descriptions of modelvue commands',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        'tool',
        nargs='?',
        choices=list(TOOL_DESCRIPTIONS.keys()),
        help='Tool name to describe (omit to list all tools)',
    )
    args = parser.parse_args(argv)

    if args.tool:
        print(f"\n--- {args.tool} ---")
        print(dedent(TOOL_DESCRIPTIONS[args.tool]))
    else:
        print("\nMODELVUE — model configuration viewers\n")
        for tool, desc in TOOL_DESCRIPTIONS.items():
            first_line = dedent(desc).strip().splitlines()[0]
            print(f"  {tool:<14} {first_line}")
        print(dedent(GENERAL))
        print("Run 'modelvue-info <tool>' for details on a specific command.")


if __name__ == "__main__":
    main()
