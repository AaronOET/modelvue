#!/usr/bin/env python
"""
Top-level `modelvue` command. Running `modelvue -h` (or with no arguments)
prints the same information shown by `modelvue-info`.

    modelvue PATH [options]        detect the model type, then run its tool
    modelvue TOOL PATH [options]   run one tool explicitly (dfmvue, epanetvue, ...)
    modelvue TOOL                  describe one tool
    modelvue detect PATH           only print the detected tool
"""

import importlib
import importlib.metadata
import sys

from . import describe
from .detection import candidates, detect


def _version():
    try:
        return importlib.metadata.version('modelvue')
    except importlib.metadata.PackageNotFoundError:
        from . import __version__
        return __version__


def _run(tool, args):
    try:
        mod = importlib.import_module(f".{tool}", __package__)
    except ImportError as e:
        print(f"modelvue: {tool} needs a package that is not installed: {e}", file=sys.stderr)
        return 1
    mod.main(args, prog=tool)
    return 0


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)

    if not args or args[0] in ('-h', '--help'):
        describe.main([])
        return 0

    if args[0] in ('-v', '--version'):
        print(f"modelvue {_version()}")
        return 0

    if args[0] == 'detect':
        if len(args) < 2:
            print("usage: modelvue detect PATH", file=sys.stderr)
            return 2
        found = candidates(args[1])
        if not found:
            print(f"no supported model found at {args[1]}", file=sys.stderr)
            return 1
        for tool, files in found.items():
            print(f"{tool}\t{len(files)} file(s), e.g. {files[0]}")
        return 0

    if args[0] in describe.TOOL_DESCRIPTIONS:
        if len(args) == 1:
            describe.main([args[0]])
            return 0
        return _run(args[0], args[1:])

    try:
        tool = detect(args[0])
    except ValueError as e:
        print(f"modelvue: {e}", file=sys.stderr)
        print("Run 'modelvue -h' to see available tools.", file=sys.stderr)
        return 2
    print(f"[modelvue: reading with {tool}]", file=sys.stderr)
    return _run(tool, args)


if __name__ == "__main__":
    sys.exit(main())
