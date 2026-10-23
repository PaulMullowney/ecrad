#!/usr/bin/env python3
"""Merge the per-part outputs of ecrad_ifs_blocked run with ECRAD_NPARTS=N.

Each part writes a full-size output file but only computes the columns of its
own blocks. Block b (0-based) belongs to part (b // chunk) % nparts, the same
rule as ecrad_ifs_driver_blocked.F90, and every variable with a "column"
dimension is taken from the owning part.

usage: merge_parts.py --nparts N --nproma P [--chunk C] OUTPUT PART0 ... PARTN-1
"""
import argparse
import shutil

import netCDF4
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--nparts', type=int, required=True)
    parser.add_argument('--nproma', type=int, required=True)
    parser.add_argument('--chunk', type=int, default=8)
    parser.add_argument('output')
    parser.add_argument('parts', nargs='+')
    args = parser.parse_args()
    if len(args.parts) != args.nparts:
        parser.error(f'expected {args.nparts} part files, got {len(args.parts)}')

    shutil.copyfile(args.parts[0], args.output)
    out = netCDF4.Dataset(args.output, 'r+')
    ncol = out.dimensions['column'].size
    owner = (np.arange(ncol) // args.nproma // args.chunk) % args.nparts

    parts = [netCDF4.Dataset(p) for p in args.parts]
    for path, part in zip(args.parts, parts):
        if part.dimensions['column'].size != ncol:
            raise SystemExit(f'{path}: column dimension differs from {args.parts[0]}')

    for name, var in out.variables.items():
        if 'column' not in var.dimensions:
            continue
        axis = var.dimensions.index('column')
        merged = var[:]
        for ipart in range(1, args.nparts):
            mask = owner == ipart
            if not mask.any():
                continue
            index = [slice(None)] * var.ndim
            index[axis] = mask
            merged[tuple(index)] = parts[ipart].variables[name][:][tuple(index)]
        var[:] = merged

    for part in parts:
        part.close()
    out.close()


if __name__ == '__main__':
    main()
