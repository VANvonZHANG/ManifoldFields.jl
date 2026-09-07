#!/usr/bin/env python
"""Authoring script: convert the local MPAS oQU480 grid+data pair into a
single UGRID-conformant NetCDF file (examples/data/oQU480.ugrid.nc).

Maintainer tool: run inside the `geo` conda environment (uxarray 2026.7.0).
Tutorial users never need this script; they fetch the packed artifact.

    /home/zhangfan/anaconda3/envs/geo/bin/python examples/scripts/prepare_ugrid.py \
        [--grid PATH] [--data PATH] [--out PATH]

Defaults point at the UXarray-Tutorial copies already on disk in this workspace.
"""
import argparse
from pathlib import Path

import uxarray as ux

HERE = Path(__file__).resolve().parent
DEFAULTS = dict(
    grid="/home/zhangfan/Project/20260328_HEMCOManifold/UXarray-Tutorial/notebooks/oQU480.grid.nc",
    data="/home/zhangfan/Project/20260328_HEMCOManifold/UXarray-Tutorial/notebooks/oQU480.data.nc",
    out=HERE.parent / "data" / "oQU480.ugrid.nc",
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid", default=DEFAULTS["grid"])
    ap.add_argument("--data", default=DEFAULTS["data"])
    ap.add_argument("--out", type=Path, default=DEFAULTS["out"])
    args = ap.parse_args()

    uxds = ux.open_dataset(args.grid, args.data)
    out = uxds.uxgrid.to_xarray()  # UGRID conventions (default grid_format="ugrid")
    for name in uxds.data_vars:
        out[name] = uxds[name].variable  # data vars share the face dimension
        out[name].attrs.update(uxds[name].attrs)  # keep units/long_name etc.
        out[name].attrs["mesh"] = "grid_topology"  # UGRID data-var -> topology pointer
    # The MPAS source carries no units attr for bottomDepth; state the known unit
    # (ocean bottom depth, meters) so the written file is self-describing.
    if "bottomDepth" in out.variables:
        out["bottomDepth"].attrs.setdefault("units", "m")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_netcdf(args.out)
    print(f"wrote {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")
    print("variables:", sorted(out.variables))


if __name__ == "__main__":
    main()
