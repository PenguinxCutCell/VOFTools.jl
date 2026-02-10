# VOFTools.jl

A Julia translation of VOFTools (Version 5, January 2020) — analytical and
geometrical tools for 2D/3D Volume-of-Fluid (VOF) methods on general convex
and Cartesian grids.

This repository contains:
- `src/` — package source (core algorithms, mesh generators, utilities)
- `test/` — unit tests that match the original Fortran test suite

## Quick start

Activate the project and load the module (development/local use):

```bash
julia --project=. -e 'using Pkg; Pkg.activate("."); include("src/VOFTools.jl"); using .VOFTools'
```

Run the test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Usage example (from a Julia REPL or script):

```julia
using Pkg; Pkg.activate("."); include("src/VOFTools.jl"); using .VOFTools
poly = cubicmesh()
vol = toolv3d(poly)
println("cube volume = ", vol)
```

## Notes

- This package is a direct translation of the VOFTools Fortran implementation.
- See `src/meshes3d.jl` and `src/meshes2d.jl` for many provided test meshes.

## License and credits

All credits go to the original authors of VOFTools (J. Lopez, J. Hernandez, and collaborators).
Translated from VOFTools by J. Lopez and J. Hernandez. Original Fortran code
and authorship details remain with the original authors (see source headers).
