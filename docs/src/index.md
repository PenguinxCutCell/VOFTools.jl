```@meta
CurrentModule = VOFTools
```

# VOFTools

## Installation

```julia
using Pkg
Pkg.add("VOFTools")
```

For local development:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Overview

`VOFTools.jl` provides analytical and geometrical tools for 2D/3D VOF workflows
on convex and Cartesian-style grids.

The package includes:

- Polygon/polyhedron geometry types (`Polygon2D`, `Polyhedron3D`)
- Geometric operators (volume/area computation, truncation by planes/lines)
- Volume conservation enforcement routines
- A large set of predefined 2D/3D test meshes

## Quick Example

```julia
using VOFTools

poly3d = cubicmesh()
v0 = toolv3d(poly3d)

# Truncate with x + y + z - 1 = 0
n = (1.0, 1.0, 1.0)
d = sqrt(sum(x -> x^2, n))
icontn, icontp = inte3d!(poly3d, -1.0 / d, n[1] / d, n[2] / d, n[3] / d)
v1 = toolv3d(poly3d)

println((icontn, icontp, v0, v1))
```

## Main Sections

- [Geometry types](@ref geometry-types)
- [2D tools](@ref tools2d)
- [3D tools](@ref tools3d)
- [Mesh generators](@ref mesh-generators)
- [Reference](@ref reference)

## Performance

By using task-local workspaces and preallocating necessary arrays, this implementation minimizes memory allocations and optimizes performance.
Allocations free and fast execution times are expected for the provided benchmark cases (see `benchmark/profile_case.jl`).

## Reference

Joaquín López, Julio Hernández, Pablo Gómez, Claudio Zanzi, Rosendo Zamora, VOFTools 5: An extension to non-convex geometries of calculation tools for volume of fluid methods, Computer Physics Communications, Volume 252, 2020, 107277, ISSN 0010-4655, https://doi.org/10.1016/j.cpc.2020.107277.