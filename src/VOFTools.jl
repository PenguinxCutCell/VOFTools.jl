"""
    VOFTools

A Julia package of analytical and geometrical tools for 2D/3D Volume of Fluid
(VOF) methods in general convex grids and Cartesian geometry.

Translated from VOFTools 5 (January 2020) by J. Lopez and J. Hernandez.

# References
- [1] J. Lopez and J. Hernandez, Analytical and geometrical tools for 3D volume
      of fluid methods in general grids, J. Comput. Phys. 227 (2008) 5939-5948.
- [2] J. Lopez, P. Gomez, J. Hernandez, F. Faura, A two-grid adaptive volume of
      fluid approach for dendritic solidification, Comput. Fluids 86 (2013) 326-342.
- [3] J. Lopez and J. Hernandez, P. Gomez, F. Faura, A new volume conservation
      enforcement method for PLIC reconstruction in general convex grids, J.
      Comput. Phys. 316 (2016) 338-359.
"""
module VOFTools

using LinearAlgebra

export Polygon2D, Polyhedron3D
export toolv2d, toolv3d
export inte2d!, inte3d!
export enforv2d, enforv2d!, enforv3d, enforv3d!
export enforv2dsz, enforv3dsz
export cppol2d, cppol2d!, cppol3d, cppol3d!
export restore2d!, restore3d!
export newpol2d!, newpol3d!
export eqsol3d, newton3d
export dist2d, dist3d
export initf2d, initf3d

# 2D mesh generators
export squaremesh, hexagomesh, trianglemesh, quadranglemesh, pentagonmesh, hexagonmesh
export ncquadranglemesh, ncpentagonmesh, nchexagonmesh, ncshexagonmesh
export nchollowedsquare, ncmultisquare

# 3D mesh generators
export cubicmesh, hexahemesh, tetramesh, dodecamesh, icosamesh, complexmesh
export distortcubicmesh, ncpentapyramid, nccubicpyramid
export ncscubicmesh, nchexahemesh, ncdodecamesh, ncicosamesh
export nchollowedcube, drilledcube, zigzagmesh, voftoolslogo

include("types.jl")
include("utils.jl")
include("tools3d.jl")
include("tools2d.jl")
include("meshes3d.jl")
include("meshes2d.jl")

end # module VOFTools
