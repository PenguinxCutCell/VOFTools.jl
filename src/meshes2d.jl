# ---------------------------------------------------------------------------
#  2D mesh generators for VOFTools.jl
#  Translated from mesh.f (Version 5, January 2020)
#  Copyright (C) 2020 J. Lopez and J. Hernandez
# ---------------------------------------------------------------------------

"""
    squaremesh() -> Polygon2D

Create a unit square `[0,1]^2` polygon.
"""
function squaremesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 4; ntp = 4
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:] = [0.0, 0.0]
    vertp[2,:] = [1.0, 0.0]
    vertp[3,:] = [1.0, 1.0]
    vertp[4,:] = [0.0, 1.0]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    hexagomesh() -> Polygon2D

Create a regular hexagonal polygon.
"""
function hexagomesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 6; ntp = 6
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:] = [0.5,       0.0]
    vertp[2,:] = [0.9330127, 0.25]
    vertp[3,:] = [0.9330127, 0.75]
    vertp[4,:] = [0.5,       1.0]
    vertp[5,:] = [0.066987298, 0.75]
    vertp[6,:] = [0.066987298, 0.25]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    trianglemesh() -> Polygon2D

Create a triangular cell.
"""
function trianglemesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 3; ntp = 3
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:] = [0.0,  0.0]
    vertp[2,:] = [0.72, 0.13]
    vertp[3,:] = [1.0,  1.0]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    quadranglemesh() -> Polygon2D

Create a quadrangular cell.
"""
function quadranglemesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 4; ntp = 4
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:] = [0.0,  0.0]
    vertp[2,:] = [1.0,  0.13]
    vertp[3,:] = [0.72, 1.0]
    vertp[4,:] = [0.13, 0.56]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    pentagonmesh() -> Polygon2D

Create a pentagonal cell.
"""
function pentagonmesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 5; ntp = 5
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:] = [0.0,  0.15]
    vertp[2,:] = [0.09, 0.0]
    vertp[3,:] = [1.0,  0.13]
    vertp[4,:] = [0.16, 1.0]
    vertp[5,:] = [0.04, 0.77]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    hexagonmesh() -> Polygon2D

Create an irregular hexagonal cell.
"""
function hexagonmesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 6; ntp = 6
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:] = [0.0,  0.0]
    vertp[2,:] = [0.66, 0.03]
    vertp[3,:] = [1.0,  0.22]
    vertp[4,:] = [0.9,  0.77]
    vertp[5,:] = [0.72, 1.0]
    vertp[6,:] = [0.33, 0.86]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

# ========================== New 2D Mesh Generators ===========================

"""
    ncquadranglemesh() -> Polygon2D

Create a non-convex quadrangular cell (4 vertices).
"""
function ncquadranglemesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 4; ntp = 4
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:] = [0.0,  0.0]
    vertp[2,:] = [1.0,  0.13]
    vertp[3,:] = [0.72, 1.0]
    vertp[4,:] = [0.75, 0.56]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    ncpentagonmesh() -> Polygon2D

Create a non-convex pentagonal cell (5 vertices).
"""
function ncpentagonmesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 5; ntp = 5
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:] = [0.0,  0.0]
    vertp[2,:] = [0.49, 0.22]
    vertp[3,:] = [1.0,  0.13]
    vertp[4,:] = [0.16, 1.0]
    vertp[5,:] = [0.04, 0.77]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    nchexagonmesh() -> Polygon2D

Create a non-convex hexagonal cell (6 vertices).
"""
function nchexagonmesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 6; ntp = 6
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:] = [0.0,  0.0]
    vertp[2,:] = [0.33, 0.43]
    vertp[3,:] = [1.0,  0.22]
    vertp[4,:] = [0.9,  0.77]
    vertp[5,:] = [0.72, 1.0]
    vertp[6,:] = [0.33, 0.86]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    ncshexagonmesh() -> Polygon2D

Create a non-convex stellated hexagonal cell (12 vertices).
"""
function ncshexagonmesh()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 12; ntp = 12
    for i in 1:ntv; ipv[i] = i; end
    vertp[1,:]  = [0.0,  0.2]
    vertp[2,:]  = [0.25, 0.25]
    vertp[3,:]  = [0.4,  0.0]
    vertp[4,:]  = [0.55, 0.3]
    vertp[5,:]  = [0.8,  0.35]
    vertp[6,:]  = [0.65, 0.55]
    vertp[7,:]  = [0.8,  0.75]
    vertp[8,:]  = [0.5,  0.7]
    vertp[9,:]  = [0.35, 0.9]
    vertp[10,:] = [0.25, 0.7]
    vertp[11,:] = [-0.1, 0.55]
    vertp[12,:] = [0.2,  0.45]
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    nchollowedsquare() -> Polygon2D

Create a non-convex hollowed square (10 vertices): outer square contour
followed by inner square hole with reversed winding.
Vertex 5 repeats vertex 1 to close the outer contour;
vertex 10 repeats vertex 6 to close the inner contour.
"""
function nchollowedsquare()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 10; ntp = 10
    for i in 1:ntv; ipv[i] = i; end
    # Outer square
    vertp[1,:] = [0.0,  0.0]
    vertp[2,:] = [1.0,  0.0]
    vertp[3,:] = [1.0,  1.0]
    vertp[4,:] = [0.0,  1.0]
    vertp[5,:] = vertp[1,:]   # close outer contour
    # Inner square (hole, reversed winding)
    vertp[6,:]  = [0.25, 0.25]
    vertp[7,:]  = [0.25, 0.75]
    vertp[8,:]  = [0.75, 0.75]
    vertp[9,:]  = [0.75, 0.25]
    vertp[10,:] = vertp[6,:]  # close inner contour
    return Polygon2D(vertp, ipv, ntp, ntv)
end

"""
    ncmultisquare() -> Polygon2D

Create a non-convex multi-square cell (15 vertices): outer square, inner
hollow, and a small liquid island inside the hollow.
Vertex 5 = vertex 1, vertex 10 = vertex 6, vertex 15 = vertex 11.
"""
function ncmultisquare()
    nv = NV_DEFAULT
    ipv = zeros(Int, nv)
    vertp = zeros(Float64, nv, 2)
    ntv = 15; ntp = 15
    for i in 1:ntv; ipv[i] = i; end
    # Outer square
    vertp[1,:] = [0.0,  0.0]
    vertp[2,:] = [1.0,  0.0]
    vertp[3,:] = [1.0,  1.0]
    vertp[4,:] = [0.0,  1.0]
    vertp[5,:] = vertp[1,:]   # close outer contour
    # Inner hollow (reversed winding)
    vertp[6,:]  = [0.25, 0.25]
    vertp[7,:]  = [0.25, 0.75]
    vertp[8,:]  = [0.75, 0.75]
    vertp[9,:]  = [0.75, 0.25]
    vertp[10,:] = vertp[6,:]  # close inner contour
    # Liquid island inside hollow
    vertp[11,:] = [0.4, 0.4]
    vertp[12,:] = [0.6, 0.4]
    vertp[13,:] = [0.6, 0.6]
    vertp[14,:] = [0.4, 0.6]
    vertp[15,:] = vertp[11,:] # close island contour
    return Polygon2D(vertp, ipv, ntp, ntv)
end
