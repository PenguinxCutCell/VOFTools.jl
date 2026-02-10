# ---------------------------------------------------------------------------
# Type definitions for VOFTools.jl
# ---------------------------------------------------------------------------

"""
    Polygon2D

Mutable structure representing a 2D polygon used in VOF computations.

# Fields
- `vertp::Matrix{Float64}` – vertex coordinates, size `(nv_max, 2)`.
- `ipv::Vector{Int}`       – global vertex indices forming the polygon.
- `ntp::Int`               – last global vertex index (≥ ntv after truncation).
- `ntv::Int`               – total number of active vertices.
"""
mutable struct Polygon2D
    vertp::Matrix{Float64}   # (nv_max, 2)
    ipv::Vector{Int}         # length nv_max
    ntp::Int
    ntv::Int
end

"""
    Polygon2D(vertp::AbstractMatrix, ipv::AbstractVector)

Construct a `Polygon2D` from vertex coordinates and vertex‐index list.
`ntp` and `ntv` are set to `length(ipv)`.
"""
function Polygon2D(vertp::AbstractMatrix{<:Real}, ipv::AbstractVector{<:Integer})
    nv = length(ipv)
    # Copy into mutable storage with some headroom for truncation
    nv_max = max(size(vertp, 1), 2 * nv)
    v = zeros(Float64, nv_max, 2)
    idx = zeros(Int, nv_max)
    for i in 1:size(vertp, 1)
        v[i, 1] = vertp[i, 1]
        v[i, 2] = vertp[i, 2]
    end
    for i in 1:nv
        idx[i] = ipv[i]
    end
    return Polygon2D(v, idx, nv, nv)
end

# ------------------------------------------------------------------

"""
    Polyhedron3D

Mutable structure representing a 3D polyhedron used in VOF computations.

# Fields
- `vertp::Matrix{Float64}` – vertex coordinates, size `(nv_max, 3)`.
- `ipv::Matrix{Int}`       – vertex index array, size `(ns_max, nv_max)`.
  `ipv[face, local_vertex]` gives the global vertex index.
- `nipv::Vector{Int}`      – number of vertices per face.
- `xns::Vector{Float64}`   – x-components of outward unit normals per face.
- `yns::Vector{Float64}`   – y-components.
- `zns::Vector{Float64}`   – z-components.
- `nts::Int`               – total number of faces.
- `ntp::Int`               – last global vertex index.
- `ntv::Int`               – total number of active vertices.
"""
mutable struct Polyhedron3D
    vertp::Matrix{Float64}   # (nv_max, 3)
    ipv::Matrix{Int}         # (ns_max, nv_max)
    nipv::Vector{Int}        # length ns_max
    xns::Vector{Float64}     # length ns_max
    yns::Vector{Float64}
    zns::Vector{Float64}
    nts::Int
    ntp::Int
    ntv::Int
end

"""
    Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)

Construct a `Polyhedron3D` with explicit sizes.  The arrays are deep-copied.
"""
function Polyhedron3D(vertp::AbstractMatrix{<:Real},
                      ipv::AbstractMatrix{<:Integer},
                      nipv::AbstractVector{<:Integer},
                      xns::AbstractVector{<:Real},
                      yns::AbstractVector{<:Real},
                      zns::AbstractVector{<:Real},
                      nts::Integer, ntp::Integer, ntv::Integer)
    ns_max = max(size(ipv, 1), 2 * nts)
    nv_max = max(size(vertp, 1), size(ipv, 2), 2 * ntv)

    v = zeros(Float64, nv_max, 3)
    for i in 1:size(vertp, 1), j in 1:min(3, size(vertp, 2))
        v[i, j] = vertp[i, j]
    end
    ip = zeros(Int, ns_max, nv_max)
    for i in 1:size(ipv, 1), j in 1:size(ipv, 2)
        ip[i, j] = ipv[i, j]
    end
    ni = zeros(Int, ns_max)
    xn = zeros(Float64, ns_max)
    yn = zeros(Float64, ns_max)
    zn = zeros(Float64, ns_max)
    for i in 1:length(nipv)
        ni[i] = nipv[i]
    end
    for i in 1:length(xns)
        xn[i] = xns[i]
        yn[i] = yns[i]
        zn[i] = zns[i]
    end
    return Polyhedron3D(v, ip, ni, xn, yn, zn, nts, ntp, ntv)
end

# Convenience: default buffer sizes
const NS_DEFAULT = 200
const NV_DEFAULT = 240
