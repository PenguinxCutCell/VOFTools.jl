# ---------------------------------------------------------------------------
#  2D subroutines for VOFTools.jl
#  Translated from voftools.f (Version 5, January 2020)
#  Copyright (C) 2016 J. Lopez and J. Hernandez
# ---------------------------------------------------------------------------

# ============================= TOOLV2D =====================================
"""
    toolv2d(poly::Polygon2D) -> Float64
    toolv2d(ipv, ntv, vertp) -> Float64

Compute the area of a 2D polygon.
"""
function toolv2d(ipv::AbstractVector{Int}, ntv::Int,
                 vertp::AbstractMatrix{Float64})
    return _signed_area_2d(ipv, ntv, vertp) / 2.0
end

function toolv2d(poly::Polygon2D)
    return toolv2d(poly.ipv, poly.ntv, poly.vertp)
end

# ============================= CPPOL2D =====================================
"""
    cppol2d(poly::Polygon2D) -> Polygon2D

Deep-copy a 2D polygon.
"""
function cppol2d(poly::Polygon2D)
    return Polygon2D(copy(poly.vertp), copy(poly.ipv), poly.ntp, poly.ntv)
end

"""
    cppol2d!(dst, src)

Copy polygon `src` into `dst` (in-place).
"""
function cppol2d!(dst::Polygon2D, src::Polygon2D)
    dst.ntp = src.ntp
    dst.ntv = src.ntv
    for iv in 1:src.ntv
        ip = src.ipv[iv]
        dst.ipv[iv] = ip
        dst.vertp[ip, 1] = src.vertp[ip, 1]
        dst.vertp[ip, 2] = src.vertp[ip, 2]
    end
    return dst
end

# ============================= NEWPOL2D ====================================
"""
    newpol2d!(ia, ipia0, ipia1, ipv0, ntp0, ntv0, vertp0, xncut, yncut)
             -> (ntp0, ntv0)

Vertex indices arrangement of the truncated polygon (in-place).
Returns updated `(ntp0, ntv0)`.
"""
function _newpol2d_impl!(ia::Vector{Int}, ipia0::Vector{Int}, ipia1::Vector{Int},
                         ipv0::Vector{Int}, ntp0::Ref{Int}, ntv0::Ref{Int},
                         vertp0::Matrix{Float64}, xncut::Vector{Float64},
                         yncut::Vector{Float64}, ipv1::Vector{Int})
    ntv1 = 0
    icut = 0

    for iv in 1:ntv0[]
        ip = ipv0[iv]
        iv2 = (iv == ntv0[]) ? 1 : iv + 1
        ip2 = ipv0[iv2]
        if ia[ip] == 1
            ntv1 += 1
            ipv1[ntv1] = ipv0[iv]
        end
        if ia[ip] != ia[ip2]
            icut += 1
            ntp0[] += 1
            ntv1 += 1
            ipv1[ntv1] = ntp0[]
            ia[ntp0[]] = 0
            if ia[ip2] == 0
                ipia0[icut] = ip2
                ipia1[icut] = ip
            else
                ipia0[icut] = ip
                ipia1[icut] = ip2
            end
            xv = vertp0[ip2, 1] - vertp0[ip, 1]
            yv = vertp0[ip2, 2] - vertp0[ip, 2]
            rmod = sqrt(xv^2 + yv^2)
            xncut[icut] = yv / rmod
            yncut[icut] = -xv / rmod
        end
    end

    ntv0[] = ntv1
    for iv in 1:ntv1
        ipv0[iv] = ipv1[iv]
    end
    return
end

function newpol2d!(ia::Vector{Int}, ipia0::Vector{Int}, ipia1::Vector{Int},
                   ipv0::Vector{Int}, ntp0::Ref{Int}, ntv0::Ref{Int},
                   vertp0::Matrix{Float64}, xncut::Vector{Float64},
                   yncut::Vector{Float64})
    ipv1 = zeros(Int, length(ipv0))
    return _newpol2d_impl!(ia, ipia0, ipia1, ipv0, ntp0, ntv0, vertp0, xncut, yncut, ipv1)
end

mutable struct _Inte2DWork
    ia::Vector{Int}
    phiv::Vector{Float64}
    ipia0::Vector{Int}
    ipia1::Vector{Int}
    xncut::Vector{Float64}
    yncut::Vector{Float64}
    ntp_ref::Base.RefValue{Int}
    ntv_ref::Base.RefValue{Int}
    ipv1::Vector{Int}
end

function _Inte2DWork(poly::Polygon2D)
    nv_store = size(poly.vertp, 1)
    nv_max = length(poly.ipv)
    return _Inte2DWork(
        zeros(Int, nv_store),
        zeros(nv_store),
        zeros(Int, 2),
        zeros(Int, 2),
        zeros(2),
        zeros(2),
        Ref(0),
        Ref(0),
        zeros(Int, nv_max),
    )
end

function _get_inte2d_work(poly::Polygon2D)
    key = (length(poly.ipv), size(poly.vertp, 1))
    tls = task_local_storage()
    cache = get!(tls, :_voftools_inte2d_work) do
        Dict{NTuple{2, Int}, _Inte2DWork}()
    end::Dict{NTuple{2, Int}, _Inte2DWork}
    return get!(cache, key) do
        _Inte2DWork(poly)
    end::_Inte2DWork
end

# ============================= INTE2D ======================================
"""
    inte2d!(poly::Polygon2D, c, xnc, ync) -> (icontn, icontp)

Truncate `poly` (in-place) by the line ``\\mathbf{x} \\cdot \\mathbf{n}_c + c = 0``.
Returns `(icontn, icontp)`.
"""
function inte2d!(poly::Polygon2D, c::Float64, xnc::Float64, ync::Float64)
    ipv0   = poly.ipv
    vertp0 = poly.vertp
    ntp0   = poly.ntp
    ntv0   = poly.ntv
    work = _get_inte2d_work(poly)::_Inte2DWork
    ia = work.ia
    phiv = work.phiv
    ipia0 = work.ipia0
    ipia1 = work.ipia1
    xncut = work.xncut
    yncut = work.yncut
    ntp_ref = work.ntp_ref
    ntv_ref = work.ntv_ref
    ipv1 = work.ipv1

    icontp = 0; icontn = 0
    tolp = 1.0e-14

    fill!(ia, 0)

    for iv in 1:ntv0
        ip = ipv0[iv]
        phiv[ip] = xnc * vertp0[ip, 1] + ync * vertp0[ip, 2] + c
        if phiv[ip] > 0.0
            ia[ip] = 1
            icontp += 1
        else
            ia[ip] = 0
            icontn += 1
        end
    end

    if icontp != 0 && icontn != 0
        fill!(ipia0, 0)
        fill!(ipia1, 0)
        fill!(xncut, 0.0)
        fill!(yncut, 0.0)
        ntp_ref[] = ntp0
        ntv_ref[] = ntv0

        _newpol2d_impl!(ia, ipia0, ipia1, ipv0, ntp_ref, ntv_ref,
                        vertp0, xncut, yncut, ipv1)
        ntp0 = ntp_ref[]; ntv0 = ntv_ref[]

        # Position of the new vertices
        for idx in 1:2
            ip  = ntp0 - 2 + idx
            ip0 = ipia0[idx]
            ip1 = ipia1[idx]
            if abs(phiv[ip1] - phiv[ip0]) < tolp
                vertp0[ip, 1] = (vertp0[ip0, 1] + vertp0[ip1, 1]) / 2.0
                vertp0[ip, 2] = (vertp0[ip0, 2] + vertp0[ip1, 2]) / 2.0
            else
                frac = phiv[ip0] / (phiv[ip1] - phiv[ip0])
                vertp0[ip, 1] = vertp0[ip0, 1] - frac * (vertp0[ip1, 1] - vertp0[ip0, 1])
                vertp0[ip, 2] = vertp0[ip0, 2] - frac * (vertp0[ip1, 2] - vertp0[ip0, 2])
            end
        end
    end

    poly.ntp = ntp0
    poly.ntv = ntv0
    return icontn, icontp
end

# ============================= RESTORE2D ===================================
"""
    restore2d!(poly::Polygon2D)

Restore the structure of a polygon after truncation.
"""
function restore2d!(poly::Polygon2D)
    tolp = 1.0e-16
    poly0 = cppol2d(poly)
    ipv0 = poly0.ipv; vertp0 = poly0.vertp
    ntp0 = poly0.ntp; ntv0 = poly0.ntv

    ivt = 0
    for iv in 1:ntv0
        ip  = ipv0[iv]
        iv0 = (iv == 1) ? ntv0 : iv - 1
        ip0 = ipv0[iv0]
        dmod = sqrt((vertp0[ip, 1] - vertp0[ip0, 1])^2 +
                     (vertp0[ip, 2] - vertp0[ip0, 2])^2)
        if dmod > tolp
            ivt += 1
            poly.ipv[ivt] = ivt
            poly.vertp[ivt, 1] = vertp0[ip, 1]
            poly.vertp[ivt, 2] = vertp0[ip, 2]
        end
    end
    poly.ntv = ivt
    poly.ntp = ivt
    return poly
end

# ============================= DIST2D ======================================
"""
    dist2d(x, y, xp, yp) -> d

Compute the exact distance from point `(xp, yp)` to a segment defined by
the two endpoints in `x[1:2]`, `y[1:2]`.
"""
function dist2d(x::AbstractVector{Float64}, y::AbstractVector{Float64},
                xp::Float64, yp::Float64)
    xnt = x[2] - x[1]
    ynt = y[2] - y[1]
    vmod = sqrt(xnt^2 + ynt^2)
    xnt /= vmod; ynt /= vmod
    xn = -ynt; yn = xnt
    c1 = -(xnt * x[1] + ynt * y[1])
    c2 =  (xnt * x[2] + ynt * y[2])
    phi1 =  xnt * xp + ynt * yp + c1
    phi2 = -xnt * xp - ynt * yp + c2
    if phi1 >= 0.0 && phi2 >= 0.0
        return abs(xn * x[1] + yn * y[1] - (xn * xp + yn * yp))
    elseif phi1 <= 0.0
        return sqrt((xp - x[1])^2 + (yp - y[1])^2)
    else
        return sqrt((xp - x[2])^2 + (yp - y[2])^2)
    end
end

# ============================= ENFORV2DSZ ==================================
"""
    enforv2dsz(dx, dy, v, vertp, xnc, ync) -> c

Solve the local volume enforcement problem for rectangular cells
using the method of Scardovelli & Zaleski (JCP 164, 2000).
"""
function enforv2dsz(dx::Float64, dy::Float64, v::Float64,
                    vertp::AbstractMatrix{Float64},
                    xnc::Float64, ync::Float64)
    cmin = 1.0e14; cmax = -1.0e14
    vt = dx * dy
    vol = v / vt
    imin = 1; imax = 1

    for i in 1:4
        ci = -(vertp[i, 1] * xnc + vertp[i, 2] * ync)
        if ci <= cmin; cmin = ci; imin = i; end
        if ci >= cmax; cmax = ci; imax = i; end
    end

    if (v / vt) > 0.5
        vol = 1.0 - v / vt
    end

    sn = abs(xnc) + abs(ync)
    xm = xnc / sn; ym = ync / sn
    xmi = xm * dx; ymi = ym * dy
    sn2 = abs(xmi) + abs(ymi)
    xm = abs(xmi) / sn2; ym = abs(ymi) / sn2

    m1 = min(xm, ym)
    m = m1
    v1 = m / (2.0 * (1.0 - m))

    if vol >= 0.0 && vol < v1
        alpha = sqrt(2.0 * m * (1.0 - m) * vol)
    else
        alpha = vol * (1.0 - m) + m / 2.0
    end

    if (v / vt) <= 0.5
        return cmin + alpha * abs(cmax - cmin)
    else
        return cmax - alpha * abs(cmax - cmin)
    end
end

# ============================= ENFORV2D ====================================
# Small in-place insertion sort to avoid closure allocations in hot paths.
function _sort_listv_desc_by_phiv!(listv::Vector{Int}, phiv::Vector{Float64}, n::Int)
    for i in 2:n
        key = listv[i]
        keyv = phiv[key]
        j = i - 1
        while j >= 1 && phiv[listv[j]] < keyv
            listv[j + 1] = listv[j]
            j -= 1
        end
        listv[j + 1] = key
    end
    return listv
end

mutable struct _Enforv2DWork
    poly0::Polygon2D
    poly1::Polygon2D
    phiv::Vector{Float64}
    ia::Vector{Int}
    listv::Vector{Int}
    ipia0::Vector{Int}
    ipia1::Vector{Int}
    xncut::Vector{Float64}
    yncut::Vector{Float64}
    ntp_ref::Base.RefValue{Int}
    ntv_ref::Base.RefValue{Int}
    ipv1::Vector{Int}
end

function _Enforv2DWork(poly::Polygon2D)
    nv_store = size(poly.vertp, 1)
    nv_max = length(poly.ipv)
    return _Enforv2DWork(
        cppol2d(poly),
        cppol2d(poly),
        zeros(nv_store),
        zeros(Int, nv_store),
        zeros(Int, nv_store),
        zeros(Int, 2),
        zeros(Int, 2),
        zeros(2),
        zeros(2),
        Ref(0),
        Ref(0),
        zeros(Int, nv_max),
    )
end

function _get_enforv2d_work(poly::Polygon2D)
    key = (length(poly.ipv), size(poly.vertp, 1))
    tls = task_local_storage()
    cache = get!(tls, :_voftools_enforv2d_work) do
        Dict{NTuple{2, Int}, _Enforv2DWork}()
    end::Dict{NTuple{2, Int}, _Enforv2DWork}
    return get!(cache, key) do
        _Enforv2DWork(poly)
    end::_Enforv2DWork
end

"""
    enforv2d!(poly::Polygon2D, v, vt, xnc, ync) -> c

Solve the local volume enforcement problem in 2D using the CIBRAVE method.
`poly` is used as scratch. Returns the constant `c` of the line `x·nc + c = 0`.
"""
function enforv2d!(poly::Polygon2D, v::Float64, vt::Float64,
                   xnc::Float64, ync::Float64)
    ipv   = poly.ipv
    vertp = poly.vertp
    ntp   = poly.ntp; ntv = poly.ntv
    tolc = 1.0e-12

    if ntp > ntv
        restore2d!(poly)
        ntp = poly.ntp; ntv = poly.ntv
    end

    work = _get_enforv2d_work(poly)::_Enforv2DWork
    poly0 = work.poly0
    phiv = work.phiv
    ia = work.ia
    listv = work.listv
    ipia0 = work.ipia0
    ipia1 = work.ipia1
    xncut = work.xncut
    yncut = work.yncut
    ntp_ref = work.ntp_ref
    ntv_ref = work.ntv_ref
    ipv1 = work.ipv1

    vaux = v

    # Compute phi and ordered list
    fill!(ia, 0)
    for iv in 1:ntv
        phiv[iv] = xnc * vertp[iv, 1] + ync * vertp[iv, 2]
        listv[iv] = iv
    end
    _sort_listv_desc_by_phiv!(listv, phiv, ntv)

    invert = 0
    xncor = xnc; yncor = ync
    cmax_loc = phiv[listv[1]]
    cmin_loc = phiv[listv[ntp]]
    imin = 1; imax = ntp
    vmin_val = 0.0; vmax_val = vt

    xnc_local = xnc; ync_local = ync

    imaxlold = ntp + 1
    local iminl, imaxl

    while true
        # Obtain tentative solution bracketing by interpolation
        phiint = phiv[listv[imin]] - (phiv[listv[imin]] - phiv[listv[imax]]) *
                 (v - vmin_val) / (vmax_val - vmin_val)
        imaxl = 0; iminl = 0
        for ip in (imin + 1):imax
            if phiv[listv[ip]] < phiint
                imaxl = ip
                iminl = ip - 1
                break
            end
        end
        if imaxl == 0
            imaxl = imax; iminl = imax - 1
        end

        cmax_loc = phiv[listv[iminl]]
        cmin_loc = phiv[listv[imaxl]]

        if (ntp - imaxl) < (iminl - 1)
            invert = 1
            caux = cmin_loc
            cmin_loc = -cmax_loc
            cmax_loc = -caux
            vaux = vt - v
            xnc_local = -xncor
            ync_local = -yncor
        else
            invert = 0
            vaux = v
            xnc_local = xncor
            ync_local = yncor
        end

        for i in 1:ntp
            if i <= iminl
                ia[listv[i]] = 1 - invert
            else
                ia[listv[i]] = invert
            end
        end

        # Copy polygon to working copy
        cppol2d!(poly0, poly)
        ipv0 = poly0.ipv; vertp0 = poly0.vertp
        ntp0 = poly0.ntp; ntv0 = poly0.ntv

        # Construction of the new polygon
        fill!(ipia0, 0)
        fill!(ipia1, 0)
        fill!(xncut, 0.0)
        fill!(yncut, 0.0)
        ntp_ref[] = ntp0
        ntv_ref[] = ntv0

        ntpini = ntp0
        _newpol2d_impl!(ia, ipia0, ipia1, ipv0, ntp_ref, ntv_ref,
                        vertp0, xncut, yncut, ipv1)
        ntp0 = ntp_ref[]; ntv0 = ntv_ref[]
        # Zero out new vertex coordinates and ia flags (v3.2 fix)
        for ip in (ntpini + 1):ntp0
            vertp0[ip, 1] = 0.0
            vertp0[ip, 2] = 0.0
            ia[ip] = 0
        end

        # Determine orientation
        if (xncut[1] * ync_local - yncut[1] * xnc_local) > 0.0
            i1 = 1; i2 = 2
        else
            i1 = 2; i2 = 1
        end

        xncs = -xnc_local; yncs = -ync_local
        ipf1 = ipia1[i1]; ipf2 = ipia1[i2]
        cf1 = -(vertp0[ipf1, 1] * xncut[i1] + vertp0[ipf1, 2] * yncut[i1])
        cf2 = -(vertp0[ipf2, 1] * xncut[i2] + vertp0[ipf2, 2] * yncut[i2])
        bet1 =  1.0 / (-yncut[i1] * xncs + xncut[i1] * yncs)
        bet2 = -1.0 / (-yncut[i2] * xncs + xncut[i2] * yncs)

        # Coefficients of the analytical equation: c2*x^2 + c1*x + c0 = 0
        c2 = (yncut[i2] * xncut[i1] - xncut[i2] * yncut[i1]) * bet1 * bet2
        c1_coef = -2.0 * (cf2 * bet2 + cf1 * bet1)
        c0 = 0.0
        ih = div(ntv0 - 2, 2)
        ip0_idx = ipv0[1]
        x1 = vertp0[ip0_idx, 1] * ia[ip0_idx]
        y1 = vertp0[ip0_idx, 2] * ia[ip0_idx]
        for i in 2:(ih + 1)
            iv  = 2 * i
            ip  = ipv0[iv]
            ip1 = ipv0[iv - 1]
            ip2 = ipv0[iv - 2]
            xv1 = vertp0[ip1, 1] * ia[ip1] - x1
            yv1 = vertp0[ip1, 2] * ia[ip1] - y1
            xv2 = vertp0[ip, 1] * ia[ip] - vertp0[ip2, 1] * ia[ip2]
            yv2 = vertp0[ip, 2] * ia[ip] - vertp0[ip2, 2] * ia[ip2]
            c0 += xv1 * yv2 - yv1 * xv2
        end
        if 2 * (ih + 1) < ntv0
            ipn0 = ipv0[ntv0 - 1]
            ipn  = ipv0[ntv0]
            xv1 = vertp0[ipn, 1] * ia[ipn] - x1
            yv1 = vertp0[ipn, 2] * ia[ipn] - y1
            xv2 = x1 - vertp0[ipn0, 1] * ia[ipn0]
            yv2 = y1 - vertp0[ipn0, 2] * ia[ipn0]
            c0 += xv1 * yv2 - yv1 * xv2
        end
        c0 -= (xncs * vertp0[ipf1, 1] + yncs * vertp0[ipf1, 2]) * cf1 * bet1 +
              (xncs * vertp0[ipf2, 1] + yncs * vertp0[ipf2, 2]) * cf2 * bet2 +
              2.0 * vaux
        c3 = 0.0

        vmaxl = (c2 * cmin_loc^2 + c1_coef * cmin_loc + c0 + 2.0 * vaux) / 2.0
        vminl = (c2 * cmax_loc^2 + c1_coef * cmax_loc + c0 + 2.0 * vaux) / 2.0
        if invert == 1
            vminll = vt - vmaxl
            vmaxl = vt - vminl
            vminl = vminll
        end
        sv = (vminl - v) * (vmaxl - v)
        if sv > 0.0 && (imax - imin) > 1 && imaxlold != imaxl
            if vmaxl > v
                vmax_val = vminl
                imax = iminl
            else
                vmin_val = vmaxl
                imin = imaxl
            end
            imaxlold = imaxl
            continue
        end

        c_sol = eqsol3d(c0, c1_coef, c2, c3, cmin_loc, cmax_loc)
        dsol = abs(c_sol - cmin_loc) + abs(c_sol - cmax_loc)
        dref = abs(cmax_loc - cmin_loc)
        if dref > 0 && (dsol / dref) > (1.0 + tolc)
            csoln, isol = newton3d(c0, c1_coef, c2, c3, cmin_loc, cmax_loc, c_sol)
            dsoln = abs(csoln - cmin_loc) + abs(csoln - cmax_loc)
            if dref > 0 && (dsoln / dref) > (1.0 + tolc)
                c_sol = csoln
            end
        end

        if invert == 0
            return -c_sol
        else
            return c_sol
        end
    end
end

"""
    enforv2d(poly::Polygon2D, v, vt, xnc, ync) -> c

Non-mutating version: uses a task-local scratch copy, then solves.
"""
function enforv2d(poly::Polygon2D, v::Float64, vt::Float64,
                  xnc::Float64, ync::Float64)
    work = _get_enforv2d_work(poly)::_Enforv2DWork
    cppol2d!(work.poly1, poly)
    return enforv2d!(work.poly1, v, vt, xnc, ync)
end

# ============================= INITF2D =====================================
mutable struct _Initf2DWork
    phiv::Vector{Float64}
    ia::Vector{Int}
    poly2::Polygon2D
    poly1::Polygon2D
    ipia0::Vector{Int}
    ipia1::Vector{Int}
    xncut::Vector{Float64}
    yncut::Vector{Float64}
    ntp_ref::Base.RefValue{Int}
    ntv_ref::Base.RefValue{Int}
    ipv1::Vector{Int}
end

function _Initf2DWork(poly::Polygon2D)
    nv_store = size(poly.vertp, 1)
    nv_max = length(poly.ipv)
    return _Initf2DWork(
        zeros(nv_store),
        zeros(Int, nv_store),
        cppol2d(poly),
        cppol2d(poly),
        zeros(Int, nv_max),
        zeros(Int, nv_max),
        zeros(nv_max),
        zeros(nv_max),
        Ref(0),
        Ref(0),
        zeros(Int, nv_max),
    )
end

function _get_initf2d_work(poly::Polygon2D)
    key = (length(poly.ipv), size(poly.vertp, 1))
    tls = task_local_storage()
    cache = get!(tls, :_voftools_initf2d_work) do
        Dict{NTuple{2, Int}, _Initf2DWork}()
    end::Dict{NTuple{2, Int}, _Initf2DWork}
    return get!(cache, key) do
        _Initf2DWork(poly)
    end::_Initf2DWork
end

"""
    initf2d(func2d, poly::Polygon2D; nc::Int=10, tol::Float64=10.0) -> Float64

Initialize the material area fraction in a polygonal cell.

`func2d(x, y)` is a user-supplied implicit function defining the interface:
- `> 0` inside the material body
- `< 0` outside
- `= 0` on the interface

`nc` is the number of sub-cells along each coordinate axis.
`tol` is a tolerance for near-interface vertices.

Returns the area fraction `vf ∈ [0, 1]`.
"""
function initf2d(func2d::F, poly::Polygon2D;
                 nc::Int=10, tol::Float64=10.0) where {F}
    ipv   = poly.ipv
    vertp = poly.vertp
    ntp   = poly.ntp; ntv = poly.ntv
    work = _get_initf2d_work(poly)::_Initf2DWork

    # Coordinate extremes and vertex tagging
    xmin = 1.0e20; xmax = -1.0e20
    ymin = 1.0e20; ymax = -1.0e20
    icontp = 0; icontn = 0

    phiv = work.phiv
    ia = work.ia

    for iv in 1:ntv
        ip = ipv[iv]
        xp = vertp[ip, 1]; yp = vertp[ip, 2]
        xmin = min(xmin, xp); xmax = max(xmax, xp)
        ymin = min(ymin, yp); ymax = max(ymax, yp)
        phiv[ip] = func2d(xp, yp)
        if phiv[ip] >= 0.0
            ia[ip] = 1; icontp += 1
        else
            ia[ip] = 0; icontn += 1
        end
    end

    dx = xmax - xmin; dy = ymax - ymin

    # Check if any vertex is near-interface
    iphi = 0
    phimin = 10.0 * max(dx, dy)
    for i in 1:ntv
        phimin = min(phimin, abs(phiv[i]))
    end
    if phimin < tol * dx
        iphi = 1
    end

    if icontp == ntv && iphi == 0
        return 1.0
    elseif icontn == ntv && iphi == 0
        return 0.0
    end

    # Total area of the original polygon
    volt = toolv2d(poly)
    ddx = dx / nc; ddy = dy / nc

    vf = 0.0
    poly2 = work.poly2
    poly1 = work.poly1
    ipia0 = work.ipia0
    ipia1 = work.ipia1
    xncut = work.xncut
    yncut = work.yncut
    ntp_ref = work.ntp_ref
    ntv_ref = work.ntv_ref
    ipv1 = work.ipv1

    for ic in 1:nc
        xc = xmin + (ic - 1) * ddx
        # Copy polygon -> poly2
        cppol2d!(poly2, poly)
        # Truncate by x-lines
        if ic > 1
            cx1 = -xc
            inte2d!(poly2, cx1, 1.0, 0.0)
        end
        cx2 = xc + ddx
        icontn_x, icontp_x = inte2d!(poly2, cx2, -1.0, 0.0)

        for jc in 1:nc
            yc = ymin + (jc - 1) * ddy
            # Copy poly2 -> poly1
            cppol2d!(poly1, poly2)
            if jc > 1
                cy1 = -yc
                icontn_y, icontp_y = inte2d!(poly1, cy1, 0.0, 1.0)
            end
            if (jc == 1) || icontp_y != 0
                cy2 = yc + ddy
                icontn_y2, icontp_y2 = inte2d!(poly1, cy2, 0.0, -1.0)
                if icontp_y2 != 0
                    # Evaluate the interface function at subcell vertices
                    icontp_sub = 0; icontn_sub = 0
                    for iv in 1:poly1.ntv
                        ip = poly1.ipv[iv]
                        xp = poly1.vertp[ip, 1]; yp = poly1.vertp[ip, 2]
                        phiv[ip] = func2d(xp, yp)
                        if phiv[ip] >= 0.0
                            ia[ip] = 1; icontp_sub += 1
                        else
                            ia[ip] = 0; icontn_sub += 1
                        end
                    end

                    if icontn_sub == 0
                        volf = toolv2d(poly1)
                        vf += volf
                    elseif icontn_sub > 0 && icontp_sub > 0
                        ntpini = poly1.ntp
                        # Truncate by the interface using newpol2d!
                        fill!(ipia0, 0)
                        fill!(ipia1, 0)
                        fill!(xncut, 0.0)
                        fill!(yncut, 0.0)
                        ntp_ref[] = poly1.ntp
                        ntv_ref[] = poly1.ntv
                        _newpol2d_impl!(ia, ipia0, ipia1, poly1.ipv, ntp_ref, ntv_ref,
                                        poly1.vertp, xncut, yncut, ipv1)
                        poly1.ntp = ntp_ref[]; poly1.ntv = ntv_ref[]

                        # Location of the new intersection points
                        for ip in (ntpini + 1):(ntpini + 2)
                            ip0 = ipia0[ip - ntpini]
                            ip1 = ipia1[ip - ntpini]
                            poly1.vertp[ip, 1] = poly1.vertp[ip0, 1] -
                                phiv[ip0] * (poly1.vertp[ip1, 1] - poly1.vertp[ip0, 1]) /
                                (phiv[ip1] - phiv[ip0])
                            poly1.vertp[ip, 2] = poly1.vertp[ip0, 2] -
                                phiv[ip0] * (poly1.vertp[ip1, 2] - poly1.vertp[ip0, 2]) /
                                (phiv[ip1] - phiv[ip0])
                        end
                        volf = toolv2d(poly1)
                        vf += volf
                    end
                end
            end
        end
    end
    vf /= volt
    return vf
end
