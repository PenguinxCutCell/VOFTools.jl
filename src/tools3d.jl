# ---------------------------------------------------------------------------
#  3D subroutines for VOFTools.jl
#  Translated from voftools.f (Version 5, January 2020)
#  Copyright (C) 2016 J. Lopez and J. Hernandez
# ---------------------------------------------------------------------------

# ============================= TOOLV3D =====================================
"""
    toolv3d(poly::Polyhedron3D) -> Float64
    toolv3d(ipv, nipv, nts, verti, xns, yns, zns) -> Float64

Compute the volume of a convex polyhedron.
"""
function toolv3d(ipv::AbstractMatrix{Int}, nipv::AbstractVector{Int},
                 nts::Int, verti::AbstractMatrix{Float64},
                 xns::AbstractVector{Float64}, yns::AbstractVector{Float64},
                 zns::AbstractVector{Float64})
    sums = 0.0
    for is in 1:nts
        if nipv[is] > 0
            sump, dnmax = _face_area_contribution(is, nipv[is], ipv, verti,
                                                   xns, yns, zns)
            if dnmax != 0.0
                cs_face = xns[is] * verti[ipv[is, 1], 1] +
                          yns[is] * verti[ipv[is, 1], 2] +
                          zns[is] * verti[ipv[is, 1], 3]
                sums += cs_face * (sump / dnmax)
            end
        end
    end
    return sums / 6.0
end

function toolv3d(poly::Polyhedron3D)
    return toolv3d(poly.ipv, poly.nipv, poly.nts, poly.vertp,
                   poly.xns, poly.yns, poly.zns)
end

# ============================= CPPOL3D =====================================
"""
    cppol3d(poly::Polyhedron3D) -> Polyhedron3D

Deep-copy a polyhedron.
"""
function cppol3d(poly::Polyhedron3D)
    return Polyhedron3D(copy(poly.vertp), copy(poly.ipv), copy(poly.nipv),
                        copy(poly.xns), copy(poly.yns), copy(poly.zns),
                        poly.nts, poly.ntp, poly.ntv)
end

"""
    cppol3d!(dst, src)

Copy all fields of `src` polyhedron into `dst` (in-place).
"""
function cppol3d!(dst::Polyhedron3D, src::Polyhedron3D)
    dst.nts = src.nts
    dst.ntv = src.ntv
    dst.ntp = src.ntp
    ntp = src.ntp
    for ip in 1:ntp, j in 1:3
        dst.vertp[ip, j] = src.vertp[ip, j]
    end
    for is in 1:src.nts
        dst.xns[is] = src.xns[is]
        dst.yns[is] = src.yns[is]
        dst.zns[is] = src.zns[is]
        dst.nipv[is] = src.nipv[is]
        for j in 1:src.nipv[is]
            dst.ipv[is, j] = src.ipv[is, j]
        end
    end
    return dst
end

# Low-level version operating on raw arrays (mirrors Fortran signature)
function cppol3d!(cs, cs0, ipv, ipv0, nipv, nipv0,
                  ntp, ntp0, nts, nts0, ntv, ntv0,
                  verti, verti0, xns, xns0, yns, yns0, zns, zns0)
    nts_val = nts0[]
    ntv_val = ntv0[]
    ntp_val = ntp0[]
    nts[] = nts_val
    ntv[] = ntv_val
    ntp[] = ntp_val
    for ip in 1:ntp_val, j in 1:3
        verti[ip, j] = verti0[ip, j]
    end
    for is in 1:nts_val
        xns[is] = xns0[is]
        yns[is] = yns0[is]
        zns[is] = zns0[is]
        nipv[is] = nipv0[is]
        cs[is] = cs0[is]
        for j in 1:nipv0[is]
            ipv[is, j] = ipv0[is, j]
        end
    end
end

# ============================= NEWPOL3D ====================================
"""
    newpol3d!(ia, ipia0, ipia1, ipv0, iscut, nipv0, ntp0, nts0, ntv0,
              xnc, xns0, ync, yns0, znc, zns0)

Vertex-index arrangement of the truncated polyhedron (in-place).
Modifies `ntp0`, `nts0`, `ntv0` (all `Ref{Int}`).

Uses the v5 edge-tracking algorithm that supports multiple new faces
(disjoint truncation regions in non-convex polyhedra).
"""
function _newpol3d_impl!(ia::Vector{Int}, ipia0::Vector{Int}, ipia1::Vector{Int},
                         ipv0::Matrix{Int}, iscut::Vector{Int},
                         nipv0::Vector{Int}, ntp0::Ref{Int}, nts0::Ref{Int},
                         ntv0::Ref{Int}, xnc::Float64, xns0::Vector{Float64},
                         ync::Float64, yns0::Vector{Float64},
                         znc::Float64, zns0::Vector{Float64},
                         ipv1::Matrix{Int}, nipv1::Vector{Int}, nedge::Vector{Int},
                         ise::Matrix{Int}, ivise::Matrix{Int},
                         ipise::Matrix{Int}, ipmark::Vector{Int})
    nv_max = size(ipv0, 2)
    fill!(nipv1, 0)
    fill!(nedge, 0)
    fill!(ipmark, 0)

    # Determination of the cut faces
    isini = 0
    niscut = 0
    nts0_val = nts0[]

    for is in 1:nts0_val
        nedge[is] = 0
        if nipv0[is] > 0
            iscut[is] = 0
            for iv in 1:nipv0[is]
                ip = ipv0[is, iv]
                iv1 = (iv == nipv0[is]) ? 1 : iv + 1
                ip1 = ipv0[is, iv1]
                if ia[ip] != ia[ip1]
                    iscut[is] = 1
                    niscut += 1
                    if isini == 0; isini = is; end
                    nedge[is] += 1
                end
            end
            if iscut[is] == 0 && ia[ipv0[is, 1]] == 0
                nipv0[is] = -nipv0[is]
            end
        end
    end

    # Disjoint regions may produce NISCUT=0 with both ICONTP and ICONTN != 0
    if niscut == 0
        nipnew = 0
        nivnew = 0
        isnew  = 0
        @goto label50
    end

    # Construction of the cut faces
    nipnew = ntp0[]
    for is in 1:nts0_val
        if iscut[is] == 1
            niv  = 0
            nint = 0
            for iv in 1:nipv0[is]
                ip = ipv0[is, iv]
                iv1 = (iv > nipv0[is] - 1) ? 1 : iv + 1
                ip1 = ipv0[is, iv1]
                if ia[ip] == 1
                    niv += 1
                    ipv1[is, niv] = ipv0[is, iv]
                end
                if ia[ip] != ia[ip1]
                    nint += 1
                    niv += 1
                    if ia[ip] == 1
                        ip1i = ip;  ip0i = ip1; itype = 2
                    else
                        ip1i = ip1; ip0i = ip;  itype = 1
                    end
                    # Search for a shared edge in previously processed faces
                    found_shared = false
                    for is1 in 1:is-1
                        for ie in 1:nedge[is1]
                            ipnew_val = ise[is1, ie]
                            ip0n = ipia0[ipnew_val]
                            ip1n = ipia1[ipnew_val]
                            if ip0n == ip0i && ip1n == ip1i
                                ise[is, nint]          = ipnew_val
                                ipv1[is, niv]          = ipnew_val
                                ivise[is, ipnew_val]   = niv
                                ipise[ipnew_val, itype] = is
                                found_shared = true
                                break
                            end
                        end
                        found_shared && break
                    end
                    if !found_shared
                        nipnew += 1
                        ipia0[nipnew] = ip0i
                        ipia1[nipnew] = ip1i
                        ipv1[is, niv]        = nipnew
                        ise[is, nint]        = nipnew
                        ivise[is, nipnew]    = niv
                        ipise[nipnew, itype] = is
                    end
                end
            end
            nipv1[is] = niv
        end
    end

    # Construction of the new faces (may be multiple for non-convex polyhedra)
    nivnew = nipnew - ntp0[]
    isnew  = nts0_val
    for ip in (ntp0[]+1):nipnew
        ipmark[ip] = 0
    end
    ivnewt = 0
    ipnew_cur = ntp0[] + 1

    @label label40
    ivnew  = 1
    ivnewt += 1
    isnew  += 1
    ipini  = ipnew_cur
    ipv0[isnew, ivnew] = ipnew_cur
    ipmark[ipnew_cur] = 1

    @label label20
    is_cur = ipise[ipnew_cur, 1]
    iv_cur = ivise[is_cur, ipnew_cur]
    iv1_cur = iv_cur - 1
    if iv1_cur == 0; iv1_cur = nipv1[is_cur]; end
    ipnew_cur = ipv1[is_cur, iv1_cur]
    if ipnew_cur != ipini
        ivnew  += 1
        ivnewt += 1
        ipv0[isnew, ivnew] = ipnew_cur
        ipmark[ipnew_cur] = 1
        if ivnewt == nivnew; @goto label30; end
        @goto label20
    end
    nipv0[isnew] = ivnew
    for ipnew_scan in (ntp0[]+2):nipnew
        if ipmark[ipnew_scan] == 0
            ipnew_cur = ipnew_scan
            @goto label40
        end
    end
    @label label30
    nipv0[isnew] = ivnew

    # Assign the vertices of the new truncated polyhedron
    @label label50
    niv = nivnew
    ntpmax = nipnew
    ntsmax = isnew
    for is in 1:nts0_val
        if nipv0[is] > 0
            if iscut[is] == 1
                nipv0[is] = nipv1[is]
                for iv in 1:nipv1[is]
                    ipv0[is, iv] = ipv1[is, iv]
                    if ia[ipv1[is, iv]] == 1
                        niv += 1
                        ia[ipv1[is, iv]] = -1
                    end
                end
            else
                if iscut[is] == 0 && nipv0[is] < 0
                    nipv0[is] = 0
                end
                for iv in 1:nipv0[is]
                    ntsmax = max(ntsmax, is)
                    if ia[ipv0[is, iv]] == 1
                        ntpmax = max(ntpmax, ipv0[is, iv])
                        niv += 1
                        ia[ipv0[is, iv]] = -1
                    end
                end
            end
        end
    end

    # Reset ia flags
    for ip in 1:ntp0[]
        if ia[ip] == -1
            ia[ip] = 1
        end
    end

    # Set normals for all new faces
    for is in (nts0_val+1):isnew
        xns0[is] = -xnc
        yns0[is] = -ync
        zns0[is] = -znc
    end

    ntv0[] = niv
    ntp0[] = ntpmax
    nts0[] = ntsmax
    return
end

function newpol3d!(ia::Vector{Int}, ipia0::Vector{Int}, ipia1::Vector{Int},
                   ipv0::Matrix{Int}, iscut::Vector{Int},
                   nipv0::Vector{Int}, ntp0::Ref{Int}, nts0::Ref{Int},
                   ntv0::Ref{Int}, xnc::Float64, xns0::Vector{Float64},
                   ync::Float64, yns0::Vector{Float64},
                   znc::Float64, zns0::Vector{Float64})
    ns_max = size(ipv0, 1)
    nv_max = size(ipv0, 2)
    ipv1 = zeros(Int, ns_max, nv_max)
    nipv1 = zeros(Int, ns_max)
    nedge = zeros(Int, ns_max)
    ise = zeros(Int, ns_max, nv_max)
    ivise = zeros(Int, ns_max, nv_max)
    ipise = zeros(Int, nv_max, 2)
    ipmark = zeros(Int, nv_max)
    return _newpol3d_impl!(ia, ipia0, ipia1, ipv0, iscut, nipv0, ntp0, nts0, ntv0,
                           xnc, xns0, ync, yns0, znc, zns0,
                           ipv1, nipv1, nedge, ise, ivise, ipise, ipmark)
end

# Reusable task-local workspace for inte3d! hot path.
mutable struct _Inte3DWork
    ia::Vector{Int}
    phiv::Vector{Float64}
    iscut::Vector{Int}
    ipia0::Vector{Int}
    ipia1::Vector{Int}
    ntp_ref::Base.RefValue{Int}
    nts_ref::Base.RefValue{Int}
    ntv_ref::Base.RefValue{Int}
    ipv1::Matrix{Int}
    nipv1::Vector{Int}
    nedge::Vector{Int}
    ise::Matrix{Int}
    ivise::Matrix{Int}
    ipise::Matrix{Int}
    ipmark::Vector{Int}
end

function _Inte3DWork(poly::Polyhedron3D)
    ns_max = size(poly.ipv, 1)
    nv_max = size(poly.ipv, 2)
    nv_store = size(poly.vertp, 1)
    return _Inte3DWork(
        zeros(Int, nv_store),
        zeros(nv_store),
        zeros(Int, ns_max),
        zeros(Int, nv_max),
        zeros(Int, nv_max),
        Ref(0), Ref(0), Ref(0),
        zeros(Int, ns_max, nv_max),
        zeros(Int, ns_max),
        zeros(Int, ns_max),
        zeros(Int, ns_max, nv_max),
        zeros(Int, ns_max, nv_max),
        zeros(Int, nv_max, 2),
        zeros(Int, nv_max),
    )
end

function _get_inte3d_work(poly::Polyhedron3D)
    key = (size(poly.ipv, 1), size(poly.ipv, 2), size(poly.vertp, 1))
    tls = task_local_storage()
    cache = get!(tls, :_voftools_inte3d_work) do
        Dict{NTuple{3, Int}, _Inte3DWork}()
    end::Dict{NTuple{3, Int}, _Inte3DWork}
    return get!(cache, key) do
        _Inte3DWork(poly)
    end::_Inte3DWork
end

# ============================= INTE3D ======================================
"""
    inte3d!(poly::Polyhedron3D, c, xnc, ync, znc) -> (icontn, icontp)

Truncate `poly` (in-place) by the plane ``\\mathbf{x} \\cdot \\mathbf{n}_c + c = 0``.
Returns `(icontn, icontp)`.
"""
function inte3d!(poly::Polyhedron3D, c::Float64,
                 xnc::Float64, ync::Float64, znc::Float64)
    ipv0  = poly.ipv
    nipv0 = poly.nipv
    vertp0 = poly.vertp
    xns0  = poly.xns
    yns0  = poly.yns
    zns0  = poly.zns
    nts0  = poly.nts
    ntp0  = poly.ntp
    ntv0  = poly.ntv
    work = _get_inte3d_work(poly)::_Inte3DWork
    ia = work.ia
    phiv = work.phiv
    iscut = work.iscut
    ipia0 = work.ipia0
    ipia1 = work.ipia1
    ntp_ref = work.ntp_ref
    nts_ref = work.nts_ref
    ntv_ref = work.ntv_ref
    ipv1 = work.ipv1
    nipv1 = work.nipv1
    nedge = work.nedge
    ise = work.ise
    ivise = work.ivise
    ipise = work.ipise
    ipmark = work.ipmark

    icontp = 0
    icontn = 0
    sump = 0.0
    sumn = 0.0

    fill!(ia, -1)

    # Distance function and values of IA
    for is in 1:nts0
        for iv in 1:nipv0[is]
            ip = ipv0[is, iv]
            if ia[ip] == -1
                phiv[ip] = xnc * vertp0[ip, 1] + ync * vertp0[ip, 2] +
                           znc * vertp0[ip, 3] + c
                if phiv[ip] > 0.0
                    ia[ip] = 1
                    icontp += 1
                    sump += phiv[ip]
                else
                    ia[ip] = 0
                    icontn += 1
                    sumn -= phiv[ip]
                end
            end
        end
    end

    if icontp != 0 && icontn != 0
        # Construction of the new polyhedron
        nts00 = nts0
        fill!(iscut, 0)
        fill!(ipia0, 0)
        fill!(ipia1, 0)
        ntp_ref[] = ntp0
        nts_ref[] = nts0
        ntv_ref[] = ntv0

        _newpol3d_impl!(ia, ipia0, ipia1, ipv0, iscut, nipv0, ntp_ref, nts_ref,
                        ntv_ref, xnc, xns0, ync, yns0, znc, zns0,
                        ipv1, nipv1, nedge, ise, ivise, ipise, ipmark)

        nts0 = nts_ref[]
        ntp0 = ntp_ref[]
        ntv0 = ntv_ref[]

        if nts0 < 0
            nts0 = -nts0
            icontn += icontp
            icontp = 0
            poly.nts = nts0
            poly.ntp = ntp0
            poly.ntv = ntv0
            return icontn, icontp
        end

        # Position of the new vertices (loop over all new faces)
        for is in (nts00+1):nts0
            for iv in 1:nipv0[is]
                ip  = ipv0[is, iv]
                ip0 = ipia0[ip]
                ip1 = ipia1[ip]
                vertp0[ip, 1] = vertp0[ip0, 1] - phiv[ip0] * (vertp0[ip1, 1] - vertp0[ip0, 1]) /
                    (phiv[ip1] - phiv[ip0])
                vertp0[ip, 2] = vertp0[ip0, 2] - phiv[ip0] * (vertp0[ip1, 2] - vertp0[ip0, 2]) /
                    (phiv[ip1] - phiv[ip0])
                vertp0[ip, 3] = vertp0[ip0, 3] - phiv[ip0] * (vertp0[ip1, 3] - vertp0[ip0, 3]) /
                    (phiv[ip1] - phiv[ip0])
            end
            # Faces with less than 3 vertices are suppressed
            if nipv0[is] < 3
                nipv0[is] = 0
            end
        end
    end

    poly.nts = nts0
    poly.ntp = ntp0
    poly.ntv = ntv0
    return icontn, icontp
end

# ============================= RESTORE3D ===================================
"""
    restore3d!(poly::Polyhedron3D)

Restore the structure of a polyhedron after truncation, removing
duplicate/coincident vertices and renumbering.
"""
function restore3d!(poly::Polyhedron3D)
    tolp = 1.0e-16
    ns_max = size(poly.ipv, 1)
    nv_max = size(poly.ipv, 2)

    # Work copy
    poly0 = cppol3d(poly)
    ipv0  = poly0.ipv
    nipv0 = poly0.nipv
    vertp0 = poly0.vertp
    xns0 = poly0.xns; yns0 = poly0.yns; zns0 = poly0.zns
    nts0 = poly0.nts; ntp0 = poly0.ntp

    ipv = poly.ipv; nipv = poly.nipv; vertp = poly.vertp
    xns = poly.xns; yns = poly.yns; zns = poly.zns

    cs0 = zeros(ns_max)
    cs  = zeros(ns_max)

    # Eliminate consecutive duplicate vertices in each face
    for is in 1:nts0
        ipv[is, 1] = ipv0[is, 1]
        ivt = 0
        for iv in 1:nipv0[is]
            ip = ipv0[is, iv]
            iv0 = (iv == 1) ? nipv0[is] : iv - 1
            ip0 = ipv0[is, iv0]
            dmod = sqrt((vertp0[ip, 1] - vertp0[ip0, 1])^2 +
                        (vertp0[ip, 2] - vertp0[ip0, 2])^2 +
                        (vertp0[ip, 3] - vertp0[ip0, 3])^2)
            if dmod > tolp
                ivt += 1
                ipv[is, ivt] = ipv0[is, iv]
            end
        end
        nipv[is] = ivt
    end

    # Copy back
    poly0_2 = cppol3d(poly)
    ipv0  = poly0_2.ipv
    nipv0 = poly0_2.nipv
    vertp0 = poly0_2.vertp
    xns0 = poly0_2.xns; yns0 = poly0_2.yns; zns0 = poly0_2.zns
    nts0 = poly0_2.nts; ntp0 = poly0_2.ntp

    # Eliminate faces with zero or only one vertex
    nts_new = 0
    for is in 1:nts0
        if nipv0[is] > 1
            nts_new += 1
            nipv[nts_new] = nipv0[is]
            for iv in 1:nipv0[is]
                ipv[nts_new, iv] = ipv0[is, iv]
            end
            xns[nts_new] = xns0[is]
            yns[nts_new] = yns0[is]
            zns[nts_new] = zns0[is]
        end
    end
    poly.nts = nts_new
    ntp = poly.ntp

    # Link coincident vertices of different faces
    for ip1 in 1:ntp-1
        for ip2 in ip1+1:ntp
            dmod = sqrt((vertp[ip1, 1] - vertp[ip2, 1])^2 +
                        (vertp[ip1, 2] - vertp[ip2, 2])^2 +
                        (vertp[ip1, 3] - vertp[ip2, 3])^2)
            if dmod <= tolp
                for is in 1:nts_new
                    for iv in 1:nipv[is]
                        if ipv[is, iv] == ip2
                            ipv[is, iv] = ip1
                        end
                    end
                end
            end
        end
    end

    # Renumber consecutively: make ntp = ntv
    ipt0_map = zeros(Int, ntp)
    ipt = 0
    for ip in 1:ntp
        ic = false
        for is in 1:nts_new
            for iv in 1:nipv[is]
                if ipv[is, iv] == ip && !ic
                    ic = true
                    ipt += 1
                    ipt0_map[ip] = ipt
                end
            end
        end
    end

    # Copy to work arrays
    poly0_3 = cppol3d(poly)
    ipv0  = poly0_3.ipv
    nipv0 = poly0_3.nipv
    vertp0 = poly0_3.vertp
    nts0 = poly0_3.nts

    poly.ntp = ipt
    poly.ntv = ipt
    for is in 1:nts0
        for iv in 1:nipv0[is]
            ip  = ipv0[is, iv]
            ip1 = ipt0_map[ip]
            ipv[is, iv] = ip1
            vertp[ip1, 1] = vertp0[ip, 1]
            vertp[ip1, 2] = vertp0[ip, 2]
            vertp[ip1, 3] = vertp0[ip, 3]
        end
    end
    return poly
end

# ============================= EQSOL3D =====================================
"""
    eqsol3d(c0, c1, c2, c3, cmin, cmax) -> csol

Solve analytically the cubic equation `c3*x^3 + c2*x^2 + c1*x + c0 = 0`,
returning the solution bracketed by `[cmin, cmax]`.
"""
function eqsol3d(c0::Float64, c1::Float64, c2::Float64, c3::Float64,
                 cmin::Float64, cmax::Float64)
    tolc  = 1.0e-12
    tolc1 = 1.0e-12
    PI = π

    d = c3; c = c2; b = c1; a = c0

    if abs(d) <= tolc && abs(c) <= tolc
        return -a / b
    end

    if abs(d) <= tolc1
        # Quadratic: c*x^2 + b*x + a = 0
        e = b^2 - 4.0 * c * a
        if e < 0.0
            q = -0.5 * b
        else
            q = -0.5 * (b + copysign(1.0, b) * sqrt(e))
        end
        csol_1 = q / c
        csol_2 = a / q
        dsol_1 = abs(csol_1 - cmin) + abs(csol_1 - cmax)
        dsol_2 = abs(csol_2 - cmin) + abs(csol_2 - cmax)
        return dsol_1 < dsol_2 ? csol_1 : csol_2
    end

    # Full cubic: x^3 + E*x^2 + F*x + G = 0
    ee = c / d
    f  = b / d
    g  = a / d
    Q = (ee^2 - 3.0 * f) / 9.0
    R = (2.0 * ee^3 - 9.0 * ee * f + 27.0 * g) / 54.0

    if R^2 < Q^3
        theta = acos(R / sqrt(Q^3)) 
        csol_1 = -2.0 * sqrt(Q) * cos(theta / 3.0)        - ee / 3.0
        csol_2 = -2.0 * sqrt(Q) * cos((theta + 2π) / 3.0) - ee / 3.0
        csol_3 = -2.0 * sqrt(Q) * cos((theta - 2π) / 3.0) - ee / 3.0
        dsol_1 = abs(csol_1 - cmin) + abs(csol_1 - cmax)
        dsol_2 = abs(csol_2 - cmin) + abs(csol_2 - cmax)
        dsol_3 = abs(csol_3 - cmin) + abs(csol_3 - cmax)
        if dsol_1 < dsol_2 && dsol_1 < dsol_3
            return csol_1
        elseif dsol_2 < dsol_1 && dsol_2 < dsol_3
            return csol_2
        else
            return csol_3
        end
    else
        P = -copysign(1.0, R) * (abs(R) + sqrt(R^2 - Q^3))^(1.0 / 3.0)
        T = abs(P) <= tolc ? 0.0 : Q / P
        return (P + T) - ee / 3.0
    end
end

# ============================= NEWTON3D ====================================
"""
    newton3d(a, b, c, d, cmin, cmax, csol_init) -> (csol, isol)

Solve `d*x^3 + c*x^2 + b*x + a = 0` using Newton-Raphson iterations
with bisection fallback, bracketed by `[cmin, cmax]`.
"""
function newton3d(a::Float64, b::Float64, c::Float64, d::Float64,
                  cmin::Float64, cmax::Float64,
                  csol_init::Float64=0.5*(cmin+cmax))
    itemax = 100
    isol = 0
    tolc = 1.0e-14
    csol = csol_init
    cl = cmin; ch = cmax
    cold = csol
    func  = d * csol^3 + c * csol^2 + b * csol + a
    dfunc = 3.0 * d * csol^2 + 2.0 * c * csol + b
    funcmin = d * cmin^3 + c * cmin^2 + b * cmin + a
    funcmax = d * cmax^3 + c * cmax^2 + b * cmax + a

    if funcmin * funcmax > 0.0
        if dfunc * funcmin > 0.0
            return cmin, -1
        else
            return cmax, 1
        end
    end

    for _ in 1:itemax
        dcsol = (dfunc != 0.0) ? func / dfunc : 0.0
        # Use bisection when Newton goes out of bounds or dfunc==0
        if dfunc == 0.0 || (cl - csol + dcsol) * (csol - dcsol - ch) < 0.0
            dcsol = 0.5 * (ch - cl)
            csol = cl + dcsol
            if abs(csol - cold) < tolc
                return csol, isol
            end
        else
            csol -= dcsol
            if abs(dcsol) < tolc
                return csol, isol
            end
        end
        func  = d * csol^3 + c * csol^2 + b * csol + a
        dfunc = 3.0 * d * csol^2 + 2.0 * c * csol + b
        if func * funcmin > 0.0
            cl = csol
        else
            ch = csol
        end
        cold = csol
    end
    return csol, isol
end

# ============================= DIST3D ======================================
"""
    dist3d(x, y, z, xp, yp, zp) -> d

Compute the exact distance from point `(xp, yp, zp)` to a polygon
with `n` vertices given by vectors `x`, `y`, `z`.
"""
function dist3d(x::AbstractVector{Float64}, y::AbstractVector{Float64},
                z::AbstractVector{Float64},
                xp::Float64, yp::Float64, zp::Float64)
    n = length(x)
    @assert length(y) == n && length(z) == n

    # Obtain normal to the polygon
    xn = 0.0; yn = 0.0; zn = 0.0; vmod = 0.0
    for i in 1:n-2
        i2 = i + 1; i3 = i2 + 1
        xv1 = x[i2] - x[i];  yv1 = y[i2] - y[i];  zv1 = z[i2] - z[i]
        xv2 = x[i3] - x[i2]; yv2 = y[i3] - y[i2]; zv2 = z[i3] - z[i2]
        xn = yv1 * zv2 - zv1 * yv2
        yn = zv1 * xv2 - xv1 * zv2
        zn = xv1 * yv2 - yv1 * xv2
        vmod = sqrt(xn^2 + yn^2 + zn^2)
        vmod > 0.0 && break
    end
    if vmod == 0.0
        return sqrt((xp - x[1])^2 + (yp - y[1])^2 + (zp - z[1])^2)
    end
    xn /= vmod; yn /= vmod; zn /= vmod
    c_plane = -(xn * x[1] + yn * y[1] + zn * z[1])

    # Edge normals, tangents, distances
    xnt = zeros(n); ynt = zeros(n); znt = zeros(n)
    xni = zeros(n); yni = zeros(n); zni = zeros(n)
    phi = zeros(n)

    for i in 1:n
        i2 = (i == n) ? 1 : i + 1
        xt = x[i2] - x[i]; yt = y[i2] - y[i]; zt = z[i2] - z[i]
        tmod = sqrt(xt^2 + yt^2 + zt^2)
        if tmod == 0.0
            phi[i] = 0.0
        else
            xnt[i] = xt / tmod; ynt[i] = yt / tmod; znt[i] = zt / tmod
            xni[i] = yn * znt[i] - zn * ynt[i]
            yni[i] = zn * xnt[i] - xn * znt[i]
            zni[i] = xn * ynt[i] - yn * xnt[i]
            ci = -(xni[i] * x[i] + yni[i] * y[i] + zni[i] * z[i])
            phi[i] = xni[i] * xp + yni[i] * yp + zni[i] * zp + ci
        end
    end

    inside = true
    d = 0.0
    for i in 1:n
        i2 = (i == n) ? 1 : i + 1
        if phi[i] < 0.0
            inside = false
            c1 = -(xnt[i] * x[i] + ynt[i] * y[i] + znt[i] * z[i])
            c2 =  (xnt[i] * x[i2] + ynt[i] * y[i2] + znt[i] * z[i2])
            phi1 = xnt[i] * xp + ynt[i] * yp + znt[i] * zp + c1
            phi2 = -(xnt[i] * xp + ynt[i] * yp + znt[i] * zp) + c2
            if phi1 >= 0.0 && phi2 >= 0.0
                t0 = xnt[i] * (xp - x[i]) + ynt[i] * (yp - y[i]) + znt[i] * (zp - z[i])
                xq = x[i] + t0 * xnt[i]
                yq = y[i] + t0 * ynt[i]
                zq = z[i] + t0 * znt[i]
                return sqrt((xp - xq)^2 + (yp - yq)^2 + (zp - zq)^2)
            elseif phi1 <= 0.0
                d = sqrt((xp - x[i])^2 + (yp - y[i])^2 + (zp - z[i])^2)
                if i != 1
                    return d
                end
            elseif phi2 <= 0.0 && phi[i2] >= 0.0
                return sqrt((xp - x[i2])^2 + (yp - y[i2])^2 + (zp - z[i2])^2)
            end
        end
    end
    if inside
        d = abs(xn * xp + yn * yp + zn * zp + c_plane)
    end
    return d
end

# ============================= ENFORV3DSZ ==================================
"""
    enforv3dsz(dx, dy, dz, v, vertp, xnc, ync, znc) -> c

Solve the local volume enforcement problem for rectangular parallelepiped
cells using the analytical method of Scardovelli & Zaleski (JCP 164, 2000).
"""
function enforv3dsz(dx::Float64, dy::Float64, dz::Float64,
                    v::Float64, vertp::AbstractMatrix{Float64},
                    xnc::Float64, ync::Float64, znc::Float64)
    tole = 1.0e-9
    cmin = 1.0e14; cmax = -1.0e14
    vt = dx * dy * dz
    vol = v / vt
    imin = 1; imax = 1

    for i in 1:8
        ci = -(vertp[i, 1] * xnc + vertp[i, 2] * ync + vertp[i, 3] * znc)
        if ci <= cmin
            cmin = ci; imin = i
        end
        if ci >= cmax
            cmax = ci; imax = i
        end
    end

    if (v / vt) <= 0.5
        vol = v / vt
    else
        vol = 1.0 - v / vt
    end

    sn = abs(xnc) + abs(ync) + abs(znc)
    xm = xnc / sn; ym = ync / sn; zm = znc / sn
    xmi = xm * dx; ymi = ym * dy; zmi = zm * dz
    sn2 = abs(xmi) + abs(ymi) + abs(zmi)
    xm = abs(xmi) / sn2; ym = abs(ymi) / sn2; zm = abs(zmi) / sn2

    m1 = min(xm, ym, zm)
    m3 = max(xm, ym, zm)
    m2 = xm + ym + zm - m1 - m3  # middle value
    m12 = m1 + m2

    v1 = m1^2 / max(6.0 * m2 * m3, tole)
    v2 = v1 + (m2 - m1) / (2.0 * m3)
    if m3 < m12
        v3 = ((3.0 * m12 - m3) * m3^2 + (m1 - 3.0 * m3) * m1^2 +
              (m2 - 3.0 * m3) * m2^2) / (6.0 * m1 * m2 * m3)
    else
        v3 = m12 / (2.0 * m3)
    end

    # Prepare cubic coefficients for regions 2-3
    a3 = 1.0; a2 = 0.0; a1 = 0.0; a0 = 0.0
    if vol >= v2 && vol < v3
        a3_raw = -1.0
        a2 = 3.0 * m12 / a3_raw
        a1 = -3.0 * (m1^2 + m2^2) / a3_raw
        a0 = (m1^3 + m2^3 - 6.0 * m1 * m2 * m3 * vol) / a3_raw
        a3 = 1.0
    elseif vol >= v3 && vol <= 0.5 && m3 < m12
        a3_raw = -2.0
        a2 = 3.0 / a3_raw
        a1 = -3.0 * (m1^2 + m2^2 + m3^2) / a3_raw
        a0 = (m1^3 + m2^3 + m3^3 - 6.0 * m1 * m2 * m3 * vol) / a3_raw
        a3 = 1.0
    end

    alpha = 0.0
    if vol >= 0.0 && vol < v1
        alpha = (6.0 * m1 * m2 * m3 * vol)^(1.0 / 3.0)
    elseif vol >= v1 && vol < v2
        alpha = 0.5 * (m1 + sqrt(m1^2 + 8.0 * m2 * m3 * (vol - v1)))
    elseif (vol >= v2 && vol < v3) || (vol >= v3 && vol <= 0.5 && m3 < m12)
        p0 = a1 / 3.0 - a2^2 / 9.0
        q0 = (a1 * a2 - 3.0 * a0) / 6.0 - a2^3 / 27.0
        theta = acos(q0 / sqrt((-p0)^3)) / 3.0
        alpha = sqrt(-p0) * (sin(theta) * sqrt(3.0) - cos(theta)) - a2 / 3.0
    else
        alpha = m3 * vol + m12 / 2.0
    end

    if (v / vt) <= 0.5
        return cmin + alpha * abs(cmax - cmin)
    else
        return cmax - alpha * abs(cmax - cmin)
    end
end

# Reusable task-local workspace for enforv3d! hot path.
mutable struct _Enforv3DWork
    poly0::Polyhedron3D
    poly1::Polyhedron3D
    cs::Vector{Float64}
    phiv::Vector{Float64}
    ia::Vector{Int}
    listv::Vector{Int}
    iscut::Vector{Int}
    ipia0::Vector{Int}
    ipia1::Vector{Int}
    ntp_ref::Base.RefValue{Int}
    nts_ref::Base.RefValue{Int}
    ntv_ref::Base.RefValue{Int}
    betxe::Matrix{Float64}
    betye::Matrix{Float64}
    betze::Matrix{Float64}
    x0::Matrix{Float64}
    y0::Matrix{Float64}
    z0::Matrix{Float64}
    sumk::Vector{Float64}
    suml::Vector{Float64}
    summ::Vector{Float64}
    ipv1::Matrix{Int}
    nipv1::Vector{Int}
    nedge::Vector{Int}
    ise::Matrix{Int}
    ivise::Matrix{Int}
    ipise::Matrix{Int}
    ipmark::Vector{Int}
end

function _Enforv3DWork(poly::Polyhedron3D)
    ns_max = size(poly.ipv, 1)
    nv_max = size(poly.ipv, 2)
    nv_store = size(poly.vertp, 1)
    return _Enforv3DWork(
        cppol3d(poly),
        cppol3d(poly),
        zeros(ns_max),
        zeros(nv_store),
        zeros(Int, nv_store),
        zeros(Int, nv_store),
        zeros(Int, ns_max),
        zeros(Int, nv_max),
        zeros(Int, nv_max),
        Ref(0), Ref(0), Ref(0),
        zeros(ns_max, nv_max),
        zeros(ns_max, nv_max),
        zeros(ns_max, nv_max),
        zeros(ns_max, nv_max),
        zeros(ns_max, nv_max),
        zeros(ns_max, nv_max),
        zeros(ns_max),
        zeros(ns_max),
        zeros(ns_max),
        zeros(Int, ns_max, nv_max),
        zeros(Int, ns_max),
        zeros(Int, ns_max),
        zeros(Int, ns_max, nv_max),
        zeros(Int, ns_max, nv_max),
        zeros(Int, nv_max, 2),
        zeros(Int, nv_max),
    )
end

function _get_enforv3d_work(poly::Polyhedron3D)
    key = (size(poly.ipv, 1), size(poly.ipv, 2), size(poly.vertp, 1))
    tls = task_local_storage()
    cache = get!(tls, :_voftools_enforv3d_work) do
        Dict{NTuple{3, Int}, _Enforv3DWork}()
    end::Dict{NTuple{3, Int}, _Enforv3DWork}
    return get!(cache, key) do
        _Enforv3DWork(poly)
    end::_Enforv3DWork
end

# ============================= ENFORV3D ====================================
"""
    enforv3d!(poly::Polyhedron3D, v, vt, xnc, ync, znc) -> c

Solve the local volume enforcement problem in 3D using the CIBRAVE method.
`poly` is used as scratch and should be a fresh copy if you need to preserve
the original. `v` is the liquid volume, `vt` is the total cell volume.
Returns the constant `c` of the plane `x·nc + c = 0`.
"""
function enforv3d!(poly::Polyhedron3D, v::Float64, vt::Float64,
                   xnc::Float64, ync::Float64, znc::Float64)
    ipv   = poly.ipv
    nipv  = poly.nipv
    vertp = poly.vertp
    xns   = poly.xns; yns = poly.yns; zns = poly.zns
    nts   = poly.nts; ntp = poly.ntp; ntv = poly.ntv
    ns_max = size(ipv, 1); nv_max = size(ipv, 2)

    if vt <= 0.0
        @warn "THE POLYHEDRON HAS NULL OR NEGATIVE VOLUME."
        return 0.0
    end

    tolc = 1.0e-12

    work = _get_enforv3d_work(poly)::_Enforv3DWork
    cs = work.cs
    phiv = work.phiv
    ia = work.ia
    listv = work.listv
    poly0 = work.poly1
    iscut = work.iscut
    ipia0 = work.ipia0
    ipia1 = work.ipia1
    ntp_ref = work.ntp_ref
    nts_ref = work.nts_ref
    ntv_ref = work.ntv_ref
    betxe = work.betxe
    betye = work.betye
    betze = work.betze
    x0 = work.x0
    y0 = work.y0
    z0 = work.z0
    sumk = work.sumk
    suml = work.suml
    summ = work.summ
    ipv1 = work.ipv1
    nipv1 = work.nipv1
    nedge = work.nedge
    ise = work.ise
    ivise = work.ivise
    ipise = work.ipise
    ipmark = work.ipmark

    # Compute face constants cs
    fill!(cs, 0.0)
    for is in 1:nts
        ip = ipv[is, 1]
        cs[is] = -(xns[is] * vertp[ip, 1] + yns[is] * vertp[ip, 2] +
                   zns[is] * vertp[ip, 3])
    end

    # Restore if previously truncated
    if ntp > ntv
        restore3d!(poly)
        nts = poly.nts; ntp = poly.ntp; ntv = poly.ntv
        # Recompute cs
        for is in 1:nts
            ip = ipv[is, 1]
            cs[is] = -(xns[is] * vertp[ip, 1] + yns[is] * vertp[ip, 2] +
                       zns[is] * vertp[ip, 3])
        end
    end

    vaux = v

    # Compute phi values and ordered list
    for iv in 1:ntv
        phiv[iv] = xnc * vertp[iv, 1] + ync * vertp[iv, 2] + znc * vertp[iv, 3]
    end

    for iv in 1:ntp
        listv[iv] = iv
    end
    sort!(view(listv, 1:ntp), lt = (i, j) -> phiv[i] > phiv[j])  # descending order

    invert = 0
    xncor = xnc; yncor = ync; zncor = znc
    imin = 1; imax = ntp
    vmin = 0.0; vmax = vt

    xnc_local = xnc; ync_local = ync; znc_local = znc
    imaxlold = ntp + 1

    local cmax_loc, cmin_loc, iminl, imaxl

    while true
        # Obtain tentative solution bracketing by interpolation
        phiint = phiv[listv[imin]] - (phiv[listv[imin]] - phiv[listv[imax]]) *
                 (v - vmin) / (vmax - vmin)
        imaxl = 0; iminl = 0
        for ip in (imin + 1):imax
            if phiv[listv[ip]] < phiint
                imaxl = ip
                iminl = ip - 1
                break
            end
        end
        if imaxl == 0 && iminl == 0
            return -phiint
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
            znc_local = -zncor
        else
            invert = 0
            vaux = v
            xnc_local = xncor
            ync_local = yncor
            znc_local = zncor
        end

        fill!(ia, 0)
        for i in 1:ntp
            if i <= iminl
                ia[listv[i]] = 1 - invert
            else
                ia[listv[i]] = invert
            end
        end

        # Copy polyhedron to working copy
        cppol3d!(poly0, poly)
        ipv0  = poly0.ipv; nipv0 = poly0.nipv; vertp0 = poly0.vertp
        xns0 = poly0.xns; yns0 = poly0.yns; zns0 = poly0.zns
        nts0 = poly0.nts; ntp0 = poly0.ntp; ntv0 = poly0.ntv

        # Construction of the new polyhedron
        nts00 = nts0
        fill!(iscut, 0)
        fill!(ipia0, 0)
        fill!(ipia1, 0)
        ntp_ref[] = ntp0; nts_ref[] = nts0; ntv_ref[] = ntv0

        _newpol3d_impl!(ia, ipia0, ipia1, ipv0, iscut, nipv0, ntp_ref, nts_ref,
                        ntv_ref, xnc_local, xns0, ync_local, yns0, znc_local, zns0,
                        ipv1, nipv1, nedge, ise, ivise, ipise, ipmark)
        nts0 = nts_ref[]; ntp0 = ntp_ref[]; ntv0 = ntv_ref[]

        if nts0 <= nts00  # disjoint regions may produce this situation
            if (imax - imaxl) > (iminl - imin)
                imaxl += 1
                iminl = imaxl - 1
            else
                iminl -= 1
                imaxl = iminl + 1
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
                znc_local = -zncor
            else
                invert = 0
                vaux = v
                xnc_local = xncor
                ync_local = yncor
                znc_local = zncor
            end
            for i in 1:ntp
                if i <= iminl
                    ia[listv[i]] = 1 - invert
                else
                    ia[listv[i]] = invert
                end
            end
            cppol3d!(poly0, poly)
            ipv0 = poly0.ipv; nipv0 = poly0.nipv; vertp0 = poly0.vertp
            xns0 = poly0.xns; yns0 = poly0.yns; zns0 = poly0.zns
            nts0 = poly0.nts; ntp0 = poly0.ntp; ntv0 = poly0.ntv
            nts00 = nts0
            fill!(iscut, 0)
            fill!(ipia0, 0)
            fill!(ipia1, 0)
            ntp_ref[] = ntp0; nts_ref[] = nts0; ntv_ref[] = ntv0
            _newpol3d_impl!(ia, ipia0, ipia1, ipv0, iscut, nipv0, ntp_ref, nts_ref,
                            ntv_ref, xnc_local, xns0, ync_local, yns0, znc_local, zns0,
                            ipv1, nipv1, nedge, ise, ivise, ipise, ipmark)
            nts0 = nts_ref[]; ntp0 = ntp_ref[]; ntv0 = ntv_ref[]
        end

        # Contributions of the new faces Γ_c (v5: loop over NTS00+1:NTS0)
        fill!(betxe, 0.0)
        fill!(betye, 0.0)
        fill!(betze, 0.0)
        fill!(x0, 0.0)
        fill!(y0, 0.0)
        fill!(z0, 0.0)

        for is in (nts00+1):nts0
            for iv in 1:nipv0[is]
                ipi = ipv0[is, iv]
                ipf = ipia1[ipi]
                ip  = ipia0[ipi]
                xv = vertp0[ipf, 1] - vertp0[ip, 1]
                yv = vertp0[ipf, 2] - vertp0[ip, 2]
                zv = vertp0[ipf, 3] - vertp0[ip, 3]
                cis = xns0[is] * xv + yns0[is] * yv + zns0[is] * zv
                if cis == 0.0
                    betxe[is, iv] = 0.0
                    betye[is, iv] = 0.0
                    betze[is, iv] = 0.0
                else
                    betxe[is, iv] = -xv / cis
                    betye[is, iv] = -yv / cis
                    betze[is, iv] = -zv / cis
                end
                coef = xns0[is] * vertp0[ipf, 1] + yns0[is] * vertp0[ipf, 2] +
                       zns0[is] * vertp0[ipf, 3]
                x0[is, iv] = vertp0[ipf, 1] + betxe[is, iv] * coef
                y0[is, iv] = vertp0[ipf, 2] + betye[is, iv] * coef
                z0[is, iv] = vertp0[ipf, 3] + betze[is, iv] * coef
            end
        end

        # Contributions of the rest of faces
        for is in 1:nts00
            for iv in 1:nipv0[is]
                ip = ipv0[is, iv]
                if ia[ip] == 1
                    x0[is, iv] = vertp0[ip, 1]
                    y0[is, iv] = vertp0[ip, 2]
                    z0[is, iv] = vertp0[ip, 3]
                    betxe[is, iv] = 0.0
                    betye[is, iv] = 0.0
                    betze[is, iv] = 0.0
                else
                    # Find matching vertex in new faces
                    found = false
                    for is1 in (nts00+1):nts0
                        for iv1 in 1:nipv0[is1]
                            ip1 = ipv0[is1, iv1]
                            if ip1 == ip
                                betxe[is, iv] = betxe[is1, iv1]
                                betye[is, iv] = betye[is1, iv1]
                                betze[is, iv] = betze[is1, iv1]
                                x0[is, iv]    = x0[is1, iv1]
                                y0[is, iv]    = y0[is1, iv1]
                                z0[is, iv]    = z0[is1, iv1]
                                found = true
                                break
                            end
                        end
                        found && break
                    end
                end
            end
        end

        # Vectors K, L and M
        for is in (nts00+1):nts0
            iscut[is] = 1
        end
        fill!(sumk, 0.0)
        fill!(suml, 0.0)
        fill!(summ, 0.0)

        for is in 1:nts0
            ax = abs(yns0[is]); bx = abs(xns0[is]); cx = abs(zns0[is])
            if ax >= bx && ax >= cx
                iproj = 2; dnmax = yns0[is]
            elseif cx >= bx && cx >= ax
                iproj = 3; dnmax = zns0[is]
            else
                iproj = 1; dnmax = xns0[is]
            end

            ih = div(nipv0[is] - 2, 2)
            for i in 2:(ih + 1)
                ip  = 2 * i
                ip1 = ip - 1
                ip2 = ip - 2
                if iproj == 1
                    yv1 = y0[is, ip1] - y0[is, 1]; zv1 = z0[is, ip1] - z0[is, 1]
                    yv2 = y0[is, ip]  - y0[is, ip2]; zv2 = z0[is, ip] - z0[is, ip2]
                    ye1 = betye[is, ip1] - betye[is, 1]
                    ze1 = betze[is, ip1] - betze[is, 1]
                    ye2 = betye[is, ip]  - betye[is, ip2]
                    ze2 = betze[is, ip]  - betze[is, ip2]
                    sumk[is] += yv1 * zv2 - zv1 * yv2
                    suml[is] += yv1 * ze2 - zv1 * ye2 - (yv2 * ze1 - zv2 * ye1)
                    summ[is] += ye1 * ze2 - ze1 * ye2
                elseif iproj == 2
                    xv1 = x0[is, ip1] - x0[is, 1]; zv1 = z0[is, ip1] - z0[is, 1]
                    xv2 = x0[is, ip]  - x0[is, ip2]; zv2 = z0[is, ip] - z0[is, ip2]
                    xe1 = betxe[is, ip1] - betxe[is, 1]
                    ze1 = betze[is, ip1] - betze[is, 1]
                    xe2 = betxe[is, ip]  - betxe[is, ip2]
                    ze2 = betze[is, ip]  - betze[is, ip2]
                    sumk[is] += zv1 * xv2 - xv1 * zv2
                    suml[is] += zv1 * xe2 - xv1 * ze2 - (zv2 * xe1 - xv2 * ze1)
                    summ[is] += ze1 * xe2 - xe1 * ze2
                else
                    xv1 = x0[is, ip1] - x0[is, 1]; yv1 = y0[is, ip1] - y0[is, 1]
                    xv2 = x0[is, ip]  - x0[is, ip2]; yv2 = y0[is, ip] - y0[is, ip2]
                    xe1 = betxe[is, ip1] - betxe[is, 1]
                    ye1 = betye[is, ip1] - betye[is, 1]
                    xe2 = betxe[is, ip]  - betxe[is, ip2]
                    ye2 = betye[is, ip]  - betye[is, ip2]
                    sumk[is] += xv1 * yv2 - yv1 * xv2
                    suml[is] += xv1 * ye2 - yv1 * xe2 - (xv2 * ye1 - yv2 * xe1)
                    summ[is] += xe1 * ye2 - ye1 * xe2
                end
            end

            # Handle odd vertex
            if 2 * (ih + 1) < nipv0[is]
                nv = nipv0[is]
                if iproj == 1
                    yv1 = y0[is, nv] - y0[is, 1]; zv1 = z0[is, nv] - z0[is, 1]
                    yv2 = y0[is, 1] - y0[is, nv-1]; zv2 = z0[is, 1] - z0[is, nv-1]
                    ye1 = betye[is, nv] - betye[is, 1]
                    ze1 = betze[is, nv] - betze[is, 1]
                    ye2 = betye[is, 1] - betye[is, nv-1]
                    ze2 = betze[is, 1] - betze[is, nv-1]
                    sumk[is] += yv1 * zv2 - zv1 * yv2
                    suml[is] += yv1 * ze2 - zv1 * ye2 - (yv2 * ze1 - zv2 * ye1)
                    summ[is] += ye1 * ze2 - ze1 * ye2
                elseif iproj == 2
                    xv1 = x0[is, nv] - x0[is, 1]; zv1 = z0[is, nv] - z0[is, 1]
                    xv2 = x0[is, 1] - x0[is, nv-1]; zv2 = z0[is, 1] - z0[is, nv-1]
                    xe1 = betxe[is, nv] - betxe[is, 1]
                    ze1 = betze[is, nv] - betze[is, 1]
                    xe2 = betxe[is, 1] - betxe[is, nv-1]
                    ze2 = betze[is, 1] - betze[is, nv-1]
                    sumk[is] += zv1 * xv2 - xv1 * zv2
                    suml[is] += zv1 * xe2 - xv1 * ze2 - (zv2 * xe1 - xv2 * ze1)
                    summ[is] += ze1 * xe2 - xe1 * ze2
                else
                    xv1 = x0[is, nv] - x0[is, 1]; yv1 = y0[is, nv] - y0[is, 1]
                    xv2 = x0[is, 1] - x0[is, nv-1]; yv2 = y0[is, 1] - y0[is, nv-1]
                    xe1 = betxe[is, nv] - betxe[is, 1]
                    ye1 = betye[is, nv] - betye[is, 1]
                    xe2 = betxe[is, 1] - betxe[is, nv-1]
                    ye2 = betye[is, 1] - betye[is, nv-1]
                    sumk[is] += xv1 * yv2 - yv1 * xv2
                    suml[is] += xv1 * ye2 - yv1 * xe2 - (xv2 * ye1 - yv2 * xe1)
                    summ[is] += xe1 * ye2 - ye1 * xe2
                end
            end

            sumk[is] /= dnmax
            suml[is] /= dnmax
            summ[is] /= dnmax
        end

        # Coefficients of the analytical equation: c3*x^3 + c2*x^2 + c1*x + c0 = 0
        c3 = 0.0; c2_coef = 0.0; c1_coef = 0.0
        for is in (nts00+1):nts0
            c3 += summ[is]
            c2_coef += suml[is]
            c1_coef += sumk[is]
        end
        c0 = 6.0 * vaux
        for is in 1:nts00
            c2_coef += summ[is] * cs[is]
            c1_coef += suml[is] * cs[is]
            c0 += sumk[is] * cs[is]
        end

        vmaxl = vaux - (c3 * cmin_loc^3 + c2_coef * cmin_loc^2 + c1_coef * cmin_loc +
                  c0) / 6.0
        vminl = vaux - (c3 * cmax_loc^3 + c2_coef * cmax_loc^2 + c1_coef * cmax_loc +
                  c0) / 6.0

        if invert == 1
            vminll = vt - vmaxl
            vmaxl = vt - vminl
            vminl = vminll
        end

        sv = (vminl - v) * (vmaxl - v)
        if sv > 0.0 && (imax - imin) > 1 && imaxlold != imaxl
            if vmaxl > v
                vmax = vminl
                imax = iminl
            else
                vmin = vmaxl
                imin = imaxl
            end
            imaxlold = imaxl
            continue
        end

        c_sol = eqsol3d(c0, c1_coef, c2_coef, c3, cmin_loc, cmax_loc)
        errv = c3 * c_sol^3 + c2_coef * c_sol^2 + c1_coef * c_sol + c0
        if abs(errv) > tolc
            c_sol, _ = newton3d(c0, c1_coef, c2_coef, c3, cmin_loc, cmax_loc, c_sol)
        end

        if invert == 0
            c_ret = -c_sol
            cmax2 = -cmin_loc
            cmin2 = -cmax_loc
        else
            xnc_local = -xnc_local
            ync_local = -ync_local
            znc_local = -znc_local
            c_ret = c_sol
            cmax2 = cmax_loc
            cmin2 = cmin_loc
        end

        # Forced adjustment
        if c_ret < cmin2 || c_ret > cmax2
            cppol3d!(poly0, poly)
            icontn_min, icontp_min = inte3d!(poly0, cmin2, xnc_local, ync_local, znc_local)
            if icontp_min == 0
                volmin = 0.0
            else
                volmin = toolv3d(poly0)
            end
            cppol3d!(poly0, poly)
            icontn_max, icontp_max = inte3d!(poly0, cmax2, xnc_local, ync_local, znc_local)
            if icontp_max == 0
                volmax = 0.0
            else
                volmax = toolv3d(poly0)
            end
            if abs(volmin - v) < abs(volmax - v)
                c_ret = cmin2
            else
                c_ret = cmax2
            end
        end

        return c_ret
    end
end

"""
    enforv3d(poly::Polyhedron3D, v, vt, xnc, ync, znc) -> c

Non-mutating version: uses a task-local scratch copy, then solves.
"""
function enforv3d(poly::Polyhedron3D, v::Float64, vt::Float64,
                  xnc::Float64, ync::Float64, znc::Float64)
    work = _get_enforv3d_work(poly)::_Enforv3DWork
    cppol3d!(work.poly0, poly)
    return enforv3d!(work.poly0, v, vt, xnc, ync, znc)
end

# Reusable task-local workspace for initf3d hot path.
mutable struct _Initf3DWork
    icheck::Vector{Int}
    phiv::Vector{Float64}
    ia::Vector{Int}
    poly2::Polyhedron3D
    poly1::Polyhedron3D
    poly0::Polyhedron3D
    icheck_sub::Vector{Int}
    iscut_sub::Vector{Int}
    ipia0_sub::Vector{Int}
    ipia1_sub::Vector{Int}
    ntp_ref_sub::Base.RefValue{Int}
    nts_ref_sub::Base.RefValue{Int}
    ntv_ref_sub::Base.RefValue{Int}
    ipv1_sub::Matrix{Int}
    nipv1_sub::Vector{Int}
    nedge_sub::Vector{Int}
    ise_sub::Matrix{Int}
    ivise_sub::Matrix{Int}
    ipise_sub::Matrix{Int}
    ipmark_sub::Vector{Int}
end

function _Initf3DWork(poly::Polyhedron3D)
    ns_max = size(poly.ipv, 1)
    nv_max = size(poly.ipv, 2)
    nv_store = size(poly.vertp, 1)
    return _Initf3DWork(
        zeros(Int, nv_store),
        zeros(nv_store),
        zeros(Int, nv_store),
        cppol3d(poly),
        cppol3d(poly),
        cppol3d(poly),
        zeros(Int, nv_store),
        zeros(Int, ns_max),
        zeros(Int, nv_max),
        zeros(Int, nv_max),
        Ref(0),
        Ref(0),
        Ref(0),
        zeros(Int, ns_max, nv_max),
        zeros(Int, ns_max),
        zeros(Int, ns_max),
        zeros(Int, ns_max, nv_max),
        zeros(Int, ns_max, nv_max),
        zeros(Int, nv_max, 2),
        zeros(Int, nv_max),
    )
end

function _get_initf3d_work(poly::Polyhedron3D)
    key = (size(poly.ipv, 1), size(poly.ipv, 2), size(poly.vertp, 1))
    tls = task_local_storage()
    cache_any = get(tls, :_voftools_initf3d_work, nothing)
    cache = if cache_any === nothing
        c = Dict{NTuple{3, Int}, _Initf3DWork}()
        tls[:_voftools_initf3d_work] = c
        c
    else
        cache_any::Dict{NTuple{3, Int}, _Initf3DWork}
    end
    work = get(cache, key, nothing)
    if work === nothing
        work = _Initf3DWork(poly)
        cache[key] = work
    end
    return work::_Initf3DWork
end

# ============================= INITF3D =====================================
"""
    initf3d(func3d, poly::Polyhedron3D; nc::Int=10, tol::Float64=10.0) -> Float64

Initialize the material volume fraction in a polyhedral cell.

`func3d(x, y, z)` is a user-supplied implicit function defining the interface:
- `> 0` inside the material body
- `< 0` outside
- `= 0` on the interface

`nc` is the number of sub-cells along each coordinate axis.
`tol` is a tolerance for near-interface vertices.

Returns the volume fraction `vf ∈ [0, 1]`.
"""
function initf3d(func3d::F, poly::Polyhedron3D,
                 nc::Int, tol::Float64) where {F}
    ipv   = poly.ipv
    nipv  = poly.nipv
    vertp = poly.vertp
    xns   = poly.xns; yns = poly.yns; zns = poly.zns
    nts   = poly.nts; ntp = poly.ntp; ntv = poly.ntv
    ns_max = size(ipv, 1); nv_max = size(ipv, 2)
    work = _get_initf3d_work(poly)::_Initf3DWork

    # Coordinate extremes and vertex tagging
    xmin = 1.0e20; xmax = -1.0e20
    ymin = 1.0e20; ymax = -1.0e20
    zmin = 1.0e20; zmax = -1.0e20
    icontp = 0; icontn = 0

    icheck = work.icheck
    phiv = work.phiv
    ia = work.ia
    fill!(icheck, 0)

    for is in 1:nts
        for iv in 1:nipv[is]
            ip = ipv[is, iv]
            if icheck[ip] == 0
                icheck[ip] = 1
                xp = vertp[ip, 1]; yp = vertp[ip, 2]; zp = vertp[ip, 3]
                xmin = min(xmin, xp); xmax = max(xmax, xp)
                ymin = min(ymin, yp); ymax = max(ymax, yp)
                zmin = min(zmin, zp); zmax = max(zmax, zp)
                phiv[ip] = func3d(xp, yp, zp)
                if phiv[ip] >= 0.0
                    ia[ip] = 1; icontp += 1
                else
                    ia[ip] = 0; icontn += 1
                end
            end
        end
    end

    dx = xmax - xmin; dy = ymax - ymin; dz = zmax - zmin

    # Check if any vertex is near-interface
    iphi = 0
    phimin = 10.0 * max(dx, dy, dz)
    for i in 1:ntv
        phimin = min(phimin, abs(phiv[i]))
    end
    if phimin < tol * dx
        iphi = 1
    end

    if iphi == 0
        if icontp == ntv
            return 1.0
        end
        if icontn == ntv
            return 0.0
        end
    end

    # Total volume of the original polyhedron
    volt = toolv3d(poly)
    ddx = dx / nc; ddy = dy / nc; ddz = dz / nc

    vf = 0.0
    poly2 = work.poly2
    poly1 = work.poly1
    poly0 = work.poly0
    icheck_sub = work.icheck_sub

    iscut_sub = work.iscut_sub
    ipia0_sub = work.ipia0_sub
    ipia1_sub = work.ipia1_sub
    ntp_ref_sub = work.ntp_ref_sub
    nts_ref_sub = work.nts_ref_sub
    ntv_ref_sub = work.ntv_ref_sub

    ipv1_sub = work.ipv1_sub
    nipv1_sub = work.nipv1_sub
    nedge_sub = work.nedge_sub
    ise_sub = work.ise_sub
    ivise_sub = work.ivise_sub
    ipise_sub = work.ipise_sub
    ipmark_sub = work.ipmark_sub

    for ic in 1:nc
        xc = xmin + (ic - 1) * ddx
        # Copy polyhedron -> poly2
        cppol3d!(poly2, poly)
        # Truncate by x-planes
        if ic > 1
            cx1 = -xc
            inte3d!(poly2, cx1, 1.0, 0.0, 0.0)
        end
        cx2 = xc + ddx
        inte3d!(poly2, cx2, -1.0, 0.0, 0.0)

        for jc in 1:nc
            yc = ymin + (jc - 1) * ddy
            # Copy poly2 -> poly1
            cppol3d!(poly1, poly2)
            if jc > 1
                cy1 = -yc
                icontn_tmp, icontp_tmp = inte3d!(poly1, cy1, 0.0, 1.0, 0.0)
            end
            if (jc == 1) || icontp_tmp != 0
                cy2 = yc + ddy
                icontn_tmp2, icontp_tmp2 = inte3d!(poly1, cy2, 0.0, -1.0, 0.0)
                if icontp_tmp2 != 0
                    for kc in 1:nc
                        zc = zmin + (kc - 1) * ddz
                        # Copy poly1 -> poly0
                        cppol3d!(poly0, poly1)
                        if kc > 1
                            cz1 = -zc
                            icontn_z, icontp_z = inte3d!(poly0, cz1, 0.0, 0.0, 1.0)
                        end
                        if (kc == 1) || icontp_z != 0
                            cz2 = zc + ddz
                            icontn_z2, icontp_z2 = inte3d!(poly0, cz2, 0.0, 0.0, -1.0)
                            if icontp_z2 != 0
                                # Subcell determination by truncation
                                icontp_sub = 0; icontn_sub = 0
                                fill!(icheck_sub, 0)
                                for is in 1:poly0.nts
                                    for iv in 1:poly0.nipv[is]
                                        ip = poly0.ipv[is, iv]
                                        if icheck_sub[ip] == 0
                                            icheck_sub[ip] = 1
                                            x = poly0.vertp[ip, 1]
                                            y = poly0.vertp[ip, 2]
                                            z = poly0.vertp[ip, 3]
                                            phiv[ip] = func3d(x, y, z)
                                            if phiv[ip] >= 0.0
                                                ia[ip] = 1; icontp_sub += 1
                                            else
                                                ia[ip] = 0; icontn_sub += 1
                                            end
                                        end
                                    end
                                end

                                if icontn_sub == 0
                                    volf = toolv3d(poly0)
                                    vf += volf
                                elseif icontn_sub > 0 && icontp_sub > 0
                                    ntsini = poly0.nts
                                    # Truncate by the interface
                                    fill!(iscut_sub, 0)
                                    fill!(ipia0_sub, 0)
                                    fill!(ipia1_sub, 0)
                                    ntp_ref_sub[] = poly0.ntp
                                    nts_ref_sub[] = poly0.nts
                                    ntv_ref_sub[] = poly0.ntv
                                    _newpol3d_impl!(ia, ipia0_sub, ipia1_sub, poly0.ipv, iscut_sub,
                                                    poly0.nipv, ntp_ref_sub, nts_ref_sub, ntv_ref_sub,
                                                    1.0, poly0.xns, 0.0, poly0.yns, 0.0, poly0.zns,
                                                    ipv1_sub, nipv1_sub, nedge_sub, ise_sub,
                                                    ivise_sub, ipise_sub, ipmark_sub)
                                    poly0.ntp = ntp_ref_sub[]; poly0.nts = nts_ref_sub[]; poly0.ntv = ntv_ref_sub[]

                                    # Location of the new intersection points (multiple new faces)
                                    if poly0.nts > ntsini
                                        is2 = poly0.nts
                                        for is in (ntsini + 1):poly0.nts
                                            sumx = 0.0; sumy = 0.0; sumz = 0.0
                                            for iv in 1:poly0.nipv[is]
                                                ip  = poly0.ipv[is, iv]
                                                ip0 = ipia0_sub[ip]
                                                ip1 = ipia1_sub[ip]
                                                poly0.vertp[ip, 1] = poly0.vertp[ip0, 1] -
                                                    phiv[ip0] * (poly0.vertp[ip1, 1] - poly0.vertp[ip0, 1]) /
                                                    (phiv[ip1] - phiv[ip0])
                                                poly0.vertp[ip, 2] = poly0.vertp[ip0, 2] -
                                                    phiv[ip0] * (poly0.vertp[ip1, 2] - poly0.vertp[ip0, 2]) /
                                                    (phiv[ip1] - phiv[ip0])
                                                poly0.vertp[ip, 3] = poly0.vertp[ip0, 3] -
                                                    phiv[ip0] * (poly0.vertp[ip1, 3] - poly0.vertp[ip0, 3]) /
                                                    (phiv[ip1] - phiv[ip0])
                                                sumx += poly0.vertp[ip, 1]
                                                sumy += poly0.vertp[ip, 2]
                                                sumz += poly0.vertp[ip, 3]
                                            end
                                            # Add centroid as new vertex
                                            poly0.ntp += 1
                                            ntp0_new = poly0.ntp
                                            poly0.vertp[ntp0_new, 1] = sumx / poly0.nipv[is]
                                            poly0.vertp[ntp0_new, 2] = sumy / poly0.nipv[is]
                                            poly0.vertp[ntp0_new, 3] = sumz / poly0.nipv[is]
                                            # Replace the new face IS with NIPV0[IS] triangular faces
                                            for iv in 1:poly0.nipv[is]
                                                is2 += 1
                                                iv2 = iv + 1
                                                if iv2 > poly0.nipv[is]; iv2 = 1; end
                                                poly0.nipv[is2] = 3
                                                poly0.ipv[is2, 1] = ntp0_new
                                                poly0.ipv[is2, 2] = poly0.ipv[is, iv]
                                                poly0.ipv[is2, 3] = poly0.ipv[is, iv2]
                                                xv1 = poly0.vertp[poly0.ipv[is2, 2], 1] - poly0.vertp[poly0.ipv[is2, 1], 1]
                                                yv1 = poly0.vertp[poly0.ipv[is2, 2], 2] - poly0.vertp[poly0.ipv[is2, 1], 2]
                                                zv1 = poly0.vertp[poly0.ipv[is2, 2], 3] - poly0.vertp[poly0.ipv[is2, 1], 3]
                                                xv2 = poly0.vertp[poly0.ipv[is2, 3], 1] - poly0.vertp[poly0.ipv[is2, 2], 1]
                                                yv2 = poly0.vertp[poly0.ipv[is2, 3], 2] - poly0.vertp[poly0.ipv[is2, 2], 2]
                                                zv2 = poly0.vertp[poly0.ipv[is2, 3], 3] - poly0.vertp[poly0.ipv[is2, 2], 3]
                                                xm = yv1 * zv2 - zv1 * yv2
                                                ym = zv1 * xv2 - xv1 * zv2
                                                zm = xv1 * yv2 - yv1 * xv2
                                                amod = sqrt(xm^2 + ym^2 + zm^2)
                                                if amod != 0.0
                                                    poly0.xns[is2] = xm / amod
                                                    poly0.yns[is2] = ym / amod
                                                    poly0.zns[is2] = zm / amod
                                                else
                                                    poly0.nipv[is2] = 0
                                                end
                                            end
                                            # Cancel the original new face IS
                                            if is2 > is
                                                poly0.nipv[is] = 0
                                            end
                                        end
                                        poly0.nts = is2
                                    end
                                    volf = toolv3d(poly0)
                                    vf += volf
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    vf /= volt
    return vf
end

@inline function initf3d(func3d::F, poly::Polyhedron3D;
                         nc::Int=10, tol::Float64=10.0) where {F}
    return initf3d(func3d, poly, nc, tol)
end
