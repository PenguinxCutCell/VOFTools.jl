# ---------------------------------------------------------------------------
#  3D mesh generators for VOFTools.jl
#  Translated from mesh.f (Version 5, January 2020)
#  Copyright (C) 2016-2020 J. Lopez and J. Hernandez
# ---------------------------------------------------------------------------

"""
    cubicmesh() -> Polyhedron3D

Create a unit cube `[0,1]^3` polyhedron.
"""
function cubicmesh()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    xns[1] =  1.0; yns[1] =  0.0; zns[1] =  0.0
    xns[2] =  0.0; yns[2] = -1.0; zns[2] =  0.0
    xns[3] =  0.0; yns[3] =  0.0; zns[3] = -1.0
    xns[4] =  0.0; yns[4] =  1.0; zns[4] =  0.0
    xns[5] =  0.0; yns[5] =  0.0; zns[5] =  1.0
    xns[6] = -1.0; yns[6] =  0.0; zns[6] =  0.0

    nts = 6; ntv = 8; ntp = ntv
    for is in 1:6; nipv[is] = 4; end

    ipv[1,1]=1; ipv[1,2]=2; ipv[1,3]=3; ipv[1,4]=4
    ipv[2,1]=2; ipv[2,2]=1; ipv[2,3]=5; ipv[2,4]=6
    ipv[3,1]=3; ipv[3,2]=2; ipv[3,3]=6; ipv[3,4]=7
    ipv[4,1]=4; ipv[4,2]=3; ipv[4,3]=7; ipv[4,4]=8
    ipv[5,1]=1; ipv[5,2]=4; ipv[5,3]=8; ipv[5,4]=5
    ipv[6,1]=6; ipv[6,2]=5; ipv[6,3]=8; ipv[6,4]=7

    vertp[1,:] = [1.0, 0.0, 1.0]
    vertp[2,:] = [1.0, 0.0, 0.0]
    vertp[3,:] = [1.0, 1.0, 0.0]
    vertp[4,:] = [1.0, 1.0, 1.0]
    vertp[5,:] = [0.0, 0.0, 1.0]
    vertp[6,:] = [0.0, 0.0, 0.0]
    vertp[7,:] = [0.0, 1.0, 0.0]
    vertp[8,:] = [0.0, 1.0, 1.0]

    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    hexahemesh() -> Polyhedron3D

Create a general hexahedral cell by truncating a cube.
"""
function hexahemesh()
    poly = cubicmesh()

    # Cutting planes passing through (0,0,0)
    c = 0.0
    xnc = 1.0; ync = -0.1; znc = -0.1
    dm = sqrt(xnc^2 + ync^2 + znc^2)
    xnc /= dm; ync /= dm; znc /= dm
    inte3d!(poly, c, xnc, ync, znc)

    xnc = -0.05; ync = 1.0; znc = -0.1
    dm = sqrt(xnc^2 + ync^2 + znc^2)
    xnc /= dm; ync /= dm; znc /= dm
    inte3d!(poly, c, xnc, ync, znc)

    xnc = -0.05; ync = -0.1; znc = 1.0
    dm = sqrt(xnc^2 + ync^2 + znc^2)
    xnc /= dm; ync /= dm; znc /= dm
    inte3d!(poly, c, xnc, ync, znc)

    # Cutting planes passing through (1,1,1)
    xnc = -1.0; ync = 0.05; znc = 0.1
    dm = sqrt(xnc^2 + ync^2 + znc^2)
    xnc /= dm; ync /= dm; znc /= dm
    c = -(xnc + ync + znc)
    inte3d!(poly, c, xnc, ync, znc)

    xnc = 0.05; ync = -1.0; znc = 0.025
    dm = sqrt(xnc^2 + ync^2 + znc^2)
    xnc /= dm; ync /= dm; znc /= dm
    c = -(xnc + ync + znc)
    inte3d!(poly, c, xnc, ync, znc)

    xnc = 0.05; ync = 0.05; znc = -1.0
    dm = sqrt(xnc^2 + ync^2 + znc^2)
    xnc /= dm; ync /= dm; znc /= dm
    c = -(xnc + ync + znc)
    inte3d!(poly, c, xnc, ync, znc)

    restore3d!(poly)
    return poly
end

"""
    tetramesh() -> Polyhedron3D

Create a tetrahedral cell.
"""
function tetramesh()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    nts = 4; ntv = 4; ntp = ntv
    for is in 1:4; nipv[is] = 3; end

    ipv[1,1]=1; ipv[1,2]=3; ipv[1,3]=2
    ipv[2,1]=3; ipv[2,2]=1; ipv[2,3]=4
    ipv[3,1]=2; ipv[3,2]=3; ipv[3,3]=4
    ipv[4,1]=1; ipv[4,2]=2; ipv[4,3]=4

    vertp[1,:] = [0.0, 0.0, 0.0]
    vertp[2,:] = [0.91, 0.24, 1.0]
    vertp[3,:] = [0.72, 0.16, 0.07]
    vertp[4,:] = [1.0, 1.0, 1.0]

    for is in 1:nts
        ip1 = ipv[is,1]; ip2 = ipv[is,2]; ip3 = ipv[is,3]
        xv1 = vertp[ip2,1]-vertp[ip1,1]; yv1 = vertp[ip2,2]-vertp[ip1,2]; zv1 = vertp[ip2,3]-vertp[ip1,3]
        xv2 = vertp[ip3,1]-vertp[ip2,1]; yv2 = vertp[ip3,2]-vertp[ip2,2]; zv2 = vertp[ip3,3]-vertp[ip2,3]
        xn = yv1*zv2-zv1*yv2; yn = zv1*xv2-xv1*zv2; zn = xv1*yv2-yv1*xv2
        dm = sqrt(xn^2+yn^2+zn^2)
        xns[is]=xn/dm; yns[is]=yn/dm; zns[is]=zn/dm
    end

    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    dodecamesh() -> Polyhedron3D

Create a dodecahedral cell.
"""
function dodecamesh()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    a = 1.0 / sqrt(3.0)
    b = sqrt((3.0 - sqrt(5.0)) / 6.0)
    c = sqrt((3.0 + sqrt(5.0)) / 6.0)

    nts = 12; ntv = 20; ntp = ntv

    vertp[1,:]  = [ a,  a,  a];   vertp[2,:]  = [ a,  a, -a]
    vertp[3,:]  = [ a, -a,  a];   vertp[4,:]  = [ a, -a, -a]
    vertp[5,:]  = [-a,  a,  a];   vertp[6,:]  = [-a,  a, -a]
    vertp[7,:]  = [-a, -a,  a];   vertp[8,:]  = [-a, -a, -a]
    vertp[9,:]  = [ b,  c, 0.0];  vertp[10,:] = [-b,  c, 0.0]
    vertp[11,:] = [ b, -c, 0.0];  vertp[12,:] = [-b, -c, 0.0]
    vertp[13,:] = [ c, 0.0,  b];  vertp[14,:] = [ c, 0.0, -b]
    vertp[15,:] = [-c, 0.0,  b];  vertp[16,:] = [-c, 0.0, -b]
    vertp[17,:] = [0.0,  b,  c];  vertp[18,:] = [0.0, -b,  c]
    vertp[19,:] = [0.0,  b, -c];  vertp[20,:] = [0.0, -b, -c]

    for is in 1:12; nipv[is] = 5; end

    faces = [
        1 9 10 5 17;  1 17 18 3 13;  13 3 11 4 14;  10 6 16 15 5;
        4 20 19 2 14;  8 12 7 15 16;  1 13 14 2 9;   9 2 19 6 10;
        17 5 15 7 18;  7 12 11 3 18;  8 16 6 19 20;  8 20 4 11 12
    ]
    for is in 1:12, iv in 1:5
        ipv[is, iv] = faces[is, iv]
    end

    for is in 1:nts
        ip1 = ipv[is,1]; ip2 = ipv[is,2]; ip3 = ipv[is,3]
        xv1 = vertp[ip2,1]-vertp[ip1,1]; yv1 = vertp[ip2,2]-vertp[ip1,2]; zv1 = vertp[ip2,3]-vertp[ip1,3]
        xv2 = vertp[ip3,1]-vertp[ip2,1]; yv2 = vertp[ip3,2]-vertp[ip2,2]; zv2 = vertp[ip3,3]-vertp[ip2,3]
        xn = yv1*zv2-zv1*yv2; yn = zv1*xv2-xv1*zv2; zn = xv1*yv2-yv1*xv2
        dm = sqrt(xn^2+yn^2+zn^2)
        xns[is]=xn/dm; yns[is]=yn/dm; zns[is]=zn/dm
    end

    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    icosamesh() -> Polyhedron3D

Create an icosahedral cell.
"""
function icosamesh()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    t = (1.0 + sqrt(5.0)) / 2.0
    a = sqrt(1.0 + t^2)

    nts = 20; ntv = 12; ntp = ntv

    vertp[1,:]  = [ t/a,  1.0/a, 0.0];   vertp[2,:]  = [-t/a,  1.0/a, 0.0]
    vertp[3,:]  = [ t/a, -1.0/a, 0.0];   vertp[4,:]  = [-t/a, -1.0/a, 0.0]
    vertp[5,:]  = [ 1.0/a, 0.0,  t/a];   vertp[6,:]  = [ 1.0/a, 0.0, -t/a]
    vertp[7,:]  = [-1.0/a, 0.0,  t/a];   vertp[8,:]  = [-1.0/a, 0.0, -t/a]
    vertp[9,:]  = [0.0,  t/a,  1.0/a];   vertp[10,:] = [0.0, -t/a,  1.0/a]
    vertp[11,:] = [0.0,  t/a, -1.0/a];   vertp[12,:] = [0.0, -t/a, -1.0/a]

    for is in 1:20; nipv[is] = 3; end

    faces = [
        1 9 5;   1 6 11;  3 5 10;  3 12 6;
        2 7 9;   2 11 8;  4 10 7;  4 8 12;
        1 11 9;  2 9 11;  3 10 12; 4 12 10;
        5 3 1;   6 1 3;   7 2 4;   8 4 2;
        9 7 5;   10 5 7;  11 6 8;  12 8 6
    ]
    for is in 1:20, iv in 1:3
        ipv[is, iv] = faces[is, iv]
    end

    for is in 1:nts
        ip1 = ipv[is,1]; ip2 = ipv[is,2]; ip3 = ipv[is,3]
        xv1 = vertp[ip2,1]-vertp[ip1,1]; yv1 = vertp[ip2,2]-vertp[ip1,2]; zv1 = vertp[ip2,3]-vertp[ip1,3]
        xv2 = vertp[ip3,1]-vertp[ip2,1]; yv2 = vertp[ip3,2]-vertp[ip2,2]; zv2 = vertp[ip3,3]-vertp[ip2,3]
        xn = yv1*zv2-zv1*yv2; yn = zv1*xv2-xv1*zv2; zn = xv1*yv2-yv1*xv2
        dm = sqrt(xn^2+yn^2+zn^2)
        xns[is]=xn/dm; yns[is]=yn/dm; zns[is]=zn/dm
    end

    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    complexmesh() -> Polyhedron3D

Create the complex 18-face, 32-vertex polyhedron from the VOFTools test suite.
"""
function complexmesh()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    ntp = 32; ntv = ntp; nts = 18

    coords = [
        0.542827704131611521 0.103161810115955432 0.036812175979487702;
        0.600855902479279003 0.086471258131132003 0.056436950730949877;
        0.580277964506970667 0.163039123851804080 0.063539333285191152;
        0.595353203042033874 0.109843129311963592 0.046191353733568030;
        0.543760908872757964 0.094099097938159251 0.040511888762036471;
        0.576819361018633514 0.166401017225910497 0.070695345095658793;
        0.563208132688422847 0.155406657743345611 0.108888651199260125;
        0.539375361775393247 0.123540417024526422 0.108841211139604044;
        0.567496472062291590 0.128649916760601390 0.113035361390195654;
        0.600028673063590978 0.085143233872209040 0.059111439263878407;
        0.526213727798502173 0.075636612108421930 0.077638126137438007;
        0.585556321506447985 0.087406427365641859 0.092250691458415329;
        0.559954651117444135 0.133408124832590513 0.115142136143722096;
        0.524564746594361031 0.076122995267018601 0.083047795609934638;
        0.564747090248686634 0.145410923059303421 0.110587477743786994;
        0.514395820566891815 0.112406169868730282 0.074847147170632553;
        0.561075920034260656 0.154023426807431363 0.110265478973935432;
        0.600086961628965465 0.099753327716150253 0.051150970749093465;
        0.523529840650993616 0.075758256340824406 0.081652754048494577;
        0.580306952541070453 0.163300101039615730 0.064077552617180802;
        0.560784564570740995 0.154556861208633767 0.110001263915814051;
        0.562766593882476296 0.155209023571326599 0.109311173701962666;
        0.535792240519188834 0.158071087967216972 0.055846886029540778;
        0.562461680790533380 0.155100229495437003 0.109460804789442229;
        0.528789432867762255 0.128164972986322317 0.044046149723425451;
        0.535717736632420283 0.158416280041318303 0.056565341810804984;
        0.527966627100669772 0.148427477131785279 0.097455495587370988;
        0.531197945823301376 0.154922874892573614 0.054094280101616717;
        0.519703120392784435 0.149093981287893751 0.082950196610910618;
        0.530093766604715744 0.156612444958861563 0.058283874375187901;
        0.518229286434464531 0.144572125378080146 0.084271163208895244;
        0.520436917533184662 0.148367611160416940 0.087663831740394063
    ]
    for i in 1:32, j in 1:3
        vertp[i, j] = coords[i, j]
    end

    nipv_vals = [5,9,4,5,4,6,6,10,3,7,5,5,6,5,6,3,4,3]
    for is in 1:18; nipv[is] = nipv_vals[is]; end

    face_data = [
        (1, [5,1,4,18,2]),
        (2, [18,20,6,7,15,9,12,10,2]),
        (3, [10,11,5,2]),
        (4, [23,26,6,20,3]),
        (5, [20,18,4,3]),
        (6, [4,1,25,28,23,3]),
        (7, [11,19,16,25,1,5]),
        (8, [26,30,29,32,27,21,24,22,7,6]),
        (9, [22,15,7]),
        (10, [27,32,31,16,19,14,8]),
        (11, [14,12,9,13,8]),
        (12, [13,17,21,27,8]),
        (13, [15,22,24,17,13,9]),
        (14, [12,14,19,11,10]),
        (15, [31,29,30,28,25,16]),
        (16, [24,21,17]),
        (17, [28,30,26,23]),
        (18, [31,32,29])
    ]
    for (is, verts) in face_data
        for (iv, v) in enumerate(verts)
            ipv[is, iv] = v
        end
    end

    for is in 1:nts
        ip1 = ipv[is,1]; ip2 = ipv[is,2]; ip3 = ipv[is,3]
        xv1 = vertp[ip2,1]-vertp[ip1,1]; yv1 = vertp[ip2,2]-vertp[ip1,2]; zv1 = vertp[ip2,3]-vertp[ip1,3]
        xv2 = vertp[ip3,1]-vertp[ip2,1]; yv2 = vertp[ip3,2]-vertp[ip2,2]; zv2 = vertp[ip3,3]-vertp[ip2,3]
        xn = yv1*zv2-zv1*yv2; yn = zv1*xv2-xv1*zv2; zn = xv1*yv2-yv1*xv2
        dm = sqrt(xn^2+yn^2+zn^2)
        xns[is]=xn/dm; yns[is]=yn/dm; zns[is]=zn/dm
    end

    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

# ---------------------------------------------------------------------------
#  Helper: compute outward face normals from cross products
# ---------------------------------------------------------------------------
function _compute_normals!(xns, yns, zns, ipv, nipv, vertp, nts)
    for is in 1:nts
        ip1 = ipv[is,1]; ip2 = ipv[is,2]; ip3 = ipv[is,3]
        xv1 = vertp[ip2,1]-vertp[ip1,1]; yv1 = vertp[ip2,2]-vertp[ip1,2]; zv1 = vertp[ip2,3]-vertp[ip1,3]
        xv2 = vertp[ip3,1]-vertp[ip2,1]; yv2 = vertp[ip3,2]-vertp[ip2,2]; zv2 = vertp[ip3,3]-vertp[ip2,3]
        xn = yv1*zv2-zv1*yv2; yn = zv1*xv2-xv1*zv2; zn = xv1*yv2-yv1*xv2
        dm = sqrt(xn^2+yn^2+zn^2)
        xns[is]=xn/dm; yns[is]=yn/dm; zns[is]=zn/dm
    end
end

# ---------------------------------------------------------------------------
#  Helper: stellate a polyhedron (used by ncscubicmesh, ncdodecamesh, ncicosamesh)
# ---------------------------------------------------------------------------
function _stellate(base_poly::Polyhedron3D, a::Float64)
    ns = size(base_poly.ipv, 1); nv = size(base_poly.ipv, 2)
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    # Copy base data
    nts_base = base_poly.nts
    ntp0 = base_poly.ntp
    for ip in 1:ntp0, j in 1:3
        vertp[ip, j] = base_poly.vertp[ip, j]
    end

    nts0 = nts_base
    # For each base face: add apex and create triangular sub-faces
    for is in 1:nts_base
        xc = 0.0; yc = 0.0; zc = 0.0
        for iv in 1:base_poly.nipv[is]
            ip = base_poly.ipv[is, iv]
            xc += vertp[ip, 1]; yc += vertp[ip, 2]; zc += vertp[ip, 3]
        end
        ntp0 += 1
        vertp[ntp0, 1] = xc / base_poly.nipv[is] + a * base_poly.xns[is]
        vertp[ntp0, 2] = yc / base_poly.nipv[is] + a * base_poly.yns[is]
        vertp[ntp0, 3] = zc / base_poly.nipv[is] + a * base_poly.zns[is]
        for iv in 1:base_poly.nipv[is]
            iv2 = iv == base_poly.nipv[is] ? 1 : iv + 1
            nts0 += 1
            nipv[nts0] = 3
            ipv[nts0, 1] = base_poly.ipv[is, iv]
            ipv[nts0, 2] = base_poly.ipv[is, iv2]
            ipv[nts0, 3] = ntp0
        end
    end

    # Shift new faces to beginning (discard original faces)
    nts = nts0 - nts_base
    for is in 1:nts
        nipv[is] = nipv[nts_base + is]
        for iv in 1:nipv[is]
            ipv[is, iv] = ipv[nts_base + is, iv]
        end
    end
    ntp = ntp0; ntv = ntp

    _compute_normals!(xns, yns, zns, ipv, nipv, vertp, nts)

    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

# ---------------------------------------------------------------------------
#  Helper: extrude a 2D polygon profile into a 3D letter shape
#  Used by voftoolslogo
# ---------------------------------------------------------------------------
function _extrude_letter!(ipv, nipv, vertp, xns, yns, zns,
                          nts::Int, ntp::Int, front_verts::Vector{Tuple{Float64,Float64}},
                          depth::Float64)
    ntsi = nts + 1
    ntpi = ntp + 1

    # Front face vertices (at z = depth)
    for (x, y) in front_verts
        ntp += 1
        vertp[ntp, 1] = x; vertp[ntp, 2] = y; vertp[ntp, 3] = depth
    end

    # Front face
    nts += 1
    ntpl = ntp - ntpi + 1
    nipv[nts] = ntpl
    for iv in 1:ntpl
        ipv[nts, iv] = ntpi + iv - 1
    end

    # Back face vertices (at z = 0)
    for i in 1:ntpl
        ntp += 1
        vertp[ntp, 1] = vertp[ntpi + i - 1, 1]
        vertp[ntp, 2] = vertp[ntpi + i - 1, 2]
        vertp[ntp, 3] = 0.0
    end

    # Side wall quad faces
    for is2 in 2:(ntpl + 1)
        nts += 1
        nipv[nts] = 4
        if is2 == ntpl + 1
            ipv[nts, 1] = (ntpi - 1) + 1
        else
            ipv[nts, 1] = (ntpi - 1) + is2
        end
        ipv[nts, 2] = (ntpi - 1) + is2 - 1
        ipv[nts, 3] = (ntpi - 1) + is2 + ntpl - 1
        if is2 == ntpl + 1
            ipv[nts, 4] = (ntpi - 1) + 1 + ntpl
        else
            ipv[nts, 4] = (ntpi - 1) + is2 + ntpl
        end
    end

    # Back face (reversed winding)
    nts += 1
    nipv[nts] = nipv[ntsi]
    for iv in 1:nipv[nts]
        ipv[nts, iv] = ntpi + ntpl * 2 - iv
    end

    return nts, ntp
end

# ========================== New 3D Mesh Generators ===========================

"""
    distortcubicmesh() -> Polyhedron3D

Create a distorted unit cube with vertices 4 and 7 at y=0.25, giving 7 faces
(including two triangular faces).
"""
function distortcubicmesh()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    xns[1] =  1.0; yns[1] =  0.0; zns[1] =  0.0
    xns[2] =  0.0; yns[2] = -1.0; zns[2] =  0.0
    xns[3] =  0.0; yns[3] =  0.0; zns[3] = -1.0
    xns[5] =  0.0; yns[5] =  0.0; zns[5] =  1.0
    xns[6] = -1.0; yns[6] =  0.0; zns[6] =  0.0

    nts = 7; ntv = 8; ntp = ntv

    nipv[1]=4; ipv[1,1]=1; ipv[1,2]=2; ipv[1,3]=3; ipv[1,4]=4
    nipv[2]=4; ipv[2,1]=2; ipv[2,2]=1; ipv[2,3]=5; ipv[2,4]=6
    nipv[3]=4; ipv[3,1]=3; ipv[3,2]=2; ipv[3,3]=6; ipv[3,4]=7
    nipv[4]=3; ipv[4,1]=4; ipv[4,2]=3; ipv[4,3]=7
    nipv[5]=4; ipv[5,1]=1; ipv[5,2]=4; ipv[5,3]=8; ipv[5,4]=5
    nipv[6]=4; ipv[6,1]=6; ipv[6,2]=5; ipv[6,3]=8; ipv[6,4]=7
    nipv[7]=3; ipv[7,1]=4; ipv[7,2]=7; ipv[7,3]=8

    vertp[1,:] = [1.0, 0.0, 1.0]
    vertp[2,:] = [1.0, 0.0, 0.0]
    vertp[3,:] = [1.0, 1.0, 0.0]
    vertp[4,:] = [1.0, 0.25, 1.0]
    vertp[5,:] = [0.0, 0.0, 1.0]
    vertp[6,:] = [0.0, 0.0, 0.0]
    vertp[7,:] = [0.0, 0.25, 0.0]
    vertp[8,:] = [0.0, 1.0, 1.0]

    # Compute normals for triangular faces 4 and 7
    for is in (4, 7)
        ip1 = ipv[is,1]; ip2 = ipv[is,2]; ip3 = ipv[is,3]
        xv1 = vertp[ip2,1]-vertp[ip1,1]; yv1 = vertp[ip2,2]-vertp[ip1,2]; zv1 = vertp[ip2,3]-vertp[ip1,3]
        xv2 = vertp[ip3,1]-vertp[ip2,1]; yv2 = vertp[ip3,2]-vertp[ip2,2]; zv2 = vertp[ip3,3]-vertp[ip2,3]
        xn = yv1*zv2-zv1*yv2; yn = zv1*xv2-xv1*zv2; zn = xv1*yv2-yv1*xv2
        dm = sqrt(xn^2+yn^2+zn^2)
        xns[is]=xn/dm; yns[is]=yn/dm; zns[is]=zn/dm
    end

    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    ncpentapyramid() -> Polyhedron3D

Create a non-convex pentagonal pyramid with 6 faces and 6 vertices.
"""
function ncpentapyramid()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    nts = 6; ntv = 6; ntp = ntv

    vertp[1,:] = [0.04, 0.77, 0.0]
    vertp[2,:] = [0.0,  0.0,  0.0]
    vertp[3,:] = [0.49, 0.22, 0.0]
    vertp[4,:] = [1.0,  0.13, 0.0]
    vertp[5,:] = [0.16, 1.0,  0.0]
    vertp[6,:] = [0.1,  0.5,  1.0]

    nipv[1]=5; ipv[1,1]=1; ipv[1,2]=5; ipv[1,3]=4; ipv[1,4]=3; ipv[1,5]=2
    nipv[2]=3; ipv[2,1]=1; ipv[2,2]=2; ipv[2,3]=6
    nipv[3]=3; ipv[3,1]=2; ipv[3,2]=3; ipv[3,3]=6
    nipv[4]=3; ipv[4,1]=3; ipv[4,2]=4; ipv[4,3]=6
    nipv[5]=3; ipv[5,1]=4; ipv[5,2]=5; ipv[5,3]=6
    nipv[6]=3; ipv[6,1]=5; ipv[6,2]=1; ipv[6,3]=6

    _compute_normals!(xns, yns, zns, ipv, nipv, vertp, nts)
    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    nccubicpyramid() -> Polyhedron3D

Create a non-convex cell: a unit cube with a pyramid subtracted (vertex 9 at center).
9 faces, 9 vertices.
"""
function nccubicpyramid()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    nts = 9; ntv = 9; ntp = ntv

    nipv[1]=4; ipv[1,1]=1; ipv[1,2]=2; ipv[1,3]=3; ipv[1,4]=4
    nipv[2]=4; ipv[2,1]=2; ipv[2,2]=1; ipv[2,3]=5; ipv[2,4]=6
    nipv[3]=4; ipv[3,1]=3; ipv[3,2]=2; ipv[3,3]=6; ipv[3,4]=7
    nipv[4]=3; ipv[4,1]=4; ipv[4,2]=9; ipv[4,3]=8
    nipv[5]=4; ipv[5,1]=1; ipv[5,2]=4; ipv[5,3]=8; ipv[5,4]=5
    nipv[6]=4; ipv[6,1]=6; ipv[6,2]=5; ipv[6,3]=8; ipv[6,4]=7
    nipv[7]=3; ipv[7,1]=7; ipv[7,2]=8; ipv[7,3]=9
    nipv[8]=3; ipv[8,1]=3; ipv[8,2]=7; ipv[8,3]=9
    nipv[9]=3; ipv[9,1]=4; ipv[9,2]=3; ipv[9,3]=9

    vertp[1,:] = [1.0, 0.0, 1.0]
    vertp[2,:] = [1.0, 0.0, 0.0]
    vertp[3,:] = [1.0, 1.0, 0.0]
    vertp[4,:] = [1.0, 1.0, 1.0]
    vertp[5,:] = [0.0, 0.0, 1.0]
    vertp[6,:] = [0.0, 0.0, 0.0]
    vertp[7,:] = [0.0, 1.0, 0.0]
    vertp[8,:] = [0.0, 1.0, 1.0]
    vertp[9,:] = [0.5, 0.5, 0.5]

    _compute_normals!(xns, yns, zns, ipv, nipv, vertp, nts)
    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    ncscubicmesh() -> Polyhedron3D

Create a small stellated cubic mesh: each cube face is replaced by a pyramid
of triangular faces pushed outward by `a=1.0`. 24 triangular faces, 14 vertices.
"""
function ncscubicmesh()
    return _stellate(cubicmesh(), 1.0)
end

"""
    nchexahemesh() -> Polyhedron3D

Create a non-convex hexahedron (6 quad faces, 8 vertices) with vertices 1,2
moved inward to (0.5, 0.75, z).
"""
function nchexahemesh()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    nts = 6; ntv = 8; ntp = ntv

    nipv[1]=4; ipv[1,1]=1; ipv[1,2]=2; ipv[1,3]=3; ipv[1,4]=4
    nipv[2]=4; ipv[2,1]=2; ipv[2,2]=1; ipv[2,3]=5; ipv[2,4]=6
    nipv[3]=4; ipv[3,1]=3; ipv[3,2]=2; ipv[3,3]=6; ipv[3,4]=7
    nipv[4]=4; ipv[4,1]=4; ipv[4,2]=3; ipv[4,3]=7; ipv[4,4]=8
    nipv[5]=4; ipv[5,1]=1; ipv[5,2]=4; ipv[5,3]=8; ipv[5,4]=5
    nipv[6]=4; ipv[6,1]=6; ipv[6,2]=5; ipv[6,3]=8; ipv[6,4]=7

    vertp[1,:] = [0.5, 0.75, 1.0]
    vertp[2,:] = [0.5, 0.75, 0.0]
    vertp[3,:] = [1.0, 1.0,  0.0]
    vertp[4,:] = [1.0, 1.0,  1.0]
    vertp[5,:] = [0.0, 0.0,  1.0]
    vertp[6,:] = [0.0, 0.0,  0.0]
    vertp[7,:] = [0.0, 1.0,  0.0]
    vertp[8,:] = [0.0, 1.0,  1.0]

    _compute_normals!(xns, yns, zns, ipv, nipv, vertp, nts)
    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    ncdodecamesh() -> Polyhedron3D

Create a small stellated dodecahedron: each of the 12 pentagonal faces is
replaced by 5 triangular faces pushed outward by `a=0.5`.
60 triangular faces, 32 vertices.
"""
function ncdodecamesh()
    return _stellate(dodecamesh(), 0.5)
end

"""
    ncicosamesh() -> Polyhedron3D

Create a small stellated icosahedron: each of the 20 triangular faces is
replaced by 3 triangular faces pushed outward by `a=1.0`.
60 triangular faces, 32 vertices.
"""
function ncicosamesh()
    return _stellate(icosamesh(), 1.0)
end

"""
    nchollowedcube() -> Polyhedron3D

Create a unit cube with a half-length cubic hollow in its center.
12 faces (6 outer + 6 inner with reversed normals), 16 vertices.
"""
function nchollowedcube()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    # Outer cube normals
    xns[1] =  1.0; yns[1] =  0.0; zns[1] =  0.0
    xns[2] =  0.0; yns[2] = -1.0; zns[2] =  0.0
    xns[3] =  0.0; yns[3] =  0.0; zns[3] = -1.0
    xns[4] =  0.0; yns[4] =  1.0; zns[4] =  0.0
    xns[5] =  0.0; yns[5] =  0.0; zns[5] =  1.0
    xns[6] = -1.0; yns[6] =  0.0; zns[6] =  0.0

    nts = 12; ntv = 16; ntp = ntv

    # Outer faces (same as cubicmesh)
    nipv[1]=4; ipv[1,1]=1; ipv[1,2]=2; ipv[1,3]=3; ipv[1,4]=4
    nipv[2]=4; ipv[2,1]=2; ipv[2,2]=1; ipv[2,3]=5; ipv[2,4]=6
    nipv[3]=4; ipv[3,1]=3; ipv[3,2]=2; ipv[3,3]=6; ipv[3,4]=7
    nipv[4]=4; ipv[4,1]=4; ipv[4,2]=3; ipv[4,3]=7; ipv[4,4]=8
    nipv[5]=4; ipv[5,1]=1; ipv[5,2]=4; ipv[5,3]=8; ipv[5,4]=5
    nipv[6]=4; ipv[6,1]=6; ipv[6,2]=5; ipv[6,3]=8; ipv[6,4]=7

    # Inner faces: reversed normals and reversed winding, vertices offset by 8
    for is in 1:6
        xns[is+6] = -xns[is]; yns[is+6] = -yns[is]; zns[is+6] = -zns[is]
        nipv[is+6] = nipv[is]
        for iv in 1:4
            iv2 = 4 - iv + 1
            ipv[is+6, iv] = ipv[is, iv2] + 8
        end
    end

    # Outer vertices (same as cubicmesh)
    vertp[1,:] = [1.0, 0.0, 1.0]
    vertp[2,:] = [1.0, 0.0, 0.0]
    vertp[3,:] = [1.0, 1.0, 0.0]
    vertp[4,:] = [1.0, 1.0, 1.0]
    vertp[5,:] = [0.0, 0.0, 1.0]
    vertp[6,:] = [0.0, 0.0, 0.0]
    vertp[7,:] = [0.0, 1.0, 0.0]
    vertp[8,:] = [0.0, 1.0, 1.0]

    # Inner vertices (inset by 0.25 from each wall)
    d14 = 0.25
    vertp[9,:]  = [1.0-d14, 0.0+d14, 1.0-d14]
    vertp[10,:] = [1.0-d14, 0.0+d14, 0.0+d14]
    vertp[11,:] = [1.0-d14, 1.0-d14, 0.0+d14]
    vertp[12,:] = [1.0-d14, 1.0-d14, 1.0-d14]
    vertp[13,:] = [0.0+d14, 0.0+d14, 1.0-d14]
    vertp[14,:] = [0.0+d14, 0.0+d14, 0.0+d14]
    vertp[15,:] = [0.0+d14, 1.0-d14, 0.0+d14]
    vertp[16,:] = [0.0+d14, 1.0-d14, 1.0-d14]

    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    drilledcube() -> Polyhedron3D

Create a non-simply connected polyhedron: a unit cube with a rectangular tunnel
drilled through it. 12 faces (6 outer + 6 inner), 16 vertices.
"""
function drilledcube()
    poly = cubicmesh()
    ipv = poly.ipv; nipv = poly.nipv; vertp = poly.vertp
    xns = poly.xns; yns = poly.yns; zns = poly.zns
    nts = poly.nts; ntp = poly.ntp; ntv = poly.ntv

    # Inner tunnel normals
    xns[nts+1] = -1.0; yns[nts+1] =  0.0; zns[nts+1] =  0.0
    xns[nts+2] =  0.0; yns[nts+2] =  1.0; zns[nts+2] =  0.0
    xns[nts+3] =  0.0; yns[nts+3] =  0.0; zns[nts+3] =  1.0
    xns[nts+4] =  0.0; yns[nts+4] = -1.0; zns[nts+4] =  0.0
    xns[nts+5] =  0.0; yns[nts+5] =  0.0; zns[nts+5] = -1.0
    xns[nts+6] =  1.0; yns[nts+6] =  0.0; zns[nts+6] =  0.0

    nipv[nts+1]=4; ipv[nts+1,4]=ntp+1; ipv[nts+1,3]=ntp+2; ipv[nts+1,2]=ntp+3; ipv[nts+1,1]=ntp+4
    nipv[nts+2]=4; ipv[nts+2,4]=ntp+2; ipv[nts+2,3]=ntp+1; ipv[nts+2,2]=ntp+5; ipv[nts+2,1]=ntp+6
    nipv[nts+3]=4; ipv[nts+3,4]=ntp+3; ipv[nts+3,3]=ntp+2; ipv[nts+3,2]=ntp+6; ipv[nts+3,1]=ntp+7
    nipv[nts+4]=4; ipv[nts+4,4]=ntp+4; ipv[nts+4,3]=ntp+3; ipv[nts+4,2]=ntp+7; ipv[nts+4,1]=ntp+8
    nipv[nts+5]=4; ipv[nts+5,4]=ntp+1; ipv[nts+5,3]=ntp+4; ipv[nts+5,2]=ntp+8; ipv[nts+5,1]=ntp+5
    nipv[nts+6]=4; ipv[nts+6,4]=ntp+6; ipv[nts+6,3]=ntp+5; ipv[nts+6,2]=ntp+8; ipv[nts+6,1]=ntp+7

    d02 = 0.25; d12 = 0.75
    vertp[ntp+1,:] = [d12, 0.0, d12]
    vertp[ntp+2,:] = [d12, 0.0, d02]
    vertp[ntp+3,:] = [d12, 1.0, d02]
    vertp[ntp+4,:] = [d12, 1.0, d12]
    vertp[ntp+5,:] = [d02, 0.0, d12]
    vertp[ntp+6,:] = [d02, 0.0, d02]
    vertp[ntp+7,:] = [d02, 1.0, d02]
    vertp[ntp+8,:] = [d02, 1.0, d12]

    nts += 6; ntv += 8; ntp = ntv
    poly.nts = nts; poly.ntp = ntp; poly.ntv = ntv
    return poly
end

"""
    zigzagmesh(; nzigs::Int=5, doff::Float64=0.1) -> Polyhedron3D

Create a zig-zag shaped polyhedron. `nzigs` controls the number of zig-zag
sections, `doff` is the vertical offset.
"""
function zigzagmesh(; nzigs::Int=5, doff::Float64=0.1)
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    ntp = 4 * nzigs; ntv = ntp; nts = 2 * nzigs + 2

    for iv in 1:nzigs
        vertp[iv, 1] = Float64(iv - 1)
        vertp[iv, 2] = doff + mod(iv - 1, 2)
        vertp[iv, 3] = 0.0

        vertp[iv + nzigs, 1] = Float64(nzigs - (iv - 1) - 1)
        vertp[iv + nzigs, 2] = Float64(mod(nzigs - (iv - 1) - 1, 2))
        vertp[iv + nzigs, 3] = 0.0

        vertp[iv + 2*nzigs, 1] = Float64(iv - 1)
        vertp[iv + 2*nzigs, 2] = doff + mod(iv - 1, 2)
        vertp[iv + 2*nzigs, 3] = 1.0

        vertp[iv + 3*nzigs, 1] = Float64(nzigs - (iv - 1) - 1)
        vertp[iv + 3*nzigs, 2] = Float64(mod(nzigs - (iv - 1) - 1, 2))
        vertp[iv + 3*nzigs, 3] = 1.0
    end

    # Front face
    nipv[1] = 2 * nzigs
    nipv[2] = 2 * nzigs
    for iv in 1:(2*nzigs)
        ipv[1, iv] = iv
        ipv[2, iv] = 4 * nzigs - (iv - 1)
    end

    # End-cap face 3
    iv = 1
    nipv[3] = 4
    ipv[3, 1] = iv
    ipv[3, 2] = 2 * nzigs - (iv - 1)
    ipv[3, 3] = 4 * nzigs - (iv - 1)
    ipv[3, 4] = iv + 2 * nzigs

    # End-cap face 4
    iv = nzigs
    nipv[4] = 4
    ipv[4, 4] = iv
    ipv[4, 3] = 2 * nzigs - (iv - 1)
    ipv[4, 2] = 4 * nzigs - (iv - 1)
    ipv[4, 1] = iv + 2 * nzigs

    # Side wall quads
    for iv in 1:(nzigs - 1)
        nipv[iv + 4] = 4
        ipv[iv + 4, 1] = iv
        ipv[iv + 4, 2] = iv + 2 * nzigs
        ipv[iv + 4, 3] = iv + 2 * nzigs + 1
        ipv[iv + 4, 4] = iv + 1

        nipv[iv + 4 + (nzigs - 1)] = 4
        ipv[iv + 4 + (nzigs - 1), 1] = 2 * nzigs - (iv - 1)
        ipv[iv + 4 + (nzigs - 1), 2] = 2 * nzigs - (iv - 1) - 1
        ipv[iv + 4 + (nzigs - 1), 3] = 2 * nzigs - (iv - 1) - 1 + 2 * nzigs
        ipv[iv + 4 + (nzigs - 1), 4] = 2 * nzigs - (iv - 1) + 2 * nzigs
    end

    _compute_normals!(xns, yns, zns, ipv, nipv, vertp, nts)
    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end

"""
    voftoolslogo() -> Polyhedron3D

Create a complex extruded 3D "VOFTools" text logo.
"""
function voftoolslogo()
    ns = NS_DEFAULT; nv = NV_DEFAULT
    ipv  = zeros(Int, ns, nv)
    nipv = zeros(Int, ns)
    vertp = zeros(Float64, nv, 3)
    xns = zeros(Float64, ns); yns = zeros(Float64, ns); zns = zeros(Float64, ns)

    nts = 0; ntp = 0
    depth = 1.0

    # V letter
    factor = 1.0; disp = 0.0
    verts = [(disp+0.0*factor, 5.0*factor), (disp+2.0*factor, 0.0*factor),
             (disp+4.0*factor, 0.0*factor), (disp+6.0*factor, 5.0*factor),
             (disp+5.0*factor, 5.0*factor), (disp+3.0*factor, 1.0*factor),
             (disp+1.0*factor, 5.0*factor)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # O letter (outer)
    disp = 6.2; factor = 1.0
    verts = [(disp+0.0, 0.0), (disp+4.0, 0.0), (disp+4.0, 5.0),
             (disp+3.0, 4.0), (disp+3.0, 1.0), (disp+1.0, 1.0)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # O letter (inner)
    verts = [(disp+0.0, 0.0), (disp+1.0, 1.0), (disp+1.0, 4.0),
             (disp+3.0, 4.0), (disp+4.0, 5.0), (disp+0.0, 5.0)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # F letter
    disp += 4.6; factor = 1.0
    verts = [(disp+0.0, 0.0), (disp+1.0, 0.0), (disp+1.0, 2.0),
             (disp+3.0, 2.0), (disp+3.0, 3.0), (disp+1.0, 3.0),
             (disp+1.0, 4.0), (disp+3.4, 4.0), (disp+3.4, 5.0),
             (disp+0.0, 5.0)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # T letter
    disp += 3.45; factor = 1.0
    verts = [(disp+2.0, 0.0), (disp+3.0, 0.0), (disp+3.0, 4.0),
             (disp+4.65, 4.0), (disp+4.65, 5.0), (disp+0.35, 5.0),
             (disp+0.35, 4.0), (disp+2.0, 4.0)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # First 'o' letter (outer)
    disp += 4.2; factor = 0.7
    verts = [(disp+0.0*factor, 0.0*factor), (disp+4.0*factor, 0.0*factor),
             (disp+4.0*factor, 5.0*factor), (disp+3.0*factor, 4.0*factor),
             (disp+3.0*factor, 1.0*factor), (disp+1.0*factor, 1.0*factor)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # First 'o' letter (inner)
    verts = [(disp+0.0*factor, 0.0*factor), (disp+1.0*factor, 1.0*factor),
             (disp+1.0*factor, 4.0*factor), (disp+3.0*factor, 4.0*factor),
             (disp+4.0*factor, 5.0*factor), (disp+0.0*factor, 5.0*factor)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # Second 'o' letter (outer)
    disp += 3.3; factor = 0.7
    verts = [(disp+0.0*factor, 0.0*factor), (disp+4.0*factor, 0.0*factor),
             (disp+4.0*factor, 5.0*factor), (disp+3.0*factor, 4.0*factor),
             (disp+3.0*factor, 1.0*factor), (disp+1.0*factor, 1.0*factor)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # Second 'o' letter (inner)
    verts = [(disp+0.0*factor, 0.0*factor), (disp+1.0*factor, 1.0*factor),
             (disp+1.0*factor, 4.0*factor), (disp+3.0*factor, 4.0*factor),
             (disp+4.0*factor, 5.0*factor), (disp+0.0*factor, 5.0*factor)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # 'l' letter
    disp += 2.7; factor = 1.0
    verts = [(disp+1.0, 0.0), (disp+2.0, 0.0), (disp+2.0, 5.0),
             (disp+0.0, 5.0), (disp+0.0, 4.0), (disp+1.0, 4.0)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    # 's' letter
    disp += 2.5; factor = 0.7
    verts = [(disp+0.0*factor, 0.0*factor), (disp+4.0*factor, 0.0*factor),
             (disp+4.0*factor, 3.0*factor), (disp+1.0*factor, 3.0*factor),
             (disp+1.0*factor, 4.0*factor), (disp+4.0*factor, 4.0*factor),
             (disp+4.0*factor, 5.0*factor), (disp+0.0*factor, 5.0*factor),
             (disp+0.0*factor, 2.0*factor), (disp+3.0*factor, 2.0*factor),
             (disp+3.0*factor, 1.0*factor), (disp+0.0*factor, 1.0*factor)]
    nts, ntp = _extrude_letter!(ipv, nipv, vertp, xns, yns, zns, nts, ntp, verts, depth)

    ntv = ntp
    _compute_normals!(xns, yns, zns, ipv, nipv, vertp, nts)
    return Polyhedron3D(vertp, ipv, nipv, xns, yns, zns, nts, ntp, ntv)
end
