using Test
using VOFTools

@testset "VOFTools.jl" begin

# ====================================================================
#  2D TESTS (matching test2d.f)
# ====================================================================
@testset "2D: Square mesh - toolv2d" begin
    poly = squaremesh()
    vt = toolv2d(poly)
    @test vt ≈ 1.0 atol=1e-12
end

@testset "2D: enforv2d on square" begin
    poly = squaremesh()
    vt = toolv2d(poly)
    f = 0.5
    v = f * vt
    xnc = 0.70710678
    ync = 0.70710678
    c = enforv2d(poly, v, vt, xnc, ync)
    # The solution should be a reasonable constant
    @test isfinite(c)
    # Verify: truncate a copy and check volume
    poly0 = cppol2d(poly)
    icontn, icontp = inte2d!(poly0, c, xnc, ync)
    vtrunc = toolv2d(poly0)
    @test vtrunc ≈ v atol=1e-6
end

@testset "2D: enforv2dsz on square" begin
    poly = squaremesh()
    vt = toolv2d(poly)
    f = 0.5
    v = f * vt
    xnc = 0.70710678
    ync = 0.70710678
    dx = abs(poly.vertp[1,1] - poly.vertp[2,1])
    dy = abs(poly.vertp[2,2] - poly.vertp[3,2])
    c = enforv2dsz(dx, dy, v, poly.vertp, xnc, ync)
    @test isfinite(c)
    # Verify by truncation
    poly0 = cppol2d(poly)
    inte2d!(poly0, c, xnc, ync)
    vtrunc = toolv2d(poly0)
    @test vtrunc ≈ v atol=1e-6
end

@testset "2D: dist2d" begin
    poly = squaremesh()
    vt = toolv2d(poly)
    f = 0.5
    v = f * vt
    xnc = 0.70710678; ync = 0.70710678
    c = enforv2d(poly, v, vt, xnc, ync)
    poly0 = cppol2d(poly)
    inte2d!(poly0, c, xnc, ync)
    ntp0 = poly0.ntp
    x_seg = [poly0.vertp[ntp0-1, 1], poly0.vertp[ntp0, 1]]
    y_seg = [poly0.vertp[ntp0-1, 2], poly0.vertp[ntp0, 2]]
    xp = 0.0; yp = 0.0
    d = dist2d(x_seg, y_seg, xp, yp)
    @test d >= 0.0
    @test isfinite(d)
end

@testset "2D: Other mesh types" begin
    for (name, meshfn) in [("hexago", hexagomesh), ("triangle", trianglemesh),
                            ("quadrangle", quadranglemesh), ("pentagon", pentagonmesh),
                            ("hexagon", hexagonmesh)]
        poly = meshfn()
        vt = toolv2d(poly)
        @test vt > 0.0
        f = 0.5; v = f * vt
        xnc = 0.70710678; ync = 0.70710678
        c = enforv2d(poly, v, vt, xnc, ync)
        @test isfinite(c)
    end
end

# ====================================================================
#  3D TESTS (matching test3d.f)
# ====================================================================
@testset "3D: Cubic mesh - toolv3d" begin
    poly = cubicmesh()
    vt = toolv3d(poly)
    @test vt ≈ 1.0 atol=1e-12
end

@testset "3D: enforv3d on cube" begin
    poly = cubicmesh()
    vt = toolv3d(poly)
    f = 0.5
    v = f * vt
    xnc = 0.57735027; ync = 0.57735027; znc = 0.57735027
    c = enforv3d(poly, v, vt, xnc, ync, znc)
    @test isfinite(c)
    # Verify by truncation
    poly0 = cppol3d(poly)
    inte3d!(poly0, c, xnc, ync, znc)
    vtrunc = toolv3d(poly0)
    @test vtrunc ≈ v atol=1e-6
end

@testset "3D: enforv3dsz on cube" begin
    poly = cubicmesh()
    vt = toolv3d(poly)
    f = 0.5; v = f * vt
    xnc = 0.57735027; ync = 0.57735027; znc = 0.57735027
    dx = abs(poly.vertp[1,1] - poly.vertp[5,1])
    dy = abs(poly.vertp[1,2] - poly.vertp[4,2])
    dz = abs(poly.vertp[2,3] - poly.vertp[1,3])
    c = enforv3dsz(dx, dy, dz, v, poly.vertp, xnc, ync, znc)
    @test isfinite(c)
    poly0 = cppol3d(poly)
    inte3d!(poly0, c, xnc, ync, znc)
    vtrunc = toolv3d(poly0)
    @test vtrunc ≈ v atol=1e-6
end

@testset "3D: dist3d" begin
    poly = cubicmesh()
    vt = toolv3d(poly)
    f = 0.5; v = f * vt
    xnc = 0.57735027; ync = 0.57735027; znc = 0.57735027
    c = enforv3d(poly, v, vt, xnc, ync, znc)
    poly0 = cppol3d(poly)
    inte3d!(poly0, c, xnc, ync, znc)
    n = poly0.nipv[poly0.nts]
    xv = [poly0.vertp[poly0.ipv[poly0.nts, i], 1] for i in 1:n]
    yv = [poly0.vertp[poly0.ipv[poly0.nts, i], 2] for i in 1:n]
    zv = [poly0.vertp[poly0.ipv[poly0.nts, i], 3] for i in 1:n]
    xp = 0.0; yp = 0.0; zp = 0.0
    d = dist3d(xv, yv, zv, xp, yp, zp)
    @test d >= 0.0
    @test isfinite(d)
end

@testset "3D: Tetrahedron mesh" begin
    poly = tetramesh()
    vt = toolv3d(poly)
    @test vt > 0.0
    f = 0.5; v = f * vt
    xnc = 0.57735027; ync = 0.57735027; znc = 0.57735027
    c = enforv3d(poly, v, vt, xnc, ync, znc)
    @test isfinite(c)
end

@testset "3D: Other mesh types" begin
    for (name, meshfn) in [("dodeca", dodecamesh), ("icosa", icosamesh)]
        poly = meshfn()
        vt = toolv3d(poly)
        @test vt > 0.0
        f = 0.5; v = f * vt
        xnc = 0.57735027; ync = 0.57735027; znc = 0.57735027
        c = enforv3d(poly, v, vt, xnc, ync, znc)
        @test isfinite(c)
    end
end

@testset "3D: eqsol3d and newton3d" begin
    # Test a simple case: x^3 - 1 = 0 → x = 1
    csol = eqsol3d(-1.0, 0.0, 0.0, 1.0, 0.0, 2.0)
    @test csol ≈ 1.0 atol=1e-8
    # Newton refinement
    csol2, isol = newton3d(-1.0, 0.0, 0.0, 1.0, 0.0, 2.0, 0.5)
    @test csol2 ≈ 1.0 atol=1e-12
end

# ====================================================================
#  INITF TESTS (new in v3.2)
# ====================================================================
@testset "2D: initf2d circle on square" begin
    # Circle with radius 0.325 centered at (0.5, 0.5)
    # Exact area = π * 0.325^2 ≈ 0.331831
    func2d1(x, y) = -((x - 0.5)^2 + (y - 0.5)^2 - 0.325^2)
    poly = squaremesh()
    vf = initf2d(func2d1, poly; nc=10, tol=10.0)
    exact_frac = π * 0.325^2 / 1.0  # area fraction in unit square
    @test vf ≈ exact_frac atol=0.01
    @test 0.0 < vf < 1.0
end

@testset "2D: initf2d ellipse on square" begin
    # Ellipse: semi-major 0.6, semi-minor 0.2, center (0.5,0.5)
    func2d2(x, y) = 1.0 - ((x - 0.5) / 0.6)^2 - ((y - 0.5) / 0.2)^2
    poly = squaremesh()
    vf = initf2d(func2d2, poly; nc=10, tol=10.0)
    # Fortran reference (nc=10): 0.33703704811172669
    @test vf ≈ 0.33703704811172669 rtol=1e-6
    @test 0.0 < vf < 1.0
end

@testset "2D: initf2d fully inside" begin
    # Function that is positive everywhere on a unit square
    func_all_in(x, y) = 1.0
    poly = squaremesh()
    vf = initf2d(func_all_in, poly; nc=5, tol=10.0)
    @test vf ≈ 1.0 atol=1e-10
end

@testset "2D: initf2d fully outside" begin
    func_all_out(x, y) = -1.0
    poly = squaremesh()
    vf = initf2d(func_all_out, poly; nc=5, tol=10.0)
    @test vf ≈ 0.0 atol=1e-10
end

@testset "3D: initf3d sphere on cube" begin
    # Sphere with radius 0.6 centered at (0.5, 0.5, 0.5)
    func3d1(x, y, z) = -((x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2 - 0.6^2)
    poly = cubicmesh()
    vf = initf3d(func3d1, poly; nc=10, tol=10.0)
    # Fortran reference (nc=10): 0.78833321559055702
    @test vf ≈ 0.78833321559055702 rtol=1e-6
    @test 0.0 < vf < 1.0
end

@testset "3D: initf3d fully inside" begin
    func_all_in(x, y, z) = 1.0
    poly = cubicmesh()
    vf = initf3d(func_all_in, poly; nc=5, tol=10.0)
    @test vf ≈ 1.0 atol=1e-10
end

@testset "3D: initf3d fully outside" begin
    func_all_out(x, y, z) = -1.0
    poly = cubicmesh()
    vf = initf3d(func_all_out, poly; nc=5, tol=10.0)
    @test vf ≈ 0.0 atol=1e-10
end

end  # @testset "VOFTools.jl"
