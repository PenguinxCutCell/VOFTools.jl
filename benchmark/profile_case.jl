#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using VOFTools
using Profile

function build_case(case::Symbol)
    if case == :toolv3d
        poly = cubicmesh()
        return () -> toolv3d(poly)
    elseif case == :toolv2d
        poly = squaremesh()
        return () -> toolv2d(poly)
    elseif case == :enforv3d
        poly = cubicmesh()
        vt = toolv3d(poly)
        v = 0.5 * vt
        xnc = 0.57735027
        ync = 0.57735027
        znc = 0.57735027
        return () -> enforv3d(poly, v, vt, xnc, ync, znc)
    elseif case == :enforv2d
        poly = squaremesh()
        vt = toolv2d(poly)
        v = 0.5 * vt
        xnc = 0.57735027
        ync = 0.57735027
        return () -> enforv2d(poly, v, vt, xnc, ync)
    elseif case == :initf3d
        poly = cubicmesh()
        func3d(x, y, z) = -((x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2 - 0.6^2)
        return () -> initf3d(func3d, poly; nc=10, tol=10.0)
    elseif case == :initf2d
        poly = squaremesh()
        func2d(x, y) = -((x - 0.5)^2 + (y - 0.5)^2 - 0.6^2)
        return () -> initf2d(func2d, poly; nc=10, tol=10.0)
    else
        error("Unknown case: $(case). Use one of :toolv3d, :toolv2d, :enforv3d, :enforv2d, :initf3d, :initf2d.")
    end
end

function benchmark_profile(; case::Symbol=:enforv3d, n_eval::Int=5000, profile_n::Int=1000)
    f = build_case(case)

    # Warm-up compilation before measuring runtime/allocations.
    f()

    elapsed = @elapsed begin
        for _ in 1:n_eval
            f()
        end
    end

    total_alloc_bytes = @allocated begin
        for _ in 1:n_eval
            f()
        end
    end

    println("Case: $(case)")
    println("Runs: $(n_eval)")
    println("Total time (s): ", elapsed)
    println("Avg time (us/run): ", elapsed * 1e6 / n_eval)
    println("Total allocated (bytes): ", total_alloc_bytes)
    println("Avg allocated (bytes/run): ", total_alloc_bytes / n_eval)

    Profile.clear()
    @profile begin
        for _ in 1:profile_n
            f()
        end
    end

    println("\nCPU profile (flat, by sample count):")
    Profile.print(format=:flat, sortedby=:count, mincount=5)

    return nothing
end

case = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :enforv3d
n_eval = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5000
profile_n = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1000

benchmark_profile(; case=case, n_eval=n_eval, profile_n=profile_n)
