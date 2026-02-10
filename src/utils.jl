# ---------------------------------------------------------------------------
# Shared utility functions for VOFTools.jl
# ---------------------------------------------------------------------------

"""
    _signed_area_2d(ipv, ntv, vertp)

Compute twice the signed area of a 2D polygon using the optimised
even-odd pairing scheme from the original VOFTools code.
"""
function _signed_area_2d(ipv::AbstractVector{Int}, ntv::Int,
                         vertp::AbstractMatrix{Float64})
    sums = 0.0
    ih = div(ntv - 2, 2)
    for i in 2:(ih + 1)
        ip  = 2 * i
        ip1 = ip - 1
        ip2 = ip - 2
        xv1 = vertp[ipv[ip1], 1] - vertp[ipv[1], 1]
        yv1 = vertp[ipv[ip1], 2] - vertp[ipv[1], 2]
        xv2 = vertp[ipv[ip], 1]  - vertp[ipv[ip2], 1]
        yv2 = vertp[ipv[ip], 2]  - vertp[ipv[ip2], 2]
        sums += xv1 * yv2 - yv1 * xv2
    end
    if 2 * (ih + 1) < ntv
        xv1 = vertp[ipv[ntv], 1]   - vertp[ipv[1], 1]
        yv1 = vertp[ipv[ntv], 2]   - vertp[ipv[1], 2]
        xv2 = vertp[ipv[1], 1]     - vertp[ipv[ntv - 1], 1]
        yv2 = vertp[ipv[1], 2]     - vertp[ipv[ntv - 1], 2]
        sums += xv1 * yv2 - yv1 * xv2
    end
    return sums
end

"""
    _face_area_contribution(is, nipv_is, ipv, verti, xns, yns, zns)

Compute the area contribution of a single face `is` of a 3D polyhedron,
projected onto the largest normal component, using the same even-odd
pairing as the original Fortran.
"""
function _face_area_contribution(is::Int, nipv_is::Int,
                                 ipv::AbstractMatrix{Int},
                                 verti::AbstractMatrix{Float64},
                                 xns::AbstractVector{Float64},
                                 yns::AbstractVector{Float64},
                                 zns::AbstractVector{Float64})
    sump = 0.0
    # Choose projection direction
    ax = abs(xns[is]); ay = abs(yns[is]); az = abs(zns[is])
    if ay >= ax && ay >= az
        iproj = 2; dnmax = yns[is]
    elseif az >= ax && az >= ay
        iproj = 3; dnmax = zns[is]
    else
        iproj = 1; dnmax = xns[is]
    end

    ih = div(nipv_is - 2, 2)
    for i in 2:(ih + 1)
        ip  = 2 * i
        ip1 = ip - 1
        ip2 = ip - 2
        if iproj == 1
            yv1 = verti[ipv[is, ip1], 2] - verti[ipv[is, 1], 2]
            zv1 = verti[ipv[is, ip1], 3] - verti[ipv[is, 1], 3]
            yv2 = verti[ipv[is, ip],  2] - verti[ipv[is, ip2], 2]
            zv2 = verti[ipv[is, ip],  3] - verti[ipv[is, ip2], 3]
            sump += yv1 * zv2 - zv1 * yv2
        elseif iproj == 2
            xv1 = verti[ipv[is, ip1], 1] - verti[ipv[is, 1], 1]
            zv1 = verti[ipv[is, ip1], 3] - verti[ipv[is, 1], 3]
            xv2 = verti[ipv[is, ip],  1] - verti[ipv[is, ip2], 1]
            zv2 = verti[ipv[is, ip],  3] - verti[ipv[is, ip2], 3]
            sump += zv1 * xv2 - xv1 * zv2
        else
            xv1 = verti[ipv[is, ip1], 1] - verti[ipv[is, 1], 1]
            yv1 = verti[ipv[is, ip1], 2] - verti[ipv[is, 1], 2]
            xv2 = verti[ipv[is, ip],  1] - verti[ipv[is, ip2], 1]
            yv2 = verti[ipv[is, ip],  2] - verti[ipv[is, ip2], 2]
            sump += xv1 * yv2 - yv1 * xv2
        end
    end
    if 2 * (ih + 1) < nipv_is
        nv = nipv_is
        if iproj == 1
            yv1 = verti[ipv[is, nv], 2]   - verti[ipv[is, 1], 2]
            zv1 = verti[ipv[is, nv], 3]   - verti[ipv[is, 1], 3]
            yv2 = verti[ipv[is, 1],  2]   - verti[ipv[is, nv - 1], 2]
            zv2 = verti[ipv[is, 1],  3]   - verti[ipv[is, nv - 1], 3]
            sump += yv1 * zv2 - zv1 * yv2
        elseif iproj == 2
            xv1 = verti[ipv[is, nv], 1]   - verti[ipv[is, 1], 1]
            zv1 = verti[ipv[is, nv], 3]   - verti[ipv[is, 1], 3]
            xv2 = verti[ipv[is, 1],  1]   - verti[ipv[is, nv - 1], 1]
            zv2 = verti[ipv[is, 1],  3]   - verti[ipv[is, nv - 1], 3]
            sump += zv1 * xv2 - xv1 * zv2
        else
            xv1 = verti[ipv[is, nv], 1]   - verti[ipv[is, 1], 1]
            yv1 = verti[ipv[is, nv], 2]   - verti[ipv[is, 1], 2]
            xv2 = verti[ipv[is, 1],  1]   - verti[ipv[is, nv - 1], 1]
            yv2 = verti[ipv[is, 1],  2]   - verti[ipv[is, nv - 1], 2]
            sump += xv1 * yv2 - yv1 * xv2
        end
    end

    return sump, dnmax
end
