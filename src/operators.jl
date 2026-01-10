module operators

using ..fields
using ..grid

function d_ucosφ_dλ!(out::Field, u::Field; g::Grid)
    nx, ny = g.nx, g.ny
    dλ = g.dλ
    cosφ = g.η.cosφ   # latitudes aligned with u in your grid choice

    @inbounds for j in 1:ny
        cj = cosφ[j]
        for i in 1:nx
            im1 = (i == 1) ? nx : i - 1   # periodic in longitude
            out.data[i,j] = (u.data[i,j]*cj - u.data[im1,j]*cj) / dλ
        end
    end
    return out
end

function ∂φ!(df::Field, f::Field; g::Grid)

end

function div!(df::Field, f::Field; g::Grid)

end

function grad!(df::Field, f::Field; g::Grid)

end

end