module operators

using ..fields
using ..grid

function div_horizontal_velocity!(div::Field, u::Field, v::Field; grid::Grid)
    div, u, v = div.data, u.data, v.data
    
    nx, ny = grid.nx, grid.ny
    
    @inbounds for j ∈ 2:ny+1
        inv_a_cosφ = grid.η.inv_a_cosφ(j)
        cosφ = grid.v.cosφ(j)
        for i ∈ 2:nx+1
            u_λ = (u[i,j] - u[i-1,j]) / grid.dλ
            v_φ = (v[i,j]*cosφ - v[i,j-1]*cosφ) / grid.dφ
            div[i,j] = inv_a_cosφ * u_λ + inv_a_cosφ * v_φ
        end
    end
    return nothing
end

function ∂φ!(df::Field, f::Field; g::Grid)

end

function div!(df::Field, f::Field; g::Grid)

end

function grad!(df::Field, f::Field; g::Grid)

end

end