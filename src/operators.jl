module operators

using ..fields
using ..grid

function div_horizontal_velocity!(∇⋅u::Field, u::Field, v::Field; grid::Grid)
    div, u, v = ∇⋅u.data, u.data, v.data

    nx, ny = grid.nx, grid.ny
    
    @inbounds for j ∈ 2:ny+1
        inv_a_cosφ = grid.η.inv_a_cosφ[j]
        cosφ = grid.v.cosφ[j]
        cosφm1 = grid.v.cosφ[j-1]
        for i ∈ 2:nx+1
            u_λ = (u[i,j] - u[i-1,j]) / grid.dλ
            v_φ = (v[i,j]*cosφ - v[i,j-1]*cosφm1) / grid.dφ
            div[i,j] = inv_a_cosφ * u_λ + inv_a_cosφ * v_φ
        end
    end
    return nothing
end

function zonal_mass!(∂η_∂λ::Field, η::Field; grid::Grid)
    partial, η = ∂η_∂λ.data, η.data

    nx, ny = grid.nx, grid.ny

    @inbounds for j ∈ 2:ny+1
        for i ∈ 2:nx+1
            ∂η_∂λ[i,j] = (η[i,j] - η[i-1,j]) / grid.dλ
        end
    end
end

function meridional_mass!(∂η_∂φ::Field, η::Field; grid::Grid)
    partial, η = ∂η_∂φ.data, η.data

    nx, ny = grid.nx, grid.ny

    @inbounds for j ∈ 2:ny+1
        for i ∈ 2:nx+1
            ∂η_∂λ[i,j] = (η[i,j] - η[i,j-1]) / grid.dφ
        end
    end
end

end