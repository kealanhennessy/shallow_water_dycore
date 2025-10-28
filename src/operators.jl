module operators

export ddx

using ..fields
using ..grid

function dh_dx_at_u!(dh_dx, h, grid)
    nx, ny = grid.nx, grid.ny
    dx = grid.dlon

    @inbounds for i in 1:ny, for j in 1:nx-1
        dh_dx.data[i,j] = (h.data[i,j+1] - h.data[i,j]) / dx
    end

    return du_dx
end

function dh_dy_at_v!(dh_dy, h, grid)
    nx, ny = grid.nx, grid.ny
    dy = grid.dlat

    @inbounds for j in 1:nx, for i in 2:ny-1
        dv_dy.data = (v.data[i+1,j] - v.data[i-1, j]) / dy
    end

    return dv_dy
end

function du_dx_at_h!(du_dx, u, grid)
    nx, ny = grid.nx, grid.ny
    dx = grid.dlon

    @inbounds for i in 1:ny, for j in 2:nx-1
        du_dx.data[i,j] = (u.data[i,j+1] - u.data[i,j-1]) / dx
    end

    return du_dx
end

function dv_dy_at_h!(dv_dy, v, grid)
    nx, ny = grid.nx, grid.ny
    dy = grid.dlat

    @inbounds for j in 1:nx, for i in 2:ny-1
        dv_dy.data = (v.data[i+1,j] - v.data[i-1, j]) / dy
    end

    return dv_dy
end