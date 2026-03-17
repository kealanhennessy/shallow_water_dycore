module operators
 
using ..fields
using ..grid
 
export ∂η_∂λ!, ∂η_∂φ!, ∂u_∂λ!, ∂vcosφ_∂φ!, interp_v_to_u!, interp_u_to_v!
 
# ------------------------------------------------------------------
# Index conventions (HALO = 1)
#
# η, u interior : i ∈ 2:nx+1,   j ∈ 2:ny+1
# v   interior  : i ∈ 2:nx+1,   j ∈ 2:ny+2
#
# C-grid layout (η at cell centre, u east face, v north face):
#
#        v(i,j+1)
#           |
#   u(i-1,j)-η(i,j)-u(i,j)
#           |
#         v(i,j)
#
# Zonal periodicity is handled via halo cells: the west halo
# (i=1) mirrors i=nx+1 and the east halo (i=nx+2) mirrors i=2.
# Halo filling must be called before operators are evaluated.
#
# Polar boundary: v = 0 at j=2 (south pole) and j=ny+2 (north
# pole) in the v-array. These are the first and last interior
# v-points and must be zeroed after each RHS evaluation.
# ------------------------------------------------------------------
 
# ------------------------------------------------------------------
# ∂η/∂λ  :  η-points → u-points
#
# Computes the zonal derivative of η at each u-point using a
# centred difference across the cell face:
#
#   (∂η/∂λ)_{i,j} = (η[i,j] - η[i-1,j]) / dλ
#
# The full zonal pressure gradient in the u-equation is:
#   -(g / a cosφ) * ∂η/∂λ
# The metric factor (g / a cosφ) is applied in rhs.jl.
# ------------------------------------------------------------------
"""
Zonal derivative of η at u-points. Output sized nx × ny.
"""
function ∂η_∂λ!(out::Field, η::Field, g::Grid)
    o = out.data
    h = η.data
    nx, ny = g.nx, g.ny
 
    @inbounds for j ∈ 2:ny+1
        for i ∈ 2:nx+1
            o[i,j] = (h[i,j] - h[i-1,j]) / g.dλ
        end
    end
    return nothing
end
 
# ------------------------------------------------------------------
# ∂η/∂φ  :  η-points → v-points
#
# Computes the meridional derivative of η at each v-point:
#
#   (∂η/∂φ)_{i,j} = (η[i,j] - η[i,j-1]) / dφ
#
# The full meridional pressure gradient in the v-equation is:
#   -(g / a) * ∂η/∂φ
# The metric factor (g / a) is applied in rhs.jl.
#
# Note: the loop runs j ∈ 2:ny+1, producing output at interior
# v-points only. The polar v-points (j=2 south, j=ny+2 north)
# are set to zero by the boundary condition, not computed here.
# ------------------------------------------------------------------
"""
Meridional derivative of η at v-points. Output sized nx × ny.
"""
function ∂η_∂φ!(out::Field, η::Field, g::Grid)
    o = out.data
    h = η.data
    nx, ny = g.nx, g.ny
 
    @inbounds for j ∈ 2:ny+1
        for i ∈ 2:nx+1
            o[i,j] = (h[i,j] - h[i,j-1]) / g.dφ
        end
    end
    return nothing
end
 
# ------------------------------------------------------------------
# ∂u/∂λ  :  u-points → η-points
#
# Computes the zonal derivative of u at each η-point. Because u is
# staggered half a cell east of η, the η-point at (i,j) is bounded
# by u-points at (i,j) to the east and (i-1+1,j) = (i,j) ... so
# the correct forward difference is:
#
#   (∂u/∂λ)_{i,j} = (u[i+1,j] - u[i,j]) / dλ
#
# The full zonal flux divergence contribution in the η-equation is:
#   -(H / a cosφ) * ∂u/∂λ
# The metric factor is applied in rhs.jl.
# ------------------------------------------------------------------
"""
Zonal derivative of u at η-points. Output sized nx × ny.
"""
function ∂u_∂λ!(out::Field, u::Field, g::Grid)
    o = out.data
    ux = u.data
    nx, ny = g.nx, g.ny
 
    @inbounds for j ∈ 2:ny+1
        for i ∈ 2:nx+1
            o[i,j] = (ux[i+1,j] - ux[i,j]) / g.dλ
        end
    end
    return nothing
end
 
# ------------------------------------------------------------------
# ∂(v cosφ)/∂φ  :  v-points → η-points
#
# Computes the meridional flux divergence at each η-point. The
# cosφ weighting is applied at the v-point locations (g.v.cosφ),
# which are staggered half a cell north of η-points. The η-point
# at (i,j) is bounded by v-points at (i,j) to the south and
# (i,j+1) to the north:
#
#   (∂(v cosφ)/∂φ)_{i,j} =
#       (v[i,j+1]*cosφ[j+1] - v[i,j]*cosφ[j]) / dφ
#
# The full meridional flux divergence in the η-equation is:
#   -(H / a cosφ) * ∂(v cosφ)/∂φ
# The metric factor is applied in rhs.jl.
# ------------------------------------------------------------------
"""
Meridional derivative of v cosφ at η-points. Output sized nx × ny.
"""
function ∂vcosφ_∂φ!(out::Field, v::Field, g::Grid)
    o   = out.data
    vx  = v.data
    nx, ny = g.nx, g.ny
 
    @inbounds for j ∈ 2:ny+1
        # these are indexed differently because they lack halos
        cosφ_n = g.v.cosφ[j]       # cosφ at northern v-face (j+1 in v-array)
        cosφ_s = g.v.cosφ[j-1]     # cosφ at southern v-face (j   in v-array)
        for i ∈ 2:nx+1
            o[i,j] = (vx[i,j+1]*cosφ_n - vx[i,j]*cosφ_s) / g.dφ
        end
    end
    return nothing
end
 
# ------------------------------------------------------------------
# Interpolate v → u-points  (Coriolis term in u-equation)
#
# The Coriolis term fv in the u-equation requires v at u-point
# locations. On the C-grid, the four v-points surrounding a
# u-point at (i,j) are at (i-1,j), (i,j), (i-1,j+1), (i,j+1):
#
#   v(i-1,j+1) -- v(i,j+1)
#    |      u(i,j)    |
#   v(i-1, j) --  v(i, j)
#
#   v̄_{i,j} = 0.25 * (v[i-1,j] + v[i,j] + v[i-1,j+1] + v[i,j+1])
# ------------------------------------------------------------------
"""
Interpolate v to u-points by four-point averaging. Output sized nx × ny.
"""
function interp_v_to_u!(out::Field, v::Field, g::Grid)
    o  = out.data
    vx = v.data
    nx, ny = g.nx, g.ny
 
    @inbounds for j ∈ 2:ny+1
        for i ∈ 2:nx+1
            o[i,j] = 0.25 * (vx[i-1,j] + vx[i,j] + vx[i-1,j+1] + vx[i,j+1])
        end
    end
    return nothing
end
 
# ------------------------------------------------------------------
# Interpolate u → v-points  (Coriolis term in v-equation)
#
# The Coriolis term fu in the v-equation requires u at v-point
# locations. The four u-points surrounding a v-point at (i,j) are
# at (i,j-1), (i+1,j-1), (i,j), (i+1,j):
#
#   u(i,j  )--v(i,j)--u(i+1,j  )
#   u(i,j-1)----------u(i+1,j-1)
#
#   ū_{i,j} = 0.25 * (u[i,j-1] + u[i+1,j-1] + u[i,j] + u[i+1,j])
# ------------------------------------------------------------------
"""
Interpolate u to v-points by four-point averaging. Output sized nx × (ny+1).
"""
function interp_u_to_v!(out::Field, u::Field, g::Grid)
    o  = out.data
    ux = u.data
    nx, ny = g.nx, g.ny
 
    @inbounds for j ∈ 2:ny+2
        for i ∈ 2:nx+1
            o[i,j] = 0.25 * (ux[i,j-1] + ux[i+1,j-1] + ux[i,j] + ux[i+1,j])
        end
    end
    return nothing
end
 
end # module operators