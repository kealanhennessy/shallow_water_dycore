module rhs
 
using ..fields
using ..grid
using ..operators
using ..parameters
 
export compute_tendencies!
 
# ------------------------------------------------------------------
# Halo filling
# Zonal periodicity: copy the last interior column into the west
# halo and the first interior column into the east halo.
# Applied to η, u, and v before each RHS evaluation.
#
# Index layout (nx=4, HALO=1):
#   1     | 2  3  4  5 | 6
#   west    interior     east
#   halo                 halo
# ------------------------------------------------------------------
function fill_zonal_halos!(q::Fields, g::SuperGrid)
    H = fields.HALO
    nx = g.nx
 
    for f ∈ (q.η.data, q.u.data, q.v.data)
        # West halo ← last interior column
        f[1,   :] .= f[nx+H, :]
        # East halo ← first interior column
        f[nx+H+1, :] .= f[H+1,  :]
    end
    return nothing
end
 
# ------------------------------------------------------------------
# Meridional halo filling
# The domain is bounded at the poles, so meridional halos are filled
# using a zero-gradient (Neumann) condition — the halo cell takes
# the value of the nearest interior cell. This is applied to η only,
# since:
#   - u stencils (∂u/∂λ, interp_u_to_v) are purely zonal or reach
#     one row south, which is covered by the v BC below
#   - v has explicit polar points at j=H+1 and j=ny+1+H so no halo
#     fill is needed
#
# The zero-gradient condition ensures ∂η/∂φ = 0 at the poles,
# which is physically consistent with v = 0 there.
#
# Index layout (ny=4, HALO=1):
#   j:   1       | 2  3  4  5 |  6
#        south     interior      north
#        halo                    halo
# ------------------------------------------------------------------
function fill_meridional_halos!(q::Fields, g::SuperGrid)
    H  = fields.HALO
    ny = g.ny
 
    # η south halo ← first interior row
    q.η.data[:, H]      .= q.η.data[:, H+1]
    # η north halo ← last interior row
    q.η.data[:, ny+1+H] .= q.η.data[:, ny+1]
    return nothing
end
 
# ------------------------------------------------------------------
# Polar boundary condition
# v = 0 at both poles. In the v-array (size nx × ny+1+2*HALO),
# the south pole is the first interior row (j = HALO+1 = 2) and
# the north pole is the last interior row (j = ny+1+HALO = ny+2).
# ------------------------------------------------------------------
function apply_polar_bc!(q::Fields, g::SuperGrid)
    H  = fields.HALO
    ny = g.ny
    q.v.data[:, H+1]    .= 0.0   # south pole
    q.v.data[:, ny+1+H] .= 0.0   # north pole
    return nothing
end
 
# ------------------------------------------------------------------
# compute_rhs!
#
# Evaluates the RHS of the linearized shallow water equations on a
# spherical lat-lon Arakawa C-grid:
#
#   ∂u/∂t =  f*v̄  -  (g / a cosφ) * ∂η/∂λ
#   ∂v/∂t = -f*ū  -  (g / a)       * ∂η/∂φ
#   ∂η/∂t = -(H / a cosφ) * (∂u/∂λ + ∂(v cosφ)/∂φ)
#
# Arguments:
#   dqdt  : tendency Fields (∂u/∂t, ∂v/∂t, ∂η/∂t), overwritten
#   q   : current state Fields (u, v, η)
#   c   : Cache of intermediate fields, overwritten
#   g   : Grid
#   p   : PhysicalParams
# ------------------------------------------------------------------
"""
Compute the RHS of the linearized shallow water equations, writing
tendencies into `dqdt`. Modifies `q` halos and `c` in place.
"""
function compute_tendencies!(dqdt::Fields, q::Fields, c::Cache, g::SuperGrid, pp::PhysicalParams)
 
    # 1. Fill halos for periodicity and poles, apply polar BC
    fill_zonal_halos!(q, g)
    fill_meridional_halos!(q, g)
    apply_polar_bc!(q, g)
 
    # 2. Compute spatial derivatives into cache
    ∂η_∂λ!(   c.dη_dλ,     q.η, g)
    ∂η_∂φ!(   c.dη_dφ,     q.η, g)
    ∂u_∂λ!(   c.du_dλ,     q.u, g)
    ∂vcosφ_∂φ!(c.dvcosφ_dφ, q.v, g)
 
    # 3. Interpolate for Coriolis
    interp_v_to_u!(c.v_at_u, q.v, g)
    interp_u_to_v!(c.u_at_v, q.u, g)
 
    # 4. Assemble tendencies
    nx, ny = g.nx, g.ny
    H = fields.HALO
 
    dudt = dqdt.u.data
    dvdt = dqdt.v.data
    dηdt = dqdt.η.data
 
    # u-tendency: fv̄ - (g / a cosφ) ∂η/∂λ
    @inbounds for j ∈ H+1:ny+H
        φⱼ        = g.u.φ[j-H]           # latitude at this u-row
        f         = 2p.Ω * sin(φⱼ)
        inv_a_cosφ = g.u.inv_a_cosφ[j-H]
        for i ∈ H+1:nx+H
            dudt[i,j] = f * c.v_at_u.data[i,j] - pp.g * inv_a_cosφ * c.dη_dλ.data[i,j]
        end
    end
 
    # v-tendency: -fū - (g / a) ∂η/∂φ
    @inbounds for j ∈ H+1:ny+H
        φⱼ  = g.v.φ[j-H]                 # latitude at this v-row
        f   = 2p.Ω * sin(φⱼ)
        inv_a = g.inv_a
        for i ∈ H+1:nx+H
            dvdt[i,j] = -f * c.u_at_v.data[i,j] - pp.g * inv_a * c.dη_dφ.data[i,j]
        end
    end
 
    # η-tendency: -(H / a cosφ) * (∂u/∂λ + ∂(v cosφ)/∂φ)
    @inbounds for j ∈ H+1:ny+H
        inv_a_cosφ = g.η.inv_a_cosφ[j-H]
        for i ∈ H+1:nx+H
            dηdt[i,j] = -pp.H * inv_a_cosφ * (c.du_dλ.data[i,j] + c.dvcosφ_dφ.data[i,j])
        end
    end
 
    return nothing
end
 
end # module rhs