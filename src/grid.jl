module grid
 
export Grid, SuperGrid, build_grid
 
# ------------------------------------------------------------------
# SubGrid
# Holds the coordinates and precomputed metric weights for one of
# the three C-grid locations (η, u, v). The distinction between
# locations is carried by the field names in Grid (grid.η, grid.u,
# grid.v) and by the staggered coordinate values themselves, not
# by separate types.
#
# C-grid staggering (Arakawa C):
#   η : cell centers         (λ_i,        φ_j      )  nx    × ny
#   u : east/west faces      (λ_i + dλ/2, φ_j      )  nx    × ny
#   v : north/south faces    (λ_i,        φ_j - dφ/2) nx    × ny+1
#
# The v-grid has ny+1 points in φ, spanning both poles:
#   φ_v ∈ [-π/2, π/2]   (ny+1 points, including both poles)
#   φ_η ∈ (-π/2, π/2)   (ny points, cell centers)
#
# Coordinates:
#   λ ∈ [0, 2π)         longitude
#   φ ∈ [-π/2, π/2]     latitude
#
# Metric weights:
#   cosφ      : cos(φ), used to weight meridional fluxes
#   inv_a_cosφ: 1 / (a cosφ), used in zonal gradient / divergence
# ------------------------------------------------------------------
struct Grid
    λ::Vector{Float64}
    φ::Vector{Float64}
    cosφ::Vector{Float64}
    inv_a_cosφ::Vector{Float64}
end

# ------------------------------------------------------------------
# Grid
# Top-level grid struct. Holds global dimensions, uniform spacings,
# 1/a, and the three staggered subgrids.
#
# Array size summary (excluding halos):
#   η and u fields : nx × ny
#   v fields       : nx × (ny+1)
#
# The v-field has one extra row because v-faces bound every η cell
# from both sides in φ, including at both poles. The polar boundary
# condition v = 0 must be enforced at j = 1 (south pole) and
# j = ny+1 (north pole).
# ------------------------------------------------------------------
struct SuperGrid
    nx::Int
    ny::Int
    dλ::Float64
    dφ::Float64
    inv_a::Float64
    η::Grid
    u::Grid
    v::Grid
end
 
# ------------------------------------------------------------------
# Internal constructors
# ------------------------------------------------------------------
 
function build_η_grid(a::Float64, nx::Int, ny::Int)
    dλ = 2π / nx
    dφ = π  / ny
 
    # Cell centers: offset half a cell inward from each boundary
    λ = collect((0:nx-1) .* dλ)
    φ = collect(range(-π/2 + dφ/2, π/2 - dφ/2, length=ny))
 
    cosφ       = cos.(φ)
    inv_a_cosφ = 1.0 ./ (a .* cosφ)
 
    return Grid(λ, φ, cosφ, inv_a_cosφ), dλ, dφ
end

function build_u_grid(η::Grid, dλ::Float64)
    # u-points are offset half a cell east of η-points in λ.
    # φ is shared with η — no meridional offset.
    λ = η.λ .+ dλ/2
 
    return Grid(λ, η.φ, η.cosφ, η.inv_a_cosφ)
end

function build_v_grid(a::Float64, η::Grid)
    # v-points sit at cell faces in φ, bounding each η cell from
    # below and above. With ny η-points there are ny+1 v-points.
    # The range [-π/2, π/2] places v-points at both poles exactly,
    # where the boundary condition v = 0 is applied.
    # λ is shared with η — no zonal offset.
    φ = collect(range(-π/2, π/2, length=length(η.φ)+1))
 
    cosφ       = cos.(φ)
    inv_a_cosφ = 1.0 ./ (a .* cosφ)
 
    return Grid(η.λ, φ, cosφ, inv_a_cosφ)
end
 
# ------------------------------------------------------------------
# Public constructor
# ------------------------------------------------------------------

"""
Construct a spherical lat-lon Arakawa C-grid with `nx` zonal and
`ny` meridional cells. `a` is the planetary radius in metres.
"""
function build_supergrid(a::Float64, nx::Int, ny::Int)
    η, dλ, dφ = build_η_grid(a, nx, ny)
    u         = build_u_grid(η, dλ)
    v         = build_v_grid(a, η, dφ)
 
    return SuperGrid(nx, ny, dλ, dφ, 1.0/a, η, u, v)
end

end # module grid