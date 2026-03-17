module fields
 
using ..grid
 
export Field, Fields, allocate_fields
 
# ------------------------------------------------------------------
# Halo width
# One halo cell on each side of the domain is used to handle
# boundary conditions and periodic zonal wrapping without
# special-casing the stencil loops in operators.jl.
# ------------------------------------------------------------------
const HALO = 1
 
# ------------------------------------------------------------------
# Field
# A single 2D scalar field with physical units. The type parameter
# A is the underlying array type (e.g. Matrix{Float64}, CuArray),
# keeping this struct compatible with GPU backends without rewriting.
#
# Interior index ranges (i.e. excluding halos):
#   zonal (i)      : HALO+1 : nx+HALO      (nx points)
#   meridional (j) : HALO+1 : ny+HALO      (ny points, η and u)
#                    HALO+1 : ny+1+HALO    (ny+1 points, v)
# ------------------------------------------------------------------
struct Field{A<:AbstractMatrix{Float64}}
    data::A
    units::String
end
 
# ------------------------------------------------------------------
# Fields
# The complete model state vector q = (η, u, v) at one time level.
# Array sizes with halos:
#   η : (nx + 2*HALO) × (ny   + 2*HALO)
#   u : (nx + 2*HALO) × (ny   + 2*HALO)
#   v : (nx + 2*HALO) × (ny+1 + 2*HALO)
# ------------------------------------------------------------------
struct Fields
    η::Field
    u::Field
    v::Field
end

# ------------------------------------------------------------------
# Internal allocators
# ------------------------------------------------------------------
 
allocate_η(g::Grid) = Field(zeros(g.nx + 2*HALO, g.ny     + 2*HALO), "m")
allocate_u(g::Grid) = Field(zeros(g.nx + 2*HALO, g.ny     + 2*HALO), "m/s")
allocate_v(g::Grid) = Field(zeros(g.nx + 2*HALO, g.ny + 1 + 2*HALO), "m/s")
 
# ------------------------------------------------------------------
# Public constructors
# ------------------------------------------------------------------
 
"""
    allocate_fields(g::Grid) -> Fields
 
Allocate a zeroed `Fields` state vector sized to grid `g`,
including halo cells.
"""
function allocate_fields(g::Grid)
    return Fields(allocate_η(g), allocate_u(g), allocate_v(g))
end
 
"""
    Base.similar(f::Fields) -> Fields
 
Allocate a new `Fields` with the same array sizes and units as `f`
but with zeroed data. Used by the SSPRK3 integrator to allocate
intermediate stage arrays without needing to know the internal
structure of `Fields`.
"""
function Base.similar(f::Fields)
    return Fields(
        Field(similar(f.η.data), f.η.units),
        Field(similar(f.u.data), f.u.units),
        Field(similar(f.v.data), f.v.units),
    )
end
 
end # module fields