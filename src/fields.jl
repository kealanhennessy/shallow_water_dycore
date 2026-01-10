module fields
using ..grid

export Fields, allocate_fields

struct Field{T,N}
    data::Array{T,N}
    units::String
end

struct Fields
    η::Field{Float64,2}
    u::Field{Float64,2}
    v::Field{Float64,2}
end

# use one cell 
const halo = 1

# allocate using halo's
# actual values stored in 2:nx-1, 2:ny-1
allocate_η(grid::Grid) = Field(zeros(grid.nx + 2*halo, grid.ny + 2*halo), "m")
allocate_u(grid::Grid) = Field(zeros(grid.nx + 2*halo, grid.ny + 2*halo), "m/s")
allocate_v(grid::Grid) = Field(zeros(grid.nx + 2*halo, grid.ny + 2*halo), "m/s")

function allocate_fields(grid::Grid)
    return Fields(allocate_η(grid), allocate_u(grid), allocate_v(grid))
end

end