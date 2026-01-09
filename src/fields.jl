module fields
using ..grid

export Fields, allocate_fields

struct Field{T,N}
    data::Array{T,N}
    name::Symbol
    units::String
end

struct Fields
    η::Field{Float64,2}
    u::Field{Float64,2}
    v::Field{Float64,2}
end

# allocate using each sub-grid’s coordinate lengths
allocate_η(g::Grid) = Field(zeros(length(g.η.λ), length(g.η.φ)), :η, "m")
allocate_u(g::Grid) = Field(zeros(length(g.u.λ), length(g.u.φ)), :u, "m/s")
allocate_v(g::Grid) = Field(zeros(length(g.v.λ), length(g.v.φ)), :v, "m/s")

function allocate_fields(g::Grid)
    return Fields(allocate_η(g), allocate_u(g), allocate_v(g))
end

end