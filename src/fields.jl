module fields

using ..grid

struct Field
    data::Array{float64, 2}
    name::Symbol
    units::String
end

function allocate_field(g::Grid; name::Symbol=:unnamed, units::String="" )
    data = zeros(g.nx, g.ny)
    return Field(
        data,
        name,
        units
    )
end 

import Base: +, -, *, /, .=

+(a::Field, b::Field) = Field(a.data .+ b.data, name=:a.name, units=a.units)
-(a::Field, b::Field) = Field(a.data .- b.data, name=:a.name, units=a.units)
*(x::Number, b::Field) = Field(x .* b.data, name=:b.name, units=b.units)
/(a::Field, y::Number) = Field(a.data ./ y, name=a.name, units=a.units)

function copyField!(dest::Field, src::Field)
    dest.data .= src.data
    return dest
end

function init_fields(g::Grid)
    h = allocate_field(g, name=:h, units="m")
    u = allocate_field(g, name=:u, units="m/s")
    v = allocate_field(g, name=:v, units="m/s")
    return (h=h, u=u, v=v)
end

end