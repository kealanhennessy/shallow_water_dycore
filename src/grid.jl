module grid

export build_grid

struct Grid
    nx::Int           # number of longitude points
    ny::Int           # number of latitude points
    
    λc::Array{Float64,1}  # longitude array [rad]
    φc::Array{Float64,1}  # latitude array [rad]
    λu::Array{Float64,1}  # longitude array [rad]
    φu::Array{Float64,1}  # latitude array [rad]
    λv::Array{Float64,1}  # longitude array [rad]
    φv::Array{Float64,1}  # latitude array [rad]
    
    dλ::Float64     # grid spacing in longitude
    dφ::Float64     # grid spacing in latitude
    
    cosφc::Array{Float64,1} # cos(lat) for metric factors
    cosφu::Array{Float64,1} # cos(lat) for metric factors
    cosφv::Array{Float64,1} # cos(lat) for metric factors
    
    cell_area::Array{Float64,2} # area of each cell
end 

function build_grid(params)
    nx, ny = params.nx, params.ny

    λc = range(0, 2*π, length=ny)[begin:end-1] # chop of the end (periodic)
    φc = range(-π/2 + π/(2*ny), π/2 - π/(2*ny), length=ny) # chop off the ends (poles)

    dλ = λc[2] - λc[1]
    dφ = φc[2] - φc[1]

    λu = λc .+ dλ/2
    φu = φc

    λv = λc
    φv = φc .+ dφ/2

    cosφc = cos.(φc)
    cosφu = cos.(φu)
    cosφv = cos.(φv)

    return Grid(nx, ny, dλ, dφ, λc, φc, λu, φu, λv, φv, cosφc, cosφu, cosφv)
end

end