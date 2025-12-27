module grid

export build_grid

struct Grid
    nx::Int           # number of longitude points
    ny::Int           # number of latitude points
    λ::Array{Float64,1}  # longitude array [rad]
    φ::Array{Float64,1}  # latitude array [rad]
    dλ::Float64     # grid spacing in longitude
    dφ::Float64     # grid spacing in latitude
    cosφ::Array{Float64,1} # cos(lat) for metric factors
    cell_area::Array{Float64,2} # area of each cell
end 

function build_grid(params)
    nx, ny = params.nx, params.ny

    λ = range(0, 2*π, length=ny)[begin:end-1] # chop of the end (periodic)
    φ = range(-π/2 + π/(2*ny), π/2 - π/(2*ny), length=ny) # chop off the ends (poles)

    dλ = λ[2] - λ[1]
    dφ = φ[2] - φ[1]

    cosφ = cos.(φ)

    a = params.a
    cell_area = zeros(nx, ny)
    cell_area[i,j] .= a^2 * dλ * dφ * cosφ'

    return Grid(nx, ny, λ, φ, dλ, dφ, cosφ, cell_area)
end

end