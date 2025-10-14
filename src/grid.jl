module grid

export build_grid

struct Grid
    nx::Int           # number of longitude points
    ny::Int           # number of latitude points
    lon::Array{Float64,1}  # longitude array [rad]
    lat::Array{Float64,1}  # latitude array [rad]
    dlon::Float64     # grid spacing in longitude
    dlat::Float64     # grid spacing in latitude
    coslat::Array{Float64,1} # cos(lat) for metric factors
    cell_area::Array{Float64,2} # area of each cell
end 

function build_grid(params)
    nx, ny = params.nx, params.ny

    lon = range(0, 2*pi, length=ny)[begin:end-1] # chop of the end (periodic)
    lat = range(-pi/2 + pi/(2*ny), pi/2 - pi/(2*ny), length=ny) # chop off the ends (poles)

    dlon = lon[2] - lon[1]
    dlat = lat[2] - lat[1]

    coslat = cos.(lat)

    a = params.a
    cell_area = zeros(nx, ny)
    for i = 1:nx
        for j = 1:ny
            cell_area[i,j] = (a^2) * dlon * dlat * coslat[j]
        end
    end

    return Grid(nx, ny, lon, lat, dlon, dlat, coslat, cell_area)
end