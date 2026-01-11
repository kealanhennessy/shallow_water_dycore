module grid

export build_grid

struct ηGrid
    λ::Vector{Float64}
    φ::Vector{Float64}
    cosφ::Vector{Float64}
    inv_a_cosφ::Vector{Float64}
end

struct UGrid
    λ::Vector{Float64}
    φ::Vector{Float64}
    cosφ::Vector{Float64}
    inv_a_cosφ::Vector{Float64}
end

struct VGrid
    λ::Vector{Float64}
    φ::Vector{Float64}
    cosφ::Vector{Float64}
    inv_a_cosφ::Vector{Float64}
end

struct Grid
    nx::Int
    ny::Int
    dλ::Float64
    dφ::Float64
    inv_a::Float64
    η::ηGrid
    u::UGrid
    v::VGrid
end

# η is free surface height, located at cell centers
# u is zonal velocity, located at east/west faces
# v is meridional velocity, located at north/south faces
# λ is longitude
# φ is latitude

function build_ηgrid(a, nx, ny)
    dλ = 2π / nx
    λ = (0:nx-1) .* dλ
    φ = range(-π/2 + π/(2*ny), π/2 - π/(2*ny), length=ny)
    
    cosφ = cos.(φ)
    inv_a_cosφ = 1.0 ./ (a * cosφ)

    return ηGrid(λ, φ, cosφ, inv_a_cosφ), dλ
end

function build_ugrid(η::ηGrid, dλ)
    λ = η.λ .+ dλ/2
    φ = η.φ

    cosφ = η.cosφ
    inv_a_cosφ = η.inv_a_cosφ

    return UGrid(λ, φ, cosφ, inv_a_cosφ)
end

function build_vgrid(a, η::ηGrid, dφ)
    λ = η.λ
    φ = η.φ .+ dφ/2
    
    cosφ = cos.(φ)
    inv_a_cosφ = 1.0 ./ (a * cosφ)
    
    return VGrid(λ, φ, cosφ, inv_a_cosφ)
end

function build_grid(params)
    nx, ny = params.nx, params.ny
    a = params.a
    inv_a = 1.0 / a

    η, dλ = build_ηgrid(a, nx, ny)
    dφ = η.φ[2] - η.φ[1]

    u = build_ugrid(η, dλ)
    v = build_vgrid(a, η, dφ)

    return Grid(nx, ny, dλ, dφ, inv_a, η, u, v)
end

end