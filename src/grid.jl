module grid

export build_grid

struct ηGrid
    λ::Vector{Float64}
    φ::Vector{Float64}
    cosφ::Vector{Float64}
end

struct UGrid
    λ::Vector{Float64}
    φ::Vector{Float64}
    cosφ::Vector{Float64}
end

struct VGrid
    λ::Vector{Float64}
    φ::Vector{Float64}
    cosφ::Vector{Float64}
end

struct Grid
    nx::Int
    ny::Int
    dλ::Float64
    dφ::Float64
    η::ηGrid
    u::UGrid
    v::VGrid
end

function build_ηgrid(nx, ny)
    dλ = 2π / nx
    λ = (0:nx-1) .* dλ

    φ = range(-π/2 + π/(2*ny), π/2 - π/(2*ny), length=ny)
    cosφ = cos.(φ)

    return ηGrid(λ, φ, cosφ), dλ
end

function build_ugrid(η::ηGrid, dλ)
    λ = η.λ .+ dλ/2
    φ = η.φ
    return UGrid(λ, φ, η.cosφ)
end

function build_vgrid(η::ηGrid, dφ)
    λ = η.λ
    φ = η.φ .+ dφ/2
    cosφ = cos.(φ)
    return VGrid(λ, φ, cosφ)
end

function build_grid(params)
    nx, ny = params.nx, params.ny

    η, dλ = build_ηgrid(nx, ny)
    dφ = η.φ[2] - η.φ[1]

    u = build_ugrid(η, dλ)
    v = build_vgrid(η, dφ)

    return Grid(nx, ny, dλ, dφ, η, u, v)
end

end