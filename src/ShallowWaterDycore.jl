module shallowwaterdycore

include("grid.jl")
include("fields.jl")
include("dynamics.jl")
include("time_integration.jl")
include("boundary_conditions.jl")
include("diffusion.jl")
include("diagnostics.jl")
include("io.jl")

export run_model

function run_model(params)
    grid = build_grid(params)

    state = init_state(params)

    t = 0.0
    while t < params.tmax
        tendencies = compute_tendencies(grid, state, params)

        state = advance_state(grid, state, params, tendencies)

        apply_boundary_conditions!(grid, state params)

        t += params.dt
    end # main time stepping loop

end # run_model