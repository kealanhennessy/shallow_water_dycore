module shallowwaterdycore

include("grid.jl")
include("fields.jl")
include("dynamics.jl")
include("timestepping.jl")
include("diagnostics.jl")
include("io.jl")

export run_model

function run_model(physical_params, numerical_params)
    grid = build_supergrid(physical_params)

    state = allocate_fields(grid)

    tendencies = allocate_fields(grid)

    state_stage = allocate_fields(grid)

    tendencies_stage = allocate_fields(grid)

    cache = allocate_cache(grid)

    t = 0.0
    while t < numerical_params.tmax
        advance_state!(tendencies, state, tendencies_stage, state_stage, cache, grid, physical_params, numerical_params)

        t += numerical_params.dt
    end # main time stepping loop

end # run_model

end