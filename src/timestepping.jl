module timestepping

using ..fields
using ..grid
using ..dynamics
using ..parameters

#
# Implements the SSRK3 time stepping algorithm
#

function ssprk3_stage!(q_out::Fields, q::Fields, q_stage::Fields, dqdt_stage::Fields, α::Float64, β::Float64, dt::Float64)

    q_out.u.data .= α * q.u.data .+ β * (q_stage.u.data .+ dt .* dqdt_stage.u.data)
    q_out.v.data .= α * q.v.data .+ β * (q_stage.v.data .+ dt .* dqdt_stage.v.data)
    q_out.η.data .= α * q.η.data .+ β * (q_stage.η.data .+ dt .* dqdt_stage.η.data)
    
    return nothing

end

"""
Update the current state q using the tendencies dqdt.
"""
function advance_state!(dqdt::Fields, q::Fields, dqdt_stage::Fields, q_stage::Fields, c::Cache, g::SuperGrid, pp::PhysicalParams, np::NumericalParams)

    compute_tendencies!(dqdt, q, c, g, pp)

    α = 0
    β = 1

    ssprk3_stage!(q_stage, q, q_stage, dqdt, α, Β, np.dt)

    compute_tendencies!(dqdt_stage, q_stage, c, g, pp)

    α = 3/4
    β = 1/4

    ssprk3_stage!(q_stage, q, q_stage, dqdt_stage, α, β, np.dt) 

    compute_tendencies!(dqdt_stage, q_stage, c, g, pp)

    α = 1/3
    β = 2/3

    ssprk3_stage!(q, q, q_stage, dqdt_stage, α, β, np.dt)

    return nothing

end

end # module timestepping