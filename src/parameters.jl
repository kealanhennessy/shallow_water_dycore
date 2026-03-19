module parameters
 
export PhysicalParams, NumericalParams, default_physical, default_numerical
 
# ------------------------------------------------------------------
# PhysicalParams
# Physical constants and mean-state parameters.
#
# Fields:
#   a  : planetary radius (m)
#   Ω  : planetary rotation rate (rad/s)
#   g  : gravitational acceleration (m/s²)
#   H  : mean fluid depth for linearisation (m)
# ------------------------------------------------------------------
Base.@kwdef struct PhysicalParams
    a::Float64 = 6.37122e6
    Ω::Float64 = 7.2921e-5
    g::Float64 = 9.80616
    H::Float64 = 1000.0
end
 
# ------------------------------------------------------------------
# NumericalParams
# Time integration and I/O configuration.
# Grid dimensions (nx, ny) are excluded — they are arguments to
# build_grid and do not belong alongside physical constants.
#
# Fields:
#   dt             : timestep (s)
#   tmax           : total integration time (s)
#   cfl            : target CFL number (used for diagnostics / checks)
#   write_interval : output frequency (s)
# ------------------------------------------------------------------
Base.@kwdef struct NumericalParams
    dt::Float64             = 300.0
    tmax::Float64           = 3600.0 * 12   # 12 hours
    cfl::Float64            = 0.25
    write_interval::Float64 = 1800.0
    nx::Int                 = 144
    ny::Int                 = 72
end
 
# ------------------------------------------------------------------
# Default constructors
# Calling with no arguments returns Earth / standard test values.
# Individual fields can be overridden via keyword arguments:
#
#   p = default_physical(g = 1.62)   # Moon gravity
#   n = default_numerical(dt = 60.0)
# ------------------------------------------------------------------
 
"""
Return Earth-like physical parameters. Any field can be overridden
by keyword argument.
"""
default_physical(; kwargs...) = PhysicalParams(; kwargs...)
 
"""
Return default numerical parameters. Any field can be overridden
by keyword argument.
"""
default_numerical(; kwargs...) = NumericalParams(; kwargs...)
 
end # module parameters