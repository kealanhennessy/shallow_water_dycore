module parameters

export default_params

struct swparams
    a::float64
    omega::float64
    g::float64

    nx::int
    ny::int
    dt::float64
    tmax::float64

    h::float64
    viscosity::float64
    cfl::float64
    write_interval::float64
end # struct swparams

function default_params()
    return swparams(
        a = 6.37122e6,
        Ω = 7.2921e-5,
        g = 9.80616,
        nx = 96,
        ny = 48,
        dt = 300.0,
        tmax = 3600.0*10,  # 10 hours
        H = 1000.0,
        viscosity = 1e7,   # example biharmonic coefficient
        write_interval = 1800.0,
        CFL = 0.25
    )
end

end # module