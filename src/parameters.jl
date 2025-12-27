module parameters

export default_params

struct swparams
    a::Float64
    Ω::Float64
    g::Float64

    nx::Int
    ny::Int
    dt::Float64
    tmax::Float64

    h::Float64
    ν::Float64
    cfl::Float64
    write_interval::Float64
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
        h = 1000.0,
        ν = 1e7,   # example biharmonic coefficient
        cfl = 0.25,
        write_interval = 1800.0,
    )
end

end # module