module Planets

struct Planet

    # structure parameters
    layers::Int
    Ri::Real
    Rf::Real

    # orbital parameters
    ω::Real
    n::Real

    # rheaology parameters
    η::Real
    α::Real
    μ_f::Real
    rhea_model::Real

end

end # module
