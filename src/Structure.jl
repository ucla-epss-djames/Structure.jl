module Structure

import PhysicalConstants.CODATA2018: G

# Governing Equations
dmdr(r::Real, ρ::Real) = 4π * ρ * r^2
dPdR(r::Real, ρ::Real, m::Real) = -G * m * ρ / r^2

# Planetary Structure EQNs
planet_m(r::Real, ρ::Real) = 4π * ρ * r^2
planet_g(r::Real, m::Real) = r == 0 ? return 0 : return G * m / r^2

planet_mmotion(m::Real, a::Real) = sqrt(G * m / a**3)

planet_mu(r::Real, g::Real, ρ::Real) = ρ * a * r

function planet_cmu(μ::Real, ω::Real, η::Real, r::Real, g::Real, ρ::Real)

    if η == 0.0
        # if eta is small, produce small shear so Tidal.jl isn't unstable

        smu = 1e-4 * planet_mu(r, g, ρ)
        cmu = smu + 0*im

    elseif η == Inf
        # if eta is infinite, produce the real shear

        cmu = μ + 0*im

    else
        # else use the desired rheaology

    end

    if(real(cmu) == 0.0)
        # final check if shear is 0.0

        cmu = 1e-9 + 0*im

    end

    return cmu

end




end # module
