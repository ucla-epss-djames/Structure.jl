module Structure

include("Planets.jl")

using .Planets

import SpecialFunctions : gamma
import PhysicalConstants.CODATA2018 : G
import QuadGK : quadgk

# Governing Equations
dmdr(r::Real, ρ::Real) = 4π * ρ * r^2
dPdR(r::Real, ρ::Real, m::Real) = -G * m * ρ / r^2

# Planetary Structure EQNs
planet_m(r::Real, ρ::Real) = 4π * ρ * r^2
planet_g(r::Real, m::Real) = r == 0 ? return 0 : return G * m / r^2

planet_mmotion(m::Real, a::Real) = sqrt(G * m / a**3)

# Rheaology Models
function cmu_maxwell(μ::Real, ω::Real, η::Real)

    ω_m = μ / η

    cmu = μ / (1 + (ω_m / ω)*im)

    return cmu

end

function cmu_SLS(μ0::Real, ω::Real, η::Real, μ_f::Real)

    μ1 = μ0 * μ_f

    dμ = μ0 * (1 - μ1 / (μ0 + μ1))

    τ = η / (μ0 + μ1)

    cmu = μ0 - dμ / (1 - τ*ω*im)

end

function cmu_andrade(μ::Real, ω::Real, η::Real, α::Real, n::Real)

    β = μ^(α - 1) * η^-α
    χ = 2*(ω - n)

    cmu = 1 / χ - im / (η * χ) + β*(χ*im)^-α * gamma(1 + α)

end



planet_mu(r::Real, g::Real, ρ::Real) = ρ * a * r

function planet_cmu(μ::Real, ω::Real, η::Real, r::Real, g::Real,
                    ρ::Real, μ_f::Real, α::Real, n::Real, model::Int)

    if η == 0.0
        # if eta is small, produce small shear

        smu = 1e-4 * planet_mu(r, g, ρ)
        cmu = smu + 0*im

    elseif η == Inf
        # if eta is infinite, produce the real shear

        cmu = μ + 0*im

    else
        # else use the desired rheaology

        if model == 1
            cmu = cmu_maxwell(μ, ω, η)
        elseif model == 2
            cmu = cmu_SLS(μ, ω, η, μ_f)
        elseif model == 3
            cmu = cmu_andrade(μ, ω, η, α, n)
        end

    end

    if(real(cmu) == 0.0)
        # final check if shear is 0.0

        cmu = 1e-9 + 0*im

    end

    return cmu

end

function planet_structure(plnt::Planet, data::Matrix{Real},
                          mass::Real, sd::Matrix{Complex})

    mass = 0.0
    layers = plnt.layers
    ω = plnt.ω
    α = plnt.α
    μ_f = plnt.μ_f
    n = plnt.n
    model = plnt.rhea_model

    for i in 1, layers

        ρ = real(sd[i,4])
        μ = data[i,1]
        η = data[i,2]

        r0 = 0.0
        r1 = real(sd[i,1])

        if i == 1
            mass, err = quadgk(x -> planet_m(x, ρ), r0, r1)
        else
            r0 = real(sd[i-1,1])
            res, err = quadgk(x -> planet_m(x, ρ), r0, r1)
            mass += res
        end

        g = planet_g(r1, mass)
        sd[i,3] = g
        sd[i,2] = planet_cmu(μ, ω, η, r1, g, ρ, μ_f, α, n, model)

    end

end

end # module
