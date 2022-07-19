module Structure

using SpecialFunctions: gamma

export dmdr, dPdr
export planet_m, planet_g, planet_mmotion, planet_mu, planet_eta, planet_cmu, planet_structure

# Governing Equations
dmdr(r::Real, ρ::Real) = 4π * ρ * r^2
dPdr(r::Real, ρ::Real, gm::Real) = -gm * ρ / r^2
dPdr(ρ::Real, g::Real) = -ρ * g

# Planetary Structure EQNs
planet_m(r::Real, ρ::Real) = 4π * ρ * r^2
planet_g(r::Real, gm::Real) = r == 0 ? 0 : gm / r^2
planet_mmotion(gm::Real, a::Real) = sqrt(gm / a^3)
planet_mu(r::Real, g::Real, ρ::Real) = ρ * g * r
planet_mu(T::Real, P::Real) = (101 + P*1e-9/2.41 - (T-1650)/23.4) * 1e9


# Rheaology Models
function cmu_maxwell(μ::Real, ω::Real, η::Real)

    ω_m = μ / η

    cmu = μ / (1 - (ω_m/ω)*im)

    return cmu

end

function cmu_SLS(μ0::Real, ω::Real, η::Real, μ_f::Real)

    μ1 = μ0 * μ_f

    dμ = μ0 * (1 - μ1 / (μ0 + μ1))

    τ = η / (μ0 + μ1)

    cmu = μ0 - dμ / (1 + τ*ω*im)

end

function cmu_andrade(μ::Real, ω::Real, η::Real, α::Real; β::Real=0)

    if(β == 0) β = μ^(α - 1) * η^-α end

    J = 1 / μ - im / (η * ω) + β*(ω*im)^-α * gamma(1 + α)

    cmu = 1 / J

end


function planet_cmu(μ::Real, ω::Real, η::Real, r::Real, g::Real,
                    ρ::Real, μ_f::Real, model::Int)

    if η == 0.0
        # if eta is small, produce small shear

        smu = 1e-4 * planet_mu(r, g, ρ)
        cmu = smu + 0*im

    elseif η == -1
        # if eta is infinite, produce the real shear

        cmu = μ + 0*im

    else
        # else use the desired rheaology

        if model == 1
            cmu = cmu_maxwell(μ, ω, η)
        elseif model == 2
            cmu = cmu_SLS(μ, ω, η, μ_f)
        elseif model == 3
            cmu = cmu_andrade(μ, ω, η, μ_f)
        end

    end

    if(real(cmu) == 0.0)
        # final check if shear is 0.0

        cmu = 1e-9 + 0*im

    end

    return cmu

end

end # module
