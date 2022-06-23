module Structure

using SpecialFunctions: gamma
using PhysicalConstants.CODATA2014: G
using QuadGK: quadgk

export dmdr, dPdr
export planet_m, planet_g, planet_mmotion, planet_mu, planet_eta, planet_cmu, planet_structure

# Governing Equations
dmdr(r::Real, ρ::Real) = 4π * ρ * r^2
dPdr(r::Real, ρ::Real, m::Real) = -G.val * m * ρ / r^2
dPdr(ρ::Real, g::Real) = -ρ * g

# Planetary Structure EQNs
planet_m(r::Real, ρ::Real) = 4π * ρ * r^2
planet_g(r::Real, m::Real) = r == 0 ? 0 : G.val * m / r^2
planet_mmotion(m::Real, a::Real) = sqrt(G.val * m / a^3)
planet_mu(r::Real, g::Real, ρ::Real) = ρ * g * r

planet_eta(η0::Real, A::Real, T_m::Real, T::Real) = η0 * exp(A * (T_m/T - 1))

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

function cmu_andrade(μ::Real, ω::Real, η::Real, α::Real)

    β = μ^(α - 1) * η^-α

    J = 1 / μ - im / (η * ω) + β*(ω*im)^-α * gamma(1 + α)

    cmu = 1 / J

end


function planet_cmu(μ::Real, ω::Real, η::Real, r::Real, g::Real,
                    ρ::Real, μ_f::Real, model::Int)

    if η == 0.0
        # if eta is small, produce small shear

        smu = 1e-5 * planet_mu(r, g, ρ)
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

function planet_structure(plnt, data::Matrix)

    mass = 0.0
    r0 = 0.0
    layers = plnt.layers
    ω = plnt.ω
    μ_f = plnt.μ_f
    model = plnt.rhea_model

    sd = zeros(Complex, layers, 4)

    for i in 1:layers

        r1 = data[i,1]
        ρ = data[i,2]
        μ = data[i,3]
        η = data[i,4]

        if(i != 1) r0 = real(sd[i-1,1]) end
        res, err = quadgk(x -> planet_m(x, ρ), r0, r1)
        mass += res

        g = planet_g(r1, mass)

        sd[i,4] = ρ
        sd[i,3] = g
        sd[i,2] = planet_cmu(μ, ω, η, r1, g, ρ, μ_f, model)
        sd[i,1] = r1

    end

    return sd, mass

end

end # module
