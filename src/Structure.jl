module Structure

using SpecialFunctions: gamma

export dmdr, dPdr
export planet_m, planet_g, planet_mmotion, planet_mu, planet_cmu

"""
    dmdr(r::Real, ρ::Real)

Mass in a spherical shell.

# Arguments
- `r::Real` - radius
- `ρ::Real` - density
"""
dmdr(r::Real, ρ::Real) = 4π * ρ * r^2

"""
    dPdr(r::Real, ρ::Real, gm::Real)

Pressure due to hydrostatic equilibrium. The gravitational mass is Newton's
Gravitational constant multiplied by the mass of the body.

# Arguments
- `r::Real`  - radius
- `ρ::Real`  - density
- `gm::Real` - gravitational mass
"""
dPdr(r::Real, ρ::Real, gm::Real) = -gm * ρ / r^2
dPdr(ρ::Real, g::Real) = -ρ * g

"""
    planet_g(r::Real, gm::Real)

Gravity due to a spherical mass.

# Arguments
- `r::Real`  - radius
- `gm::Real` - gravitational mass
"""
planet_g(r::Real, gm::Real) = r == 0 ? 0 : gm / r^2

"""
    planet_mmotion(gm::Real, a::Real)

Mean motion due to a planetary body.

# Arguments
- `gm::Real` - gravitational mass
- `a::Real`  - radius of the orbit
"""
planet_mmotion(gm::Real, a::Real) = sqrt(gm / a^3)

"""
    planet_mu(r::Real, g::Real, ρ::Real)

Veritcal stress due to gravity.

# Arguments
- `r::Real` - radius
- `g::Real` - gravity
- `ρ::Real` - density
"""
planet_mu(r::Real, g::Real, ρ::Real) = ρ * g * r


"""
    cmu_maxwell(μ::Real, ω::Real, η::Real)

Maxwell rheaology using the tidal forcing frequency. Refer to Turcotte &
Schubert 2002 or Storch & Lai 2014 for a decription of this rheaology.

# Arguments
- `μ::Real` - shear modulus (rigidity)
- `ω::Real` - rotational frequency
- `η::Real` - viscosity
"""
function cmu_maxwell(μ::Real, ω::Real, η::Real)

    ω_m = μ / η

    cmu = μ / (1 - (ω_m/ω)*im)

    return cmu

end

"""
    cmu_SLS(μ0::Real, ω::Real, η::Real, μ_f::Real)

Standard Linear Solid model. Refer to Stixrude et al. 2021 or Norwick & Berry
1972 for a description of this rheaology.

# Arguments
- `μ0::Real`  - unrelaxed shear modulus
- `ω::Real`   - rotational frequency
- `η::Real`   - viscosity
- `μ_f::Real` - SLS parameter s.t. μ₁/μ₀
"""
function cmu_SLS(μ0::Real, ω::Real, η::Real, μ_f::Real)

    μ1 = μ0 * μ_f

    dμ = μ0 * (1 - μ1 / (μ0 + μ1))

    τ = η / (μ0 + μ1)

    cmu = μ0 - dμ / (1 + τ*ω*im)

end

"""
    cmu_andrade(μ::Real, ω::Real, η::Real, α::Real; β::Real=0)

Andrade model. Refer to Castillo-Rogez et al. 2011 or Dumoulin et al. 2017 for a
description of this rheaology.

# Arguments
- `μ::Real` - shear modulus
- `ω::Real` - rotational frequency
- `η::Real` - viscosity
- `α::Real` - frequency dependence (Andrade exponent)
- `β::Real` - amplitude of the transient response
"""
function cmu_andrade(μ::Real, ω::Real, η::Real, α::Real; β::Real=0)

    if(β == 0) β = μ^(α - 1) * η^-α end

    J = 1 / μ - im / (η * ω) + β*(ω*im)^-α * gamma(1 + α)

    cmu = 1 / J

end

"""
    planet_cmu(μ::Real, ω::Real, η::Real, r::Real, g::Real, ρ::Real,
               μ_f::Real, model::Int)

Calculates the complex shear modulus (CMU) based on model input desire.
Viscosity will determine if the CMU model is used or if the infinite limit is
used. If viscosity is 0.0, a small shear modulus is produced for stability.
Refer to `cmu_maxwell`, `cmu_SLS`, or `cmu_andrade` for infromation on the
models.

# Arguments
- `μ::Real`    - shear modulus
- `ω::Real`    - rotational frequency
- `η::Real`    - viscosity
- `r::Real`    - radius
- `g::Real`    - gravity
- `ρ::Real`    - density
- `μ_f::Tuple` - rheaology parameters ([1] SLS factor, [2,3] Andrade alpha, beta)
- `model::Int` - model desire ([1] Maxwell, [2] SLS, [3] Andrade)
"""
function planet_cmu(μ::Real, ω::Real, η::Real, r::Real, g::Real, ρ::Real,
                    μ_f::Tuple, model::Int)

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
            cmu = cmu_SLS(μ, ω, η, μ_f[1])
        elseif model == 3
            cmu = cmu_andrade(μ, ω, η, μ_f[2], β=μ_f[3])
        end

    end

    if(real(cmu) == 0.0)
        # final check if shear is 0.0

        cmu = 1e-9 + 0*im

    end

    return cmu

end

end # module
