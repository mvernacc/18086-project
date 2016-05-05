"""
Fluid dynamics equations.

Matt Vernacchia
18.086 Project
Spring 2016
"""

# Properties of a perfect gas.
type Gas
    # Ratio of specific heats [units: none].
    γ::Real
    # Heat capacity at constant pressure [units: joule kilogram^-1 kelvin^-1]
    c_p::Real
end

# Universal gas constant [units: joule kelvin^-1 mole^-1]
const R_univ = 8.3144598

# Functions for unpacking / packing the vector of conservative variables.
# Numerical methods for the Navier-Stokes equations typically solve for the
# vector variable U:
#   U = [ρ, ρ u, ρ v, E]
# where:
#   ρ is density [units: kilogram meter^-3]
#   u is x velocity [units: meter second^-1]
#   v is y velocity [units: meter second^-1]
#   E is total energy per unit volume [units: joule meter^-3 = pascal].


@doc """
Get the density from the U vector.

Returns:
    desnity [units: kilogram meter^-3].
""" ->
function u2ρ(U)
    return U[1]
end


@doc """
Get the internal energy from the U vector.

Returns:
    specific internal energy [units: joule kilogram^-1].
""" ->
function u2e(U)
    u = U[2] / U[1]
    v = U[3] / U[1]
    e = U[4] / U[1] - 0.5 * (u^2 + v^2)
    return e
end


@doc """
Get the static pressure from the U vector.

Returns:
    pressure [units: pascal].
""" ->
function u2p(U, gas::Gas)
    # Ideal gas law
    #   p = ρ (γ - 1) e
    # Source: https://en.wikipedia.org/wiki/Equation_of_state#Classical_ideal_gas_law
    return U[1] * (gas.γ - 1) * u2e(U)
end


@doc """
Get the absolute static temperature from the U vector.

Returns:
    temperature [units: kelvin].
""" ->
function u2T(U, gas::Gas)
    # Use the definition of heat capacity for a calorically perfect gas.
    #   e = c_v T
    return gas.γ / gas.c_p * u2e(U)
end


@doc """
Get the Mach number from the U vector.
""" ->
function u2M(U, gas::Gas)
    # Speed of sound
    a = (gas.γ * u2p(U, gas) / U[1])^0.5
    # Velocities    
    u = U[2] / U[1]
    v = U[3] / U[1]
    return (u^2 + v^2)^0.5 / a
end


@doc """
Get the stagnation pressure from the U vector.
""" ->
function u2po(U, gas::Gas)
    # Source: https://en.wikipedia.org/wiki/Stagnation_pressure#Compressible_flow
    return u2p(U, gas) * (1 + (gas.γ - 1) / 2 * u2M(U, gas)^2)^(gas.γ / (gas.γ - 1))
end


@doc """
Get the stagnation temperature from the U vector.
""" ->
function u2To(U, gas::Gas)
    # Source: https://en.wikipedia.org/wiki/Stagnation_temperature
    return u2T(U, gas) * (1 + (gas.γ - 1) / 2 * u2M(U, gas)^2)
end


@doc """
Find the U vector from pressure, temperature, and velocity.

Arguments:
    p: static pressure [units: pascal]
    T: static temperature [units: kelvin]
    u: x velocity [units: meter second^-1]
    v: y velocity[units: meter second^-1]
    gas: Gas properties.
""" ->
function pTvel2u(p, T, u, v, gas::Gas)
    e = gas.c_p / gas.γ * T
    ρ = p / (gas.γ - 1) / e
    U = [ρ, ρ * u, ρ * v, ρ * (e + (u^2 + v^2) / 2)]
end


@doc """
Find the U vector from static pressure, static temperature and
Mach number.

Arguments:
    p: static pressure [units: pascal]
    T: static temperature [units: kelvin]
    M_x: x direction Mach number [units: none]
    M_y: y direction Mach number [units: none]
    gas: Gas properties
""" ->
function pTM2u(p, T, M_x, M_y, gas)
    # Mach number
    M = (M_x^2 + M_y^2)^0.5
    # Speed of sound
    R = gas.c_p * (1 - 1 / gas.γ)
    a = (gas.γ * R * T)^0.5
    u = M_x * a
    v = M_y * a
    return pTvel2u(p, T, u, v, gas)
end


@doc """
Find the U vector from stagnation pressure, stagnation temperature and
Mach number.

Arguments:
    p_o: stagnation pressure [units: pascal]
    T_o: stagnation temperature [units: kelvin]
    M_x: x direction Mach number [units: none]
    M_y: y direction Mach number [units: none]
    gas: Gas properties
""" ->
function poToM2u(p_o, T_o, M_x, M_y, gas)
    # Mach number
    M = (M_x^2 + M_y^2)^0.5
    # Static pressure
    p = p_o / (1 + (gas.γ - 1) / 2 * M^2)^(gas.γ / (gas.γ - 1))
    # Static temperature
    T = T_o / (1 + (gas.γ - 1) / 2 * M^2)
    return pTM2u(p, T, M_x, M_y, gas)
end


# Flux functions.
@doc """
The Euler (no viscosity or thermal conductivity) flux in the x direction.
""" ->
function F_euler(U, gas::Gas)
    u = U[2] / U[1]
    v = U[3] / U[1]
    p = u2p(U, gas)
    return [U[2], U[2] * u + p, U[3] * u, (U[4] + p) * u]
end


@doc """
The Euler (no viscosity or thermal conductivity) flux in the y direction.
""" ->
function G_euler(U, gas::Gas)
    u = U[2] / U[1]
    v = U[3] / U[1]
    p = u2p(U, gas)
    return [U[3], U[2] * v, U[3] * v + p, (U[4] + p) * v]
end

