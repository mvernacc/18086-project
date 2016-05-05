# Gas properties
#
# Matt Vernacchia
# 18.086 Project
# Spring 2016


export Gas, air, helium, R_univ

# Properties of a perfect gas.
type Gas
    # Ratio of specific heats [units: none].
    γ::Real
    # Heat capacity at constant pressure [units: joule kilogram^-1 kelvin^-1]
    c_p::Real
end

const air = Gas(1.40, 1005)
const helium = Gas(5/3, 5190)

# Universal gas constant [units: joule kelvin^-1 mole^-1]
const R_univ = 8.3144598
