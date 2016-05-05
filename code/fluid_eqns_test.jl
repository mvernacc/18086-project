# Unit tests for fluid_eqns.jl
#
# Matt Vernacchia
# 18.086 Project
# Spring 2016

using Base.Test

include("fluid_eqns.jl")

air = Gas(1.40, 1005)
helium = Gas(5/3, 5190)

function test_u2e()
    @test u2e([1, 0, 0, 2]) == 2
    @test u2e([1, 2, 4, 0]) == -10
    @test u2e([10, 0, 0, 1]) == 0.1
end

function test_u2p()
    # if u = v = 0, p = (γ - 1) * E
    @test u2p([234, 0, 0, 101e3 / 0.40], air) ≈ 101e3
    @test u2p([33e-4, 0, 0, 1 / 0.40], air) ≈ 1
    @test u2p([33e-4, 0, 0, 1 / (2/3)], helium) ≈ 1

    # E = p/(γ - 1) + ρ/2 (u^2 + v^2)
    @test u2p([2, 2, 2, 1 / (2/3) + 2], helium) ≈ 1
end

function test_u2T()
end
    

test_u2e()
test_u2p()
println("Passed all tests.")
