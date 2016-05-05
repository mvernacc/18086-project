# Unit tests for solver.jl
#
# Matt Vernacchia
# 18.086 Project
# Spring 2016

using Base.Test

include("solver.jl")
include("bound.jl")

air = Gas(1.40, 1005)
helium = Gas(5/3, 5190)

function test_unif_open()
    # A uniform grid with open bounds should remain uniform and constant.
    srand(24324)
    U = 100 * rand() * ones(8, 8, 4)
    U_original = copy(U)
    ps = ProblemSpec(air, 1., 1., 1.,
        x -> x,
        x -> x,
        x -> x,
        x -> x)

    @test U ≈ MacCormack_step(U, ps)

    for i in 1:100
        U = MacCormack_step(U, ps)
        @test U ≈ U_original
    end
end

function test_unif_wall()
    # A uniform grid with walls should remain uniform and constant,
    # if the velocity is zero.
    srand(337)
    U = 100 * rand() * ones(8, 8, 4)
    U[:, :, 2:3] = 0
    U_original = copy(U)

    ps = ProblemSpec(air, 1., 1., 1.,
        x -> ghost_wall(x, [0, -1.]), # top
        x -> ghost_wall(x, [-1., 0.]), # right
        x -> ghost_wall(x, [0, 1.]), # bottom
        x -> ghost_wall(x, [1.0, 0.])) # left

    @test U ≈ MacCormack_step(U, ps)

    for i in 1:100
        U = MacCormack_step(U, ps)
        @test U ≈ U_original
    end
end

function test_stationary_T_discont()
    # The fluid has uniform pressure and uniformly zero
    # velocity. There is a temperature discontinuity. No
    # changes in the fluid state should occur.
    U = zeros(8, 8, 4)
    for i in 1:8
        for j in 1:8
            if i <= 4
                U[i, j, :] = pTvel2u(101e3, 300, 0, 0, air)
            else
                U[i, j, :] = pTvel2u(101e3, 200, 0, 0, air)
            end
        end
    end
    U_original = copy(U)

    ps = ProblemSpec(air, 1., 1., 1.,
        x -> x,
        x -> x,
        x -> x,
        x -> x)

    @test U ≈ MacCormack_step(U, ps)

    for i in 1:100
        U = MacCormack_step(U, ps)
        @test U ≈ U_original
    end
end


test_unif_open()
test_unif_wall()
test_stationary_T_discont()
println("Passed all tests.")
