# Unit tests for bound.jl
#
# Matt Vernacchia
# 18.086 Project
# Spring 2016

using Base.Test

include("cfd086.jl")
using CFD086

ps = ProblemSpec(air, 1, 1, 1,
    (x, i, j, ps) -> x + 0.1,
    (x, i, j, ps) -> x + 0.2,
    (x, i, j, ps) -> x + 0.3,
    (x, i, j, ps) -> x + 0.4,
    )

function test_pad_bounds()
    A = [1. 2. 3.;
         4. 5. 6.;
         7. 8. 9.]
    U = cat(3, A, A, A, A)

    U_pad = pad_bounds(U, ps) 

    # Dimensions
    @test size(U_pad) == (5,5,4)

    # Top bound
    @test U_pad[2, end, 1] == 3.1
    @test U_pad[3, end, 1] == 6.1
    @test U_pad[4, end, 1] == 9.1

    # Right bound
    @test U_pad[end, 2, 1] == 7.2
    @test U_pad[end, 3, 1] == 8.2
    @test U_pad[end, 4, 1] == 9.2

    # Bottom bound
    @test U_pad[2, 1, 1] == 1.3
    @test U_pad[3, 1, 1] == 4.3
    @test U_pad[4, 1, 1] == 7.3

    # Left bound
    @test U_pad[1, 2, 1] == 1.4
    @test U_pad[1, 3, 1] == 2.4
    @test U_pad[1, 4, 1] == 3.4
end

function test_pad_bounds_level2()
    A = [1. 2. 3.;
         4. 5. 6.;
         7. 8. 9.]
    U = cat(3, A, A, A, A)

    U_pad = pad_bounds(U, ps,
        level=2)

    dump(U_pad)

    # Dimensions
    @test size(U_pad) == (7,7,4)

    # Top bound
    @test U_pad[3, end-1, 1] == 3.1
    @test U_pad[4, end-1, 1] == 6.1
    @test U_pad[5, end-1, 1] == 9.1
    @test U_pad[3, end, 1] == 2.1
    @test U_pad[4, end, 1] == 5.1
    @test U_pad[5, end, 1] == 8.1

    # Right bound
    @test U_pad[end-1, 3, 1] == 7.2
    @test U_pad[end-1, 4, 1] == 8.2
    @test U_pad[end-1, 5, 1] == 9.2
    @test U_pad[end, 3, 1] == 4.2
    @test U_pad[end, 4, 1] == 5.2
    @test U_pad[end, 5, 1] == 6.2

    # Bottom bound
    @test U_pad[3, 2, 1] == 1.3
    @test U_pad[4, 2, 1] == 4.3
    @test U_pad[5, 2, 1] == 7.3
    @test U_pad[3, 1, 1] == 2.3
    @test U_pad[4, 1, 1] == 5.3
    @test U_pad[5, 1, 1] == 8.3

    # Left bound
    @test U_pad[2, 3, 1] == 1.4
    @test U_pad[2, 4, 1] == 2.4
    @test U_pad[2, 5, 1] == 3.4
    @test U_pad[1, 3, 1] == 4.4
    @test U_pad[1, 4, 1] == 5.4
    @test U_pad[1, 5, 1] == 6.4
end


function test_ghost_wall()
    @test ghost_wall([0, -1., -1., 0], [1., 0.]) ≈ [0, 1., -1., 0]
    @test ghost_wall([0, 2., -1., 0], [0., 1.]) ≈ [0, 2., 1., 0]
    @test ghost_wall([0, 100., 0., 0], [1/2^0.5, 1/2^0.5]) ≈ [0, 0., -100., 0]
end


test_pad_bounds()
test_pad_bounds_level2()
test_ghost_wall()
println("Passed all tests.")