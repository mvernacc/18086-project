# Unit tests for bound.jl
#
# Matt Vernacchia
# 18.086 Project
# Spring 2016

using Base.Test

include("bound.jl")

function test_pad_bounds()
    A = [1. 2. 3.;
         4. 5. 6.;
         7. 8. 9.]
    U = cat(3, A, A, A, A)

    U_pad = pad_bounds(U,
        x -> x + 0.1,
        x -> x + 0.2,
        x -> x + 0.3,
        x -> x + 0.4) 

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

test_pad_bounds()