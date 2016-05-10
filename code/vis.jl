# Visulaization

# Matt Vernacchia
# 18.086 Project
# Spring 2016

using PyPlot
using CFD086

export plot_U

function plot_U(U, ps::ProblemSpec)
    titles = ["Density", "x Momentum", "y Momentum", "Energy"]

    X = zeros(size(U,1), size(U,2))
    Y = zeros(size(U,1), size(U,2))

    for i in 1:size(U,1)
        for j in 1:size(U,2)
            X[i, j] = ps.x(i * ps.Δx, j * ps.Δy)
            Y[i, j] = ps.y(i * ps.Δx, j * ps.Δy)
        end
    end

    for i in 1:4
        subplot(2,2,i)
        pcolormesh(X', Y', U[:,:,i]')
        title(titles[i])
        xlabel("x [meter]")
        ylabel("y [meter]")
        colorbar()
    end
end