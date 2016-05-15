# Visulaization

# Matt Vernacchia
# 18.086 Project
# Spring 2016

using PyPlot
using CFD086

export plot_U, plot_pTM

function plot_U(U, ps::ProblemSpec)
    titles = ["Density [kg m^-3]", "x Momentum [kg m^-2 s^-1]", "y Momentum [kg m^-2 s^-1]", "Energy [Pa]"]

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
    tight_layout()
end

function plot_pTM(U, ps::ProblemSpec)
    titles = ["Pressure [kPa]", "Mach number [-]", "Temperature [K]"]

    X = zeros(size(U,1), size(U,2))
    Y = zeros(size(U,1), size(U,2))

    for i in 1:size(U,1)
        for j in 1:size(U,2)
            X[i, j] = ps.x(i * ps.Δx, j * ps.Δy)
            Y[i, j] = ps.y(i * ps.Δx, j * ps.Δy)
        end
    end

    result = zeros(size(U,1), size(U,2), 3)

     for i in 1:size(U,1)
        for j in 1:size(U,2)
            result[i, j, 1] = u2p(U[i,j,:], ps.gas) * 1e-3
            result[i, j, 2] = u2M(U[i,j,:], ps.gas)
            result[i, j, 3] = u2T(U[i,j,:], ps.gas)
        end
    end

    for i in 1:3
        subplot(2,2,i)
        pcolormesh(X', Y', result[:,:,i]')
        title(titles[i])
        xlabel("x [meter]")
        ylabel("y [meter]")
        colorbar()
    end
    tight_layout()
end