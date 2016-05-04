# Visulaization

# Matt Vernacchia
# 18.086 Project
# Spring 2016

using PyPlot

function plot_U(U)
    titles = ["Density", "x Momentum", "y Momentum", "Energy"]

    for i in 1:4
        subplot(2,2,i)
        pcolor(U[:,:,i]')
        title(titles[i])
        xlabel("x [grid cells]")
        ylabel("y [grid cells]")
        colorbar()
    end
end