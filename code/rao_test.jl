# Test case for
# Rao nozzle geometry parabolic approximation
#
# Matt Vernacchia
# 18.086 Project
# Spring 2016

include("cfd086.jl")
using CFD086
using PyPlot

p = nozzle_parameters(0.5, 0.5, 0.25, deg2rad(15), deg2rad(15), 2)
dump(p)
println(p)
x = collect(0:1e-2:p[11])
y = [nozzle_contour(x1, p) for x1 in x]
dydx = [nozzle_contour_derivative(x1, p) for x1 in x]
dump(x)
dump(y)
subplot(2,1,1)
plot(x, y)
subplot(2,1,2)
plot(x, dydx, marker="x")
show()
