"""
Shock tube to demonstrate effect of arificial diffusion.

Matt Vernacchia
18.086 Project
Spring 2016
"""

# I need to have a line in between the top docstring and include
# because of a language bug.
# See https://github.com/JuliaArchive/Markdown.jl/issues/25
unused = 0

include("cfd086.jl")
using CFD086
using PyPlot

# Inlet conditions
T = 300.
p  = 101e3
Mx = 0
My = 0
U_inlet = pTM2u(p, T, Mx, My, air)

# Boundary conditions.
# Top: solid wall.
function top_bound(U, i, j, ps)
    n = [ps.dy_dξ(i * ps.Δx, j * ps.Δy), -1]
    n = n / norm(n)
    return ghost_wall(U, n)
end
# Right: blank outlet.
function right_bound(U, i, j, ps)
    return U
end
# Bottom: solid wall.
function bottom_bound(U, i, j, ps)
    n = [-ps.dy_dξ(i * ps.Δx, j * ps.Δy), 1]
    n = n / norm(n)
    return ghost_wall(U, n)
end
# Left: Full-state inlet.
function left_bound(U, i, j, ps)
    return U_inlet
end

# Grid size
Nx = 80
Ny = 10

# Step sizes
Δx = 0.5 / Nx
Δy = 0.125 / Ny
Δt = 0.4 * Δt_cfl(U_inlet[2] / U_inlet[1], U_inlet[3] / U_inlet[1], 340, Δx, Δy)

println(Δt)


ps = ProblemSpec(air, Δt, Δx, Δy, top_bound, right_bound, bottom_bound,
    left_bound)

# initial conditions
U = zeros(Nx, Ny, 4)
for i in 1:Nx
    for j in 1:Ny
        if i < Nx/2
            U[i, j, :] = U_inlet
        else
            U[i, j, :] = pTM2u(500e3, T, 0, 0, air)
        end
    end
end

dump(U)
U_ad = copy(U)

function plot_p(U, U_ad)
     p = zeros(Nx)
    p_ad = zeros(Nx)
    for i in 1:Nx
        p[i] = u2p(U[i, 5, :], ps.gas)
        p_ad[i] = u2p(U_ad[i, 5, :], ps.gas)
    end
    plot(p * 1e-3, color="blue", label="Without AD")
    plot(p_ad * 1e-3, color="red", label="With AD")
    xlabel("i")
    ylabel("pressure [kPa]")
    ylim([0, 700])
    legend()
end

figure(figsize=(16,8))
plot_p(U, U_ad)
savefig("results/shock_tube/it0.svg")

tic()
for it in 1:20
    U = MacCormack_step(U, ps, use_ad=false)
    U_ad = MacCormack_step(U_ad, ps, use_ad=true)
       
    if it % 2 == 0
        clf()
        plot_p(U, U_ad)
        suptitle(@sprintf("Time step %d", it))
        savefig(@sprintf("results/shock_tube/it%d.svg", it))
    end
end
toc()

dump(U)
show()