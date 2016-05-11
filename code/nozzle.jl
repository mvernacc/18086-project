"""
Nozzle flow.

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
T_c = 2000.
p_c  = 5e6
Mx = 0.3
My = 0
U_inlet = pTM2u(p_c, T_c, Mx, My, air)

# Exit pressue [units: pascal]
p_e = 101e3
# Exit mach number [units: none]
# guess using fig 1-12 of Huzel and Huang
M_e = 3.5
# Exit temperature [units: Kelvin]
T_e = T_c * (p_e / p_c)^((air.γ - 1) / air.γ)

# Boundary conditions.
# Top: solid wall.
function top_bound(U, i, j, ps)
    n = [ps.dy_dξ(i * ps.Δx, j * ps.Δy), -1]
    n = n / norm(n)
    return ghost_wall(U, n)
end
# Right: pressure.
function right_bound(U, i, j, ps)
    return ghost_p(U, p_e, ps.gas)
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
Nx = 120
Ny = 40

# Nozzle shape parameters
# Chamber radius [units: meter]
r_c = 0.5
# Throat radius [units: meter]
r_t = 0.25
# Exit radius [units: meter]
r_e = 2
# Convergent angle [units: radian]
θ_1 = deg2rad(30)
# Divergent circle angle [units: radian]
θ_2 = deg2rad(25)
# Find the nozzle shape parameters
nozzle_shape_param = nozzle_parameters(r_c, r_t, θ_1, θ_2, r_e)
# Nozzle length [units: meter]
x_e = nozzle_shape_param[10]

# Step sizes
Δx = x_e / Nx
Δy = 1 / Ny
Δt = 0.4 * Δt_cfl(U_inlet[2] / U_inlet[1], U_inlet[3] / U_inlet[1], 340, Δx, Δy)

println(Δt)

# Grid shape
function y(ξ, η)
    return (η - Δy) * nozzle_contour(ξ, nozzle_shape_param)
end

function dy_dξ(ξ, η)
    return (η - Δy) * nozzle_contour_derivative(ξ, nozzle_shape_param)
end

function dy_dη(ξ, η)
    return nozzle_contour(ξ, nozzle_shape_param)
end

ps = ProblemSpec(air, Δt, Δx, Δy, top_bound, right_bound, bottom_bound,
    left_bound,
    (ξ, η) -> ξ, # x(ξ, η)
    y,
    (ξ, η) -> 1, # dx_dξ
    (ξ, η) -> 0, # dx_dη
    dy_dξ, # dy_dξ
    (ξ, η) -> 1, # dy_dη
    )

# initial conditions
U = zeros(Nx, Ny, 4)
for i in 1:Nx
    for j in 1:Ny
        U[i, j, :] = pTM2u(p_c + (p_e - p_c) * i / Nx,
            T_c + (T_e - T_c) * i / Nx,
            Mx,
            My,
            air)
    end
end

dump(U)

figure(figsize=(16,8))
# plot_U(U, ps)
# savefig("results/nozzle/t=00000.000us.svg")

tic()
for it in 1:100
    U = MacCormack_step(U, ps, use_ad=true)
    # if it % 10 == 0
    #     clf()
    #     plot_U(U, ps)
    #     suptitle(@sprintf("Nozzle, t=%5.3fus.svg", it * Δt * 1e6))
    #     savefig(@sprintf("results/nozzle/t=%5.3fus.svg", it * Δt * 1e6))
    # end
end
toc()

dump(U)
plot_U(U, ps)
figure(figsize=(16,8))
plot_pTM(U, ps)
show()