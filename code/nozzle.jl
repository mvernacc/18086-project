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
M_in = 0.3
U_inlet = poToM2u(p_c, T_c, M_in, 0, air)

# Exit pressue [units: pascal]
p_e = 5e6 / 35
# Exit mach number [units: none]
# guess using fig 1-12 of Huzel and Huang
M_e = 3.0
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
    return U
    # return ghost_p(U, p_e, ps.gas)
end
# Bottom: solid wall.
function bottom_bound(U, i, j, ps)
    n = [-ps.dy_dξ(i * ps.Δx, j * ps.Δy), 1]
    n = n / norm(n)
    return ghost_wall(U, n)
end
# Left: Pressure and temperature inlet.
function left_bound(U, i, j, ps)
    # return ghost_pT(U, p_c, T_c, ps.gas)
    # return U_inlet
    # n = [1, ps.dy_dξ(0, j * ps.Δy)]
    # n = n / norm(n)
    # return pTM2u(
    #     p_c,
    #     T_c,
    #     M_in * n[1],
    #     M_in * n[2],
    #     air)
    return ghost_poTo(U, p_c, T_c, ps.gas)
end

# Grid size
Nx = 100
Ny = 20

# Nozzle shape parameters
# Chamber straight section length [units: meter]
x_s = 0.5
# Chamber radius [units: meter]
r_c = 0.5
# Throat radius [units: meter]
r_t = 0.25
# Exit radius [units: meter]
r_e = 1
# Convergent angle [units: radian]
θ_1 = deg2rad(15)
# Divergent circle angle [units: radian]
θ_2 = deg2rad(15)
# Find the nozzle shape parameters
nozzle_shape_param = nozzle_parameters(x_s, r_c, r_t, θ_1, θ_2, r_e)
# Nozzle length [units: meter]
x_e = nozzle_shape_param[11]
# Nozzle throat position [units: meter]
x_t = nozzle_shape_param[9]

# Step sizes
Δx = x_e / Nx
Δy = 1 / Ny
Δt = 0.4 * Δt_cfl(2e3, 3e2, u2a(U_inlet, air), Δx, Δy)

println(Δt)

# Grid shape
# function x(ξ, η)
#     return -3 / (4 * x_t^2) * ξ^3 +
#         3 / (4 * x_t) * ξ^2 +
#         ξ
# end

# function dx_dξ(ξ, η)
#     return -9 / (4 * x_t^2) * ξ^2
#         6 / (4 * x_t) * ξ +
#         1
# end

function y(ξ, η)
    return (η - Δy) * nozzle_contour(ξ - Δx, nozzle_shape_param)
end

function dy_dξ(ξ, η)
    return (η - Δy) * nozzle_contour_derivative(ξ - Δx, nozzle_shape_param)
end

function dy_dη(ξ, η)
    return nozzle_contour(ξ - Δx, nozzle_shape_param)
end

ps = ProblemSpec(air, Δt, Δx, Δy, top_bound, right_bound, bottom_bound,
    left_bound,
    (ξ, η) -> ξ, # x(ξ, η)
    y,
    (ξ, η) -> 1, # dx_dξ
    (ξ, η) -> 0, # dx_dη
    dy_dξ, # dy_dξ
    dy_dη, # dy_dη
    # (ξ, η) -> 1
    )

# initial conditions
U = zeros(Nx, Ny, 4)
for i in 1:Nx
    for j in 1:Ny
        n = [1, ps.dy_dξ(i * ps.Δx, j * ps.Δy)]
        n = n / norm(n)
        U[i, j, :] = pTM2u(
            p_c + (p_e - p_c) * i / Nx,
            T_c + (T_e - T_c) * i / Nx,
            M_in * n[1],
            M_in * n[2],
            air)
    end
end

dump(U)

figure(figsize=(16,8))
# plot_U(U, ps)
# savefig("results/nozzle/t=00000.000us.svg")

tic()
for it in 1:10
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