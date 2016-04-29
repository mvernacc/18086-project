"""
Main script.

Define a problem and run the solver.

Matt Vernacchia
18.086 Project
Spring 2016
"""

# I need to have a line in between the top docstring and include
# because of a language bug.
# See https://github.com/JuliaArchive/Markdown.jl/issues/25
unused = 0

include("bound.jl")
include("fluid_eqns.jl")
include("problem_spec.jl")
include("solver.jl")
include("stability.jl")

# Gas properties
gas = Gas(1.40, 1005)

# Inlet conditions
T = 300.
p  = 101e3
u = 10.
v = 0.
U_inlet = pTvel2u(p, T, u, v, gas)

# Boundary conditions.
# Top: solid wall.
function top_bound(U)
    return ghost_wall(U, [0., -1.])
end
# Right: blank outlet.
function right_bound(U)
    return U
end
# Bottom: solid wall.
function bottom_bound(U)
    return ghost_wall(U, [0., 1.])
end
# Left: Full-state inlet.
function left_bound(U)
    return U_inlet
end

# Step sizes
Δx = 1e-2
Δy = 1e-2
Δt = 0.5 * Δt_cfl(u, v, 340, Δx, Δy)

# Grid size
Nx = 10
Ny = 10

ps = ProblemSpec(gas, Δt, Δx, Δy, top_bound, right_bound, bottom_bound,
    left_bound)

# initial conditions
U = zeros(Nx, Ny, 4)
for i in 1:Nx
    for j in 1:Ny
        U[i, j, :] = U_inlet
    end
end

dump(U)

for it in 1:10
    U = MacCormack_step(U, ps)
end

dump(U)