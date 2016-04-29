"""
Solvers.

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

function MacCormack_step(U::Array{Float64,3},
    ps::ProblemSpec)
    # Pad U with boundary ghost cells.
    U_pad = pad_bounds(U, ps.top_bound, ps.right_bound, ps.bottom_bound,
        ps.left_bound)
    U_p = zeros(size(U_pad))
    U_c = zeros(size(U_pad))

    # Predictor
    for i in 2:size(U,1)+1
        for j in 2:size(U,2)+1
            U_p[i, j, :] = squeeze(U_pad[i, j, :], (1,2)) -
                ps.Δt / ps.Δx * (F_euler(U_pad[i+1, j, :], ps.gas)
                    - F_euler(U_pad[i, j, :], ps.gas)) -
                ps.Δt / ps.Δy * (G_euler(U_pad[i, j+1, :], ps.gas)
                    - G_euler(U_pad[i, j, :], ps.gas))
        end
    end

    # Corrector
    for i in 2:size(U,1)+1
        for j in 2:size(U,2)+1
            U_c[i, j, :] = 0.5 * (squeeze(U_pad[i, j, :] + U_p[i, j, :], (1,2)) -
                ps.Δt / ps.Δx * (F_euler(U_pad[i, j, :], ps.gas)
                    - F_euler(U_pad[i-1, j, :], ps.gas)) -
                ps.Δt / ps.Δy * (G_euler(U_pad[i, j, :], ps.gas)
                    - G_euler(U_pad[i, j-1, :], ps.gas)))
        end
    end

    # Un-pad and return
    return U_c[2:end-1, 2:end-1, :]
end

