"""
Boundary conditions.

Matt Vernacchia
18.086 Project
Spring 2016
"""

# I need to have a line in between the top docstring and include
# because of a language bug.
# See https://github.com/JuliaArchive/Markdown.jl/issues/25
unused = 0

using CFD086

export ghost_zerograd, ghost_p, ghost_pT, ghost_wall, pad_bounds

# Ghost cell functions.
# These functions return the ghost cells which are used to
# take the derivatives at the boundaries.

@doc """
Zero-gradient boundary. Useful for outflow.

Arguments:
    U: 4-vector of conservative variables, at the cell on the edge of the grid.

Returns:
    4-vector of conservative variables, at the ghost cell on the outside of the boundary.
""" ->
function ghost_zerograd(U::Array{Float64, 1})
    return U
end


@doc """
Pressure boundary.

Arguments:
    U: 4-vector of conservative variables, at the cell on the edge of the grid.
    p_bound: boundary static pressure [units: pascal].
    gas: Gas properties.

Returns:
    4-vector of conservative variables, at the ghost cell on the outside of the boundary.
""" ->
function ghost_p(U::Array{Float64, 1}, p_bound::Float64, gas::Gas)
    e = u2e(U)
    ρ_bound = p_bound / (gas.γ - 1) / e
    u = U[2] / U[1]
    v = U[3] / U[1]
    return [ρ_bound, ρ_bound * u, ρ_bound * v, ρ_bound * (e + (u^2 + v^2) / 2)]
end


@doc """
Pressure and temperature boundary.

Arguments:
    U: 4-vector of conservative variables, at the cell on the edge of the grid.
    p_bound: boundary static pressure [units: pascal].
    T_bound: boundary static temperature [units: kelvin].
    gas: Gas properties.

Returns:
    4-vector of conservative variables, at the ghost cell on the outside of the boundary.
""" ->
function ghost_pT(U::Array{Float64, 1}, p_bound::Float64, T_bound::Float64, gas::Gas)
    u = U[2] / U[1]
    v = U[3] / U[1]
    return pTvel2u(p_bound, T_bound, u, v, gas)
end


@doc """
Solid wall boundary.
Arguments:
    U: 4-vector of conservative variables, at the cell on the edge of the grid.
    n: 2-vector normal to wall, unit norm.

Returns:
    4-vector of conservative variables, at the ghost cell on the outside of the boundary.
""" ->
function ghost_wall(U::Array{Float64, 1}, n::Array{Float64, 1})
    # Reflect momentum vector across the wall.
    # Source: http://mathworld.wolfram.com/Reflection.html
    mom_ref = U[2:3] - 2 * (U[2:3] ⋅ n) * n
    return [U[1], mom_ref[1], mom_ref[2], U[4]]
end


@doc """
Pad a grid with boundary ghost cells.
"""
function pad_bounds(U::Array{Float64,3}, ps::ProblemSpec;
    level=1)

    if level > 1
        U = pad_bounds(U, ps,
                level=level - 1)
    end

    U_pad = zeros(size(U, 1) + 2, size(U, 2) + 2, 4)
    for i in level:size(U, 1) - level + 1
        U_pad[i+1, end, :] = ps.top_bound(squeeze(U[i, end - 2 * level + 2, :], (1,2)),
            i, size(U,2) - 2 * level + 2, ps)
    end
    for j in level:size(U, 2) - level + 1
        U_pad[end, j+1, :] = ps.right_bound(squeeze(U[end - 2 * level + 2, j, :], (1,2)),
            size(U,1) - 2 * level + 2, j, ps)
    end
    for i in level:size(U, 1)  - level + 1
        U_pad[i+1, 1, :] = ps.bottom_bound(squeeze(U[i, 2 * level - 1, :], (1,2)),
            i, 2 * level - 1, ps)
    end
    for j in level:size(U, 2) - level + 1
        U_pad[1, j+1, :] = ps.left_bound(squeeze(U[2 * level - 1, j, :], (1,2)),
            2 * level - 1, j, ps)
    end
    U_pad[2:end-1, 2:end-1, :] = U
    return U_pad
end