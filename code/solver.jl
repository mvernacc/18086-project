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

using CFD086

export MacCormack_step

function MacCormack_step(U::Array{Float64,3},
    ps::ProblemSpec;
    use_ad::Bool=false)

    # Pad U with boundary ghost cells.
    U_pad = pad_bounds(U, ps)
    domain_i = 2:size(U,1)+1
    domain_j = 2:size(U,2)+1

    min_ρ = 1e-2
    
    F = (U, i, j) -> ps.dy_dη(i * ps.Δx, j * ps.Δy) * F_euler(U[i, j, :], ps.gas) - 
        ps.dx_dη(i * ps.Δx, j * ps.Δy) * G_euler(U[i, j, :], ps.gas)
    G = (U, i, j) ->  - ps.dy_dξ(i * ps.Δx, j * ps.Δy) * F_euler(U[i, j, :], ps.gas) +
        ps.dx_dξ(i * ps.Δx, j * ps.Δy) * G_euler(U[i, j, :], ps.gas)
    if use_ad
        F = (U, i, j) -> ps.dy_dη(i * ps.Δx, j * ps.Δy) * F_euler(U[i, j, :], ps.gas) - 
            ps.dx_dη(i * ps.Δx, j * ps.Δy) * G_euler(U[i, j, :], ps.gas) -
            ad_F(U, i, j, ps.gas) / J[i, j]
        G = (U, i, j) -> - ps.dy_dξ(i * ps.Δx, j * ps.Δy) * F_euler(U[i, j, :], ps.gas) +
            ps.dx_dξ(i * ps.Δx, j * ps.Δy) * G_euler(U[i, j, :], ps.gas) -
            ad_G(U, i, j, ps.gas) / J[i, j]
        U_pad = pad_bounds(U, ps,
            level=4)
        domain_i = 5:size(U,1) + 4
        domain_j = 5:size(U,2) + 4
    end

    # Apply coordinate transform correction
    J = ct_J(size(U_pad, 1), size(U_pad, 2), ps,
        i_off=(size(U_pad, 1) -  size(U, 1)) / 2,
        j_off=(size(U_pad, 2) -  size(U, 2)) / 2)
    U_pad = J .* U_pad

    U_p = zeros(size(U_pad))
    U_c = zeros(size(U_pad))

    # Predictor
    for i in domain_i
        for j in domain_j
            U_p[i, j, :] = squeeze(U_pad[i, j, :], (1,2)) -
                ps.Δt / (ps.Δx * J[i, j]) * (F(U_pad, i+1, j)
                    - F(U_pad, i, j)) -
                ps.Δt / (ps.Δy * J[i, j]) * (G(U_pad, i, j+1)
                    - G(U_pad, i, j))

            # Correct negative density cells
            if U_p[i, j, 1] <= 0
                U_p[i, j, 2] *= min_ρ / U_p[i, j, 1]
                U_p[i, j, 3] *= min_ρ / U_p[i, j, 1]
                U_p[i, j, 1] = min_ρ
            end
        end
    end

    # Remove padding zeros from U_p and pad U_p with boundary condition ghost cells
    if use_ad
        U_p = U_p[5:end-4, 5:end-4, :]
        U_p = pad_bounds(U_p, ps,
            level=4)
    else
        U_p = U_p[2:end-1, 2:end-1, :]
        U_p = pad_bounds(U_p, ps)
    end

    # Corrector
    for i in domain_i
        for j in domain_j
            U_c[i, j, :] = 0.5 * (squeeze(U_pad[i, j, :] + U_p[i, j, :], (1,2)) -
                ps.Δt / (ps.Δx * J[i, j]) * (F(U_p, i, j)
                    - F(U_p, i-1, j)) -
                ps.Δt / (ps.Δy * J[i, j]) * (G(U_p, i, j)
                    - G(U_p, i, j-1)))

            # Correct negative density cells
            if U_c[i, j, 1] <= 0
                U_c[i, j, 2] *= min_ρ / U_c[i, j, 1]
                U_c[i, j, 3] *= min_ρ / U_c[i, j, 1]
                U_c[i, j, 1] = min_ρ
            end
        end
    end

    # Un-pad
    if use_ad
        U_c = U_c[5:end-4, 5:end-4, :]
    else
        U_c = U_c[2:end-1, 2:end-1, :]
    end

    return U_c
end


# Artificial Diffusion
@doc """
ϕ function for artificial diffusion, x direction.

Eqn 9 in Lobao.
""" ->
function ad_ϕ_x(U, i::Integer, j::Integer)
    u2 = zeros(4)
    if i > 2
        u2 = U[i-2, j]
    end
    return U[i+1, j] - 3 * U[i, j] + 3 * U[i-1, j] - u2
end


@doc """
ϕ function for artificial diffusion, y direction.

Eqn 9 in Lobao.
""" ->
function ad_ϕ_y(U, i::Integer, j::Integer)
    u2 = zeros(4)
    if j > 2
        u2 = U[i, j-2]
    end
    return U[i, j+1] - 3 * U[i, j] + 3 * U[i, j-1] - u2
end


@doc """
ν function for artificial diffusion, x direction.

Eqn 10 in Lobao
""" ->
function ad_ν_x(U, i::Integer, j::Integer, gas::Gas)
    return abs(u2p(U[i+1, j, :], gas) - 2 * u2p(U[i, j, :], gas) + u2p(U[i-1, j, :], gas)) /
        (u2p(U[i+1, j, :], gas) + 2 * u2p(U[i, j, :], gas) + u2p(U[i-1, j, :], gas))
end


@doc """
ν function for artificial diffusion, y direction.

Eqn 10 in Lobao
""" ->
function ad_ν_y(U, i::Integer, j::Integer, gas::Gas)
    return abs(u2p(U[i, j+1, :], gas) - 2 * u2p(U[i, j, :], gas) + u2p(U[i, j-1, :], gas)) /
        (u2p(U[i, j+1, :], gas) + 2 * u2p(U[i, j, :], gas) + u2p(U[i, j-1, :], gas))
end


# Artificial diffusion constants
const k2 = 1 / 4
const k4 = 1 / 64


@doc """
ε2 function for artificial diffusion, x direction.

Eqn 11 in Lobao
""" ->
function ad_ε2_x(U, i::Integer, j::Integer, gas::Gas)
    return k2 * ad_ν_x(U, i-1, j, gas) * (norm(u2vel(U[i,j,:])) + u2a(U[i,j,:], gas))
end


@doc """
ε2 function for artificial diffusion, y direction.

Eqn 11 in Lobao
""" ->
function ad_ε2_y(U, i::Integer, j::Integer, gas::Gas)
    return k2 * ad_ν_y(U, i, j-1, gas) * (norm(u2vel(U[i,j,:])) + u2a(U[i,j,:], gas))
end


@doc """
ε4 function for artificial diffusion, x direction.

Eqn 11 in Lobao
""" ->
function ad_ε4_x(U, i::Integer, j::Integer, gas::Gas)
    return max(0, k4 - ad_ε2_x(U, i, j, gas)) * (norm(u2vel(U[i,j,:])) + u2a(U[i,j,:], gas))
end


@doc """
ε4 function for artificial diffusion, y direction.

Eqn 11 in Lobao
""" ->
function ad_ε4_y(U, i::Integer, j::Integer, gas::Gas)
    return max(0, k4 - ad_ε2_y(U, i, j, gas)) * (norm(u2vel(U[i,j,:])) + u2a(U[i,j,:], gas))
end


@doc """
d function for artificial diffusion, x direction.

Eqn 8 in Lobao
""" ->
function ad_d_x(U, i::Integer, j::Integer, gas::Gas)
    return ad_ε2_x(U, i, j, gas) * (U[i,j,:] - U[i-1,j,:]) -
        ad_ε4_x(U, i, j, gas) * ad_ϕ_x(U,i,j)
end


@doc """
d function for artificial diffusion, y direction.

Eqn 8 in Lobao
""" ->
function ad_d_y(U, i::Integer, j::Integer, gas::Gas)
    return ad_ε2_y(U, i, j, gas) * (U[i,j,:] - U[i,j-1,:]) -
        ad_ε4_y(U, i, j, gas) * ad_ϕ_y(U, i, j)
end


@doc """
Artificial diffusion in the x direction.
""" ->
function ad_F(U, i::Integer, j::Integer, gas::Gas)
    return 0.5 * squeeze(ad_d_x(U, i+1, j, gas) - ad_d_x(U, i-1, j, gas), (1,2))
end


@doc """
Artificial diffusion in the y direction.
""" ->
function ad_G(U, i::Integer, j::Integer, gas::Gas)
    return 0.5 * squeeze(ad_d_y(U, i, j+1, gas) - ad_d_y(U, i, j-1, gas), (1,2))
end


@doc """
Coordinate transform correction factor, J
""" ->
function ct_J(nx, ny, ps::ProblemSpec; i_off=0, j_off=0)
    J = zeros(nx, ny)
    for i in 1:nx
        for j in 1:ny
            ξ = (i - i_off) * ps.Δx
            η = (j - j_off) * ps.Δy
            J[i, j] = ps.dx_dξ(ξ, η) * ps.dy_dη(ξ, η) - ps.dx_dη(ξ, η) * ps.dy_dξ(ξ, η)
        end
    end
    return J
end