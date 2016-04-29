"""
Stability checks.

Matt Vernacchia
18.086 Project
Spring 2016
"""
# I need to have a line in between the top docstring and include
# because of a language bug.
# See https://github.com/JuliaArchive/Markdown.jl/issues/25
unused = 0


@doc """
Get the time step CFL limit.

Arguments:
    u: maximum expected x velocity.
    v: maximum expected y velocity.
    a: maximum expected speed of sound.
    Δx: x spatial step size.
    Δy: y spatial step size.
""" ->
function Δt_cfl(u, v, a, Δx, Δy)
    return (u / Δx + v / Δy + a * (
        1 / Δx^2 + 1 / Δy^2)^0.5 )^-1
end
