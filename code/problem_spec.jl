using CFD086

export ProblemSpec

type ProblemSpec
    # The gas properties
    gas::Gas

    # Step sizes
    Δt::Float64
    # x direction step. On a curved grid, this is actually Δξ.
    Δx::Float64
    # x direction step. On a curved grid, this is actually Δη.
    Δy::Float64

    # Boundary condition definition:
    # Functions which map edge cell U's to ghost U's
    top_bound::Function
    right_bound::Function
    bottom_bound::Function
    left_bound::Function

    # Coordinate transform for curvilinear grid.
    # All are functions of (ξ, η).
    x::Function
    y::Function
    dx_dξ::Function
    dx_dη::Function
    dy_dξ::Function
    dy_dη::Function

    # Constructor for cartesian grid.
    ProblemSpec(gas, Δt, Δx, Δy,
        top_bound, right_bound, bottom_bound, left_bound) =
        new(gas, Δt, Δx, Δy,
            top_bound, right_bound, bottom_bound, left_bound,
            (ξ, η) -> ξ, # x(ξ, η)
            (ξ, η) -> η, # y(ξ, η)
            (ξ, η) -> 1, # dx_dξ
            (ξ, η) -> 0, # dx_dη
            (ξ, η) -> 0, # dy_dξ
            (ξ, η) -> 1, # dy_dη
        )
    # General constructor
    ProblemSpec(gas, Δt, Δx, Δy,
        top_bound, right_bound, bottom_bound, left_bound,
        x, y, dx_dξ, dx_dη, dy_dξ, dy_dη) =
        new(gas, Δt, Δx, Δy,
            top_bound, right_bound, bottom_bound, left_bound,
            x, y, dx_dξ, dx_dη, dy_dξ, dy_dη)
end