using CFD086

export ProblemSpec

type ProblemSpec
    # The gas properties
    gas::Gas

    # Step sizes
    Δt::Float64
    Δx::Float64
    Δy::Float64

    # Boundary condition definition:
    # Functions which map edge cell U's to ghost U's
    top_bound::Function
    right_bound::Function
    bottom_bound::Function
    left_bound::Function
end