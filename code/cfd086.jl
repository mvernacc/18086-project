# Top-level module
#
# Matt Vernacchia
# 18.086 Project
# Spring 2016

module CFD086
include("gas.jl")
include("fluid_eqns.jl")
include("problem_spec.jl")
include("bound.jl")
include("stability.jl")
include("solver.jl")
include("vis.jl")
include("rao.jl")
end