module CVX

# package code goes here
include("expressions.jl")
include("arithmetic.jl")
include("constraints.jl")
include("problems.jl")
include("solvers/cvxpy.jl")
include("atoms.jl")

end # module
