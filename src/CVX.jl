module CVX

# package code goes here
include("expressions/expressions.jl")
include("expressions/arithmetic.jl")
include("constraints/constraints.jl")
include("problems/problems.jl")
include("utilities/utilities.jl")
include("solvers/ecos.jl")
include("atoms/affine.jl")

end # module
