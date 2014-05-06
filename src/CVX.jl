module CVX

# package code goes here
include("expressions/expressions.jl")
include("constraints/constraints.jl")
include("problems/problems.jl")
include("utilities/utilities.jl")
include("solvers/ecos.jl")
include("atoms/add_subtract.jl")
include("atoms/mul_div.jl")
include("atoms/affine.jl")
include("atoms/minmax.jl")
include("atoms/norm.jl")
include("atoms/index.jl")


end # module
