module CVX

# package code goes here
include("expressions/expressions.jl")
include("constraints/constraints.jl")
include("solution.jl")
include("problems/problems.jl")
include("utilities/utilities.jl")
include("solvers/ecos.jl")

# Atoms
include("atoms/affine/add_subtract.jl")
include("atoms/affine/dot.jl")
include("atoms/affine/index.jl")
include("atoms/affine/multiply_divide.jl")
include("atoms/affine/stack.jl")
include("atoms/affine/sum.jl")
include("atoms/affine/transpose.jl")
include("atoms/elementwise/minmax.jl")
include("atoms/norm.jl")

end # module
