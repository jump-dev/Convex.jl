module CVX

# package code goes here
include("union.jl")
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
include("atoms/affine/reshape.jl")
include("atoms/affine/transpose.jl")
include("atoms/affine/diag.jl")

include("atoms/elementwise/min.jl")
include("atoms/elementwise/max.jl")
include("atoms/elementwise/pos.jl")
include("atoms/elementwise/neg.jl")
include("atoms/elementwise/abs.jl")
include("atoms/elementwise/sqrt.jl")
include("atoms/elementwise/square.jl")
include("atoms/elementwise/square_pos.jl")
include("atoms/elementwise/qol_elementwise.jl")
include("atoms/elementwise/inv_pos.jl")

include("atoms/norm.jl")
include("atoms/quad_form.jl")
include("atoms/geo_mean.jl")
include("atoms/sum_squares.jl")
include("atoms/quad_over_lin.jl")

end # module
