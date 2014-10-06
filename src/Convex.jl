module Convex

### modeling framework
include("dcp.jl")
include("expressions.jl")
include("variable.jl")
include("constant.jl")
include("conic_form.jl")
include("constraints.jl")
include("problems.jl")
include("solution.jl")

### affine atoms
include("atoms/affine/add_subtract.jl")
include("atoms/affine/multiply_divide.jl")
include("atoms/affine/sum.jl")
include("atoms/affine/transpose.jl")
include("atoms/affine/index.jl")
include("atoms/affine/diag.jl")
include("atoms/affine/stack.jl")
include("atoms/affine/dot.jl")

### elementwise atoms
include("atoms/abs.jl")
include("atoms/maximum.jl")
include("atoms/minimum.jl")
include("atoms/max.jl")
include("atoms/min.jl")

### SOC atoms
include("atoms/second_order_cone/norm2.jl")

### SDP atoms
include("atoms/matrix_norm.jl")

### exponential atoms

### utilities
include("utilities/show.jl")
end
