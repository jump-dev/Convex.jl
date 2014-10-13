module Convex

### modeling framework
include("dcp.jl")
include("expressions.jl")
include("variable.jl")
include("constant.jl")
include("conic_form.jl")
include("constraints/constraints.jl")
include("constraints/exponential_constraints.jl")
include("constraints/sdp_constraints.jl")
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
include("atoms/affine/reshape.jl")

### elementwise atoms
include("atoms/abs.jl")
include("atoms/maximum.jl")
include("atoms/minimum.jl")
include("atoms/max.jl")
include("atoms/min.jl")
include("atoms/norm.jl")

### SOC atoms
include("atoms/second_order_cone/norm_2.jl")
include("atoms/second_order_cone/quad_over_lin.jl")
include("atoms/second_order_cone/qol_elementwise.jl")
include("atoms/second_order_cone/geo_mean.jl")

### SDP atoms
include("atoms/matrix_norm.jl")

### exponential atoms
include("atoms/exp_cone/exp.jl")
include("atoms/exp_cone/log.jl")
include("atoms/exp_cone/logsumexp.jl")

### utilities
include("utilities/show.jl")
end
