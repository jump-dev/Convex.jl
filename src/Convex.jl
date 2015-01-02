module Convex

global DEFAULT_SOLVER = nothing
### modeling framework
include("dcp.jl")
include("expressions.jl")
include("conic_form.jl")
include("variable.jl")
include("constant.jl")
include("constraints/constraints.jl")
include("constraints/signs_and_sets.jl")
include("constraints/soc_constraints.jl")
include("constraints/exp_constraints.jl")
include("constraints/sdp_constraints.jl")
include("problems.jl")
include("solver_info.jl")
include("solution.jl")

### affine atoms
include("atoms/affine/add_subtract.jl")
include("atoms/affine/multiply_divide.jl")
include("atoms/affine/sum.jl")
include("atoms/affine/transpose.jl")
include("atoms/affine/index.jl")
include("atoms/affine/diag.jl")
include("atoms/affine/antidiag.jl")
include("atoms/affine/stack.jl")
include("atoms/affine/dot.jl")
include("atoms/affine/reshape.jl")
include("atoms/affine/trace.jl")

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
include("atoms/second_order_cone/quad_form.jl")

### SDP atoms
include("atoms/sdp_cone/nuclear_norm.jl")
include("atoms/sdp_cone/operator_norm.jl")
include("atoms/sdp_cone/lambda_min_max.jl")

### exponential atoms
include("atoms/exp_cone/exp.jl")
include("atoms/exp_cone/log.jl")
include("atoms/exp_cone/logsumexp.jl")
include("atoms/exp_cone/entropy.jl")

### other elementwise atoms
include("atoms/huber.jl")
include("atoms/logdet.jl")

### utilities
include("utilities/show.jl")
end
