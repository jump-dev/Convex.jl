module Convex
import DataStructures

importall Base.Operators

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
include("solver_info.jl")
include("problems.jl")
include("solution.jl")

### affine atoms
include("atoms/affine/add_subtract.jl")
include("atoms/affine/multiply_divide.jl")
include("atoms/affine/sum.jl")
include("atoms/affine/transpose.jl")
include("atoms/affine/index.jl")
include("atoms/affine/diag.jl")
include("atoms/affine/diagm.jl")
include("atoms/affine/stack.jl")
include("atoms/affine/dot.jl")
include("atoms/affine/reshape.jl")
include("atoms/affine/trace.jl")
include("atoms/affine/conv.jl")

### LP atoms
include("atoms/abs.jl")
include("atoms/maximum.jl")
include("atoms/minimum.jl")
include("atoms/max.jl")
include("atoms/min.jl")
include("atoms/sumlargest.jl")
include("atoms/dotsort.jl")

### some LP, some SOC atoms
include("atoms/norm.jl")

### SOC atoms
include("atoms/second_order_cone/norm2.jl")
include("atoms/second_order_cone/quadoverlin.jl")
include("atoms/second_order_cone/qol_elementwise.jl")
include("atoms/second_order_cone/geomean.jl")
include("atoms/second_order_cone/quadform.jl")
include("atoms/second_order_cone/power_to_socp.jl")
include("atoms/second_order_cone/rationalnorm.jl")
include("atoms/second_order_cone/huber.jl")

### SDP atoms
include("atoms/sdp_cone/nuclearnorm.jl")
include("atoms/sdp_cone/operatornorm.jl")
include("atoms/sdp_cone/lambda_min_max.jl")
include("atoms/sdp_cone/matrixfrac.jl")
include("atoms/sdp_cone/sumlargesteigs.jl")
include("atoms/sdp_cone/nthroot_det.jl")

### exponential atoms
include("atoms/exp_cone/exp.jl")
include("atoms/exp_cone/log.jl")
include("atoms/exp_cone/logsumexp.jl")
include("atoms/exp_cone/entropy.jl")
include("atoms/exp_cone/relative_entropy.jl")

### other atoms
include("atoms/logdet.jl")

### utilities
include("utilities/show.jl")
include("utilities/iteration.jl")
include("utilities/deprecated.jl")

#Temporary workaround for memory leak (https://github.com/JuliaOpt/Convex.jl/issues/83)
function clearmemory()
    global id_to_variables = Dict{UInt64, Variable}()
    global var_to_ranges = Dict{UInt64, Tuple{Int, Int}}()
    global conic_constr_to_constr = Dict{ConicConstr, Constraint}()
    gc()
end

end
