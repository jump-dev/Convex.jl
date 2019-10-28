__precompile__()

module Convex
using OrderedCollections: OrderedDict
using LinearAlgebra
using SparseArrays

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

global DEFAULT_SOLVER = nothing
### modeling framework
include("dcp.jl")
include("expressions.jl")
# need to define `Variable` before `UniqueConicForms`
include("variable.jl")
include("conic_form.jl")
# need to define `conic_form!` for `Variable`s after `UniqueConicForms`
include("variable_conic_form.jl")
include("constant.jl")
include("constraints/constraints.jl")
include("constraints/signs_and_sets.jl")
include("constraints/soc_constraints.jl")
include("constraints/exp_constraints.jl")
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
include("atoms/affine/diagm.jl")
include("atoms/affine/stack.jl")
include("atoms/affine/dot.jl")
include("atoms/affine/reshape.jl")
include("atoms/affine/trace.jl")
include("atoms/affine/partialtrace.jl")
include("atoms/affine/partialtranspose.jl")
include("atoms/affine/conv.jl")
include("atoms/affine/real_imag.jl")
include("atoms/affine/inner_product.jl")
include("atoms/affine/conjugate.jl")
include("atoms/affine/kron.jl")

### lp atoms
include("atoms/lp_cone/abs.jl")
include("atoms/lp_cone/maximum.jl")
include("atoms/lp_cone/minimum.jl")
include("atoms/lp_cone/max.jl")
include("atoms/lp_cone/min.jl")
include("atoms/lp_cone/sumlargest.jl")
include("atoms/lp_cone/dotsort.jl")

### SOC atoms
include("atoms/second_order_cone/norm.jl") # also includes some lp atoms
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

### exponential atoms
include("atoms/exp_cone/exp.jl")
include("atoms/exp_cone/log.jl")
include("atoms/exp_cone/logsumexp.jl")
include("atoms/exp_cone/entropy.jl")
include("atoms/exp_cone/relative_entropy.jl")

### exp + sdp atoms
include("atoms/exp_+_sdp_cone/logdet.jl")

### utilities
include("utilities/tree_print.jl")
include("utilities/tree_interface.jl")
include("utilities/show.jl")
include("utilities/iteration.jl")
include("utilities/broadcast.jl")
include("problem_depot/problem_depot.jl")

# Deprecated workaround for memory leak (https://github.com/JuliaOpt/Convex.jl/issues/83)
function clearmemory()
    Base.depwarn("Convex.clearmemory() is deprecated, as the memory leak it works around has been closed (in https://github.com/JuliaOpt/Convex.jl/pull/322). This function no longer does anything and will be removed in a future Convex.jl release.", :clearmemory )
end

end
