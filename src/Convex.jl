__precompile__()

module Convex
using OrderedCollections: OrderedDict
using LinearAlgebra
using SparseArrays
using AbstractTrees: AbstractTrees, children

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

# Functions
export conv, dotsort, entropy, exp, geomean, hinge_loss, huber, inner_product, invpos
export log_perspective, logisticloss, logsumexp, matrixfrac, neg, norm2, norm_1, norm_inf, nuclearnorm
export partialtrace, partialtranspose, pos, qol_elementwise, quadform, quadoverlin, rationalnorm
export relative_entropy, scaledgeomean, sigmamax, square, sumlargest, sumlargesteigs, sumsmallest, sumsquares

# rexports from LinearAlgebra
export diag, diagm, Diagonal, dot, eigmax, eigmin, kron, logdet, norm, tr

# Constraints
export Constraint
export isposdef, ⪰, ⪯ # PSD constraints
export socp
export Constraint # useful for making abstractly-typed vectors via `Constraint[]`

# Variables
export constant, ComplexVariable, HermitianSemidefinite, Semidefinite, Variable
export curvature, evaluate, fix!, free!, monotonicity, sign, vexity
export BinVar, IntVar, ContVar, vartype, vartype!
export constraints, add_constraint!, set_value!, evaluate

# Signs
export Positive, Negative, ComplexSign, NoSign

# Problems
export add_constraints!, maximize, minimize, Problem, satisfy, solve!, solve2!


# Module level globals

"""
    DCP_WARNINGS

Controls whether or not warnings are emitted for when an expression fails to be
of disciplined convex form. To turn warnings off, run

    Convex.DCP_WARNINGS[] = false
"""
const DCP_WARNINGS = Ref(true)

"""
    MAXDEPTH

Controls depth of tree printing globally for Convex.jl; defaults to 3. Set via

    Convex.MAXDEPTH[] = 5
"""
const MAXDEPTH = Ref(3)

"""
    MAXWIDTH

Controls width of tree printing globally for Convex.jl; defaults to 15. Set via

    Convex.MAXWIDTH[] = 15
"""
const MAXWIDTH= Ref(15)

"""
    MAXDIGITS

When priting IDs of variables, only show the initial and final digits
if the full ID has more than double the number of digits specified
here.  So, with the default setting MAXDIGITS=3, any ID longer than 7
digits would be shortened; for example, ID `14656210999710729289`
would be printed as `146…289`.

This setting controls tree printing globally for Convex.jl; defaults to 3.

Set via:

    Convex.MAXDIGITS[] = 3
"""
const MAXDIGITS= Ref(3)


# where do these go?
# used so far only in `Constant`
vectorize(v::AbstractVector) = v
vectorize(v::Number) = [v]
vectorize(v::AbstractMatrix) = vec(v)

# where should these go?
function vec_triu(M)
    L = LinearIndices(size(M))
    n, m = size(M)
    inds = [ L[i,j] for i = 1:n for j = i:m ]
    return M[inds]
end

function vec_tril(M)
    L = LinearIndices(size(M))
    n, m = size(M)
    inds = [ L[i,j]  for i = 1:n for j = 1:i ]
    return M[inds]
end


include("Context.jl")
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
include("constraints/soc_constraints.jl")
include("constraints/exp_constraints.jl")
include("constraints/sdp_constraints.jl")
include("problems.jl")
include("solution.jl")
include("VectorAffineFunctionAsMatrix.jl")
include("complex.jl")
include("solve2!.jl")

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
include("atoms/sdp_cone/eig_min_max.jl")
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


# Deprecated workaround for memory leak (https://github.com/jump-dev/Convex.jl/issues/83)
function clearmemory()
    Base.depwarn("Convex.clearmemory() is deprecated, as the memory leak it works around has been closed (in https://github.com/jump-dev/Convex.jl/pull/322). This function no longer does anything and will be removed in a future Convex.jl release.", :clearmemory )
end

@deprecate lambdamin eigmin
@deprecate lambdamax eigmax

end
