__precompile__()

module Convex
using OrderedCollections: OrderedDict
using LinearAlgebra
using SparseArrays
using LDLFactorizations
using AbstractTrees: AbstractTrees, children

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

# Functions
export conv,
    dotsort, entropy, exp, geomean, hinge_loss, huber, inner_product, invpos
export log_perspective,
    logisticloss,
    logsumexp,
    matrixfrac,
    neg,
    norm2,
    norm_1,
    norm_inf,
    nuclearnorm
export partialtrace,
    partialtranspose, pos, qol_elementwise, quadform, quadoverlin, rationalnorm
export relative_entropy,
    sigmamax, square, sumlargest, sumlargesteigs, sumsmallest, sumsquares
export GeomMeanHypoCone,
    GeomMeanEpiCone,
    RelativeEntropyEpiCone,
    quantum_relative_entropy,
    quantum_entropy
export trace_logm, trace_mpower, lieb_ando

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
export add_constraints!, maximize, minimize, Problem, satisfy, solve!

# Module level globals

"""
    emit_dcp_warnings

Controls whether or not warnings are emitted for when an expression fails to be
of disciplined convex form. To turn warnings off, override the method via

    Convex.emit_dcp_warnings() = false

This will cause Julia's method invalidation to recompile any functions emitting
DCP warnings and remove them. This should be run from top-level (not within a function).
"""
emit_dcp_warnings() = true

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
const MAXWIDTH = Ref(15)

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
const MAXDIGITS = Ref(3)

# where do these go?
# used so far only in `Constant`
vectorize(v::AbstractVector) = v
vectorize(v::Number) = [v]
vectorize(v::AbstractMatrix) = vec(v)

# where should these go?
function vec_triu(M)
    L = LinearIndices(size(M))
    n, m = size(M)
    inds = [L[i, j] for i in 1:n for j in i:m]
    return M[inds]
end

function vec_tril(M)
    L = LinearIndices(size(M))
    n, m = size(M)
    inds = [L[i, j] for i in 1:n for j in 1:i]
    return M[inds]
end

include("Context.jl")
### modeling framework
include("dcp.jl")
include("expressions.jl")
include("variable.jl")
include("variable_template.jl")
include("constant.jl")
include("constraints/constraints.jl")
include("constraints/soc_constraints.jl")
include("constraints/exp_constraints.jl")
include("constraints/sdp_constraints.jl")
include("problems.jl")
include("VAFTape.jl")
include("SparseVAFTape.jl")
include("VectorAffineFunctionAsMatrix.jl")
include("complex.jl")
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
include("atoms/second_order_cone/rationalnorm.jl")
include("atoms/second_order_cone/huber.jl")

### SDP atoms
include("atoms/sdp_cone/nuclearnorm.jl")
include("atoms/sdp_cone/operatornorm.jl")
include("atoms/sdp_cone/eig_min_max.jl")
include("atoms/sdp_cone/matrixfrac.jl")
include("atoms/sdp_cone/sumlargesteigs.jl")
include("atoms/sdp_cone/geom_mean_hypocone.jl")
include("atoms/sdp_cone/geom_mean_epicone.jl")
include("atoms/sdp_cone/relative_entropy_epicone.jl")
include("atoms/sdp_cone/quantum_relative_entropy.jl")
include("atoms/sdp_cone/quantum_entropy.jl")
include("atoms/sdp_cone/trace_logm.jl")
include("atoms/sdp_cone/trace_mpower.jl")
include("atoms/sdp_cone/lieb_ando.jl")

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

end
