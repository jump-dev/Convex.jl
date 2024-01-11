module Convex

import AbstractTrees
import LDLFactorizations
import LinearAlgebra
import MathOptInterface as MOI
import OrderedCollections
import SparseArrays

export conv,
    dotsort,
    entropy,
    exp,
    geomean,
    hinge_loss,
    huber,
    inner_product,
    invpos,
    log_perspective,
    logisticloss,
    logsumexp,
    matrixfrac,
    neg,
    norm_1,
    norm_inf,
    nuclearnorm,
    partialtrace,
    partialtranspose,
    pos,
    qol_elementwise,
    quadform,
    quadoverlin,
    rationalnorm,
    relative_entropy,
    sigmamax,
    square,
    sumlargest,
    sumlargesteigs,
    sumsmallest,
    sumsquares,
    GeomMeanHypoCone,
    GeomMeanEpiCone,
    RelativeEntropyEpiCone,
    quantum_relative_entropy,
    quantum_entropy,
    trace_logm,
    trace_mpower,
    lieb_ando,
    # Constraints
    Constraint,
    isposdef,
    ⪰,
    ⪯,
    socp,
    Constraint, # useful for making abstractly-typed vectors via `Constraint[]`
    # Variables
    constant,
    ComplexVariable,
    HermitianSemidefinite,
    Semidefinite,
    Variable,
    curvature,
    evaluate,
    fix!,
    free!,
    monotonicity,
    vexity,
    problem_vexity,
    BinVar,
    IntVar,
    ContVar,
    vartype,
    vartype!,
    get_constraints,
    add_constraint!,
    set_value!,
    evaluate,
    # Signs
    Positive,
    Negative,
    ComplexSign,
    NoSign,
    # Problems
    add_constraints!,
    maximize,
    minimize,
    Problem,
    satisfy,
    solve!,
    write_to_file,
    DCPViolationError

# Imports and exports as needed to maintain backwards compatibility
for k in (
    :diag,
    :diagm,
    :Diagonal,
    :dot,
    :eigmax,
    :eigmin,
    :logdet,
    :norm,
    :norm2,
    :tr,
)
    @eval begin
        using LinearAlgebra: $k
        export $k
    end
end
using AbstractTrees: children
using Base: sign

# Module level globals

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

# using SuiteSparseGraphBLAS
#

# vec(x) = Base.vec(x)
# function vec(x::GBMatrix)
#     # Hacks to try to get `vec` to work
#     x = reshape(x, length(x), 1)
#     return x[:, 1]
# end

# blockdiag(xs...) = SparseArrays.blockdiag(xs...)::SPARSE_MATRIX

# function blockdiag(xs::GBMatrix{T,T}...) where {T}
#     N = length(xs)
#     entries = Matrix{GBMatrix{T,T}}(undef, N, N)
#     heights = size.(xs, 1)
#     for (i, x) in enumerate(xs)
#         entries[i, i] = x
#         m = size(x, 2)
#         for j in 1:(i-1)
#             entries[j, i] = GBMatrix{T,T}(heights[j], m)
#         end
#         for j in (i+1):lastindex(entries, 1)
#             entries[j, i] = GBMatrix{T,T}(heights[j], m)
#         end
#     end
#     return cat(entries)
# end
#
# const SPARSE_VECTOR{T} = GBVector{T,T}
# const SPARSE_MATRIX{T} = GBMatrix{T,T}
# spzeros(T, d) = GBVector{T,T}(d)
# spzeros(T, n, m) = GBMatrix{T,T}(n, m)
# spidentity(T, d) = GBMatrix{T,T}(LinearAlgebra.Diagonal(ones(T, d)))
# create_sparse(T, args...) = GBMatrix{T,T}(args...)

const SPARSE_VECTOR{T} = Vector{T}
const SPARSE_MATRIX{T} = SparseArrays.SparseMatrixCSC{T,Int}
spzeros(T, d) = zeros(T, d)
spzeros(T, n, m) = SparseArrays.spzeros(T, n, m)
spidentity(T, d) = SparseArrays.sparse(one(T) * LinearAlgebra.I, d, d)
function create_sparse(::Type{T}, args...) where {T}
    local result::SPARSE_MATRIX{T}
    result = SparseArrays.sparse(args...)
    return result
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
include("SparseTape.jl")
include("VectorAffineFunctionAsMatrix.jl")
include("ComplexTape.jl")
include("operate.jl")
include("complex_operate.jl")
include("real_operate.jl")
include("solution.jl")
include("MOI_wrapper.jl")

### atoms
for (root, _, files) in walkdir(joinpath(@__DIR__, "atoms"))
    for file in files
        include(joinpath(root, file))
    end
end

### utilities
include("utilities/tree_print.jl")
include("utilities/tree_interface.jl")
include("utilities/show.jl")
include("utilities/iteration.jl")
include("problem_depot/problem_depot.jl")

end
