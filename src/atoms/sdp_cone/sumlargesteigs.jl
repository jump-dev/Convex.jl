#############################################################################
# sumlargesteigs.jl
# Handles top k eigenvalues of a symmetric or Hermitian matrix
# (and imposes the constraint that its argument be symmetric or Hermitian)
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import LinearAlgebra.eigvals

### sumlargesteigs

mutable struct SumLargestEigs <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function SumLargestEigs(x::AbstractExpr, k::AbstractExpr)
        children = (x, k)
        m, n = size(x)
        if m == n
            return new(children, (1, 1))
        else
            error("sumlargesteigs can only be applied to a square matrix.")
        end
    end
end

head(io::IO, ::SumLargestEigs) = print(io, "sumlargesteigs")

function sign(x::SumLargestEigs)
    return NoSign()
end

function monotonicity(x::SumLargestEigs)
    return (Nondecreasing(), NoMonotonicity())
end

function curvature(x::SumLargestEigs)
    return ConvexVexity()
end

function evaluate(x::SumLargestEigs)
    return eigvals(evaluate(x.children[1]))[end-x.children[2]:end]
end

function sumlargesteigs(x::AbstractExpr, k::Int)
    return k == 0 ? Constant(0) : SumLargestEigs(x, Constant(k))
end

# Create the equivalent conic problem:
#   minimize sk + Tr(Z)
#   subject to
#            Z ⪰ 0
#            A ⪰ 0
#            Z + sI ⪰ A
# See Ben-Tal and Nemirovski, "Lectures on Modern Convex Optimization"
# Example 18.c
function _conic_form!(context::Context{T}, x::SumLargestEigs) where {T}
    X = x.children[1]
    k = x.children[2]
    m, n = size(X)
    if iscomplex(sign(X))
        Z = ComplexVariable(n, n)
    else
        Z = Variable(n, n)
    end
    s = Variable()
    # Note: we know the trace is real, since Z is PSD, but we need to tell Convex.jl that.
    p = minimize(
        s * k + real(tr(Z)),
        [Z - X + s * Matrix(1.0I, n, n) ⪰ 0, Z ⪰ 0, X == X'],
    )
    return conic_form!(context, p)
end
