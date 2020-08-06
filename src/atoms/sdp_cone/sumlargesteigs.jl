#############################################################################
# sumlargesteigs.jl
# Handles top k eigenvalues of a symmetric positive definite matrix
# (and imposes the constraint that its argument be PSD)
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import LinearAlgebra.eigvals

### sumlargesteigs

struct SumLargestEigs <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}

    function SumLargestEigs(x::AbstractExpr, k::AbstractExpr)
        children = (x, k)
        m, n = size(x)
        if m == n
            return new(:sumlargesteigs, hash(children), children, (1,1))
        else
            error("sumlargesteigs can only be applied to a square matrix.")
        end
    end
end

function sign(x::SumLargestEigs)
    return Positive()
end

function monotonicity(x::SumLargestEigs)
    return (Nondecreasing(), NoMonotonicity())
end

function curvature(x::SumLargestEigs)
    return ConvexVexity()
end

function evaluate(x::SumLargestEigs)
    eigvals(evaluate(x.children[1]))[end-x.children[2]:end]
end

sumlargesteigs(x::AbstractExpr, k::Int) = SumLargestEigs(x, Constant(k))

# Create the equivalent conic problem:
#   minimize sk + Tr(Z)
#   subject to
#            Z ⪰ 0
#            A ⪰ 0
#            Z + sI ⪰ A
# See Ben-Tal and Nemirovski, "Lectures on Modern Convex Optimization"
# Example 18.c
function template(x::SumLargestEigs, context::Context{T}) where T
    A = x.children[1]
    k = x.children[2]
    m, n = size(A)
    Z = Variable(n, n)
    s = Variable()
    p = minimize(s*k + tr(Z),
                    Z + s*Matrix(1.0I, n, n) - A ⪰ 0,
                    A ⪰ 0, Z ⪰ 0)
    return template(p, context)
end
