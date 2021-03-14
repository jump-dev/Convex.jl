#############################################################################
# sumlargesteigs.jl
# Handles top k eigenvalues of a symmetric or Hermitian matrix
# (and imposes the constraint that its argument be symmetric or Hermitian)
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
    return NoSign()
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

sumlargesteigs(x::AbstractExpr, k::Int) = k == 0 ? Constant(0) : SumLargestEigs(x, Constant(k))

# Create the equivalent conic problem:
#   minimize sk + Tr(Z)
#   subject to
#            Z ⪰ 0
#            A ⪰ 0
#            Z + sI ⪰ A
# See Ben-Tal and Nemirovski, "Lectures on Modern Convex Optimization"
# Example 18.c

function conic_form!(x::SumLargestEigs, unique_conic_forms)
    if !has_conic_form(unique_conic_forms, x)
        A = x.children[1]
        k = x.children[2]
        m, n = size(A)
        if sign(A) == ComplexSign()
            Z = ComplexVariable(n, n)
        else
            Z = Variable(n, n)
        end
        s = Variable()
        # The two inequality constraints have the side effect of constraining A to be symmetric,
        # since only symmetric matrices can be positive semidefinite.
        p = minimize(s*k + real(tr(Z)),
                     Z + s*Matrix(1.0I, n, n) - A ⪰ 0,
                     Z ⪰ 0)
        cache_conic_form!(unique_conic_forms, x, p)
    end
    return get_conic_form(unique_conic_forms, x)
end
