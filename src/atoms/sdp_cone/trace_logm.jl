#############################################################################
# trace_logm(X, C) returns tr(C*logm(X)) where X and C are a positive definite
# matrices and C is constant.
#
# trace_logm is concave in X.
#
# This function implements the semidefinite programming approximation given
# in the reference below.  Parameters m and k control the accuracy of this
# approximation: m is the number of quadrature nodes to use and k the number
# of square-roots to take. See reference for more details.
#
# Implementation uses the expression
#   tr(C*logm(X)) = -tr(C*D_{op}(I||X))
# where D_{op} is the operator relative entropy:
#   D_{op}(X||Y) = X^{1/2}*logm(X^{1/2} Y^{-1} X^{1/2})*X^{1/2}
#
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#
#REFERENCE
#   Ported from CVXQUAD which is based on the paper: "Lieb's concavity
#   theorem, matrix geometric means and semidefinite optimization" by Hamza
#   Fawzi and James Saunderson (arXiv:1512.03401)
#############################################################################

struct TraceLogm <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    C::AbstractMatrix
    m::Integer
    k::Integer

    function TraceLogm(
        X::AbstractExpr,
        C::AbstractMatrix,
        m::Integer,
        k::Integer,
    )
        children = (X,)
        if size(X) != size(C)
            throw(DimensionMismatch("X and C must be the same size"))
        end
        n = size(X)[1]
        if size(X) != (n, n)
            throw(DimensionMismatch("X and C must be square"))
        end
        if norm(C - C') > 1e-6
            throw(DomainError(C, "C must be Hermitian"))
        end
        if any(eigvals(Hermitian(C)) .< -1e-6)
            throw(DomainError(C, "C must be positive semidefinite"))
        end
        return new(:trace_logm, hash(children), children, (1, 1), C, m, k)
    end
end

function sign(atom::TraceLogm)
    return NoSign()
end

function monotonicity(atom::TraceLogm)
    return (NoMonotonicity(),)
end

function curvature(atom::TraceLogm)
    return ConcaveVexity()
end

function evaluate(atom::TraceLogm)
    X = evaluate(atom.children[1])
    return trace_logm(X, atom.C)
end

const MatrixOrConstant = Union{AbstractMatrix,Constant}

function trace_logm(
    X::AbstractExpr,
    C::MatrixOrConstant,
    m::Integer = 3,
    k::Integer = 3,
)
    #println("trace_logm general case")
    return TraceLogm(X, evaluate(C), m, k)
end

function trace_logm(
    X::MatrixOrConstant,
    C::MatrixOrConstant,
    m::Integer = 3,
    k::Integer = 3,
)
    #println("trace_logm constant X")
    eye = Matrix(1.0 * I, size(X))
    return -quantum_relative_entropy(C, X) + quantum_relative_entropy(C, eye)
end

function conic_form!(context::Context, atom::TraceLogm)
    X = atom.children[1]
    C = atom.C
    m = atom.m
    k = atom.k
    eye = Matrix(1.0 * I, size(X))

    is_complex = sign(X) == ComplexSign() || sign(constant(C)) == ComplexSign()
    if is_complex
        τ = ComplexVariable(size(X))
    else
        τ = Variable(size(X))
    end
    add_constraints_to_context(
        τ in RelativeEntropyEpiCone(eye, X, m, k),
        context,
    )

    # It's already a real mathematically, but need to make it a real type.
    t = real(-tr(C * τ))
    return conic_form!(context, maximize(t))
end
