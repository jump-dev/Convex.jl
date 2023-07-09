#############################################################################
# quantum_entropy returns -tr(X*log(X)) where X is a positive semidefinite.
# Note this function uses logarithm base e, not base 2, so return value is in
# units of nats, not bits.
#
# Quantum entropy is concave. This function implements the semidefinite
# programming approximation given in the reference below.  Parameters m and k
# control the accuracy of this approximation: m is the number of quadrature
# nodes to use and k the number of square-roots to take. See reference for
# more details.
#
# Implementation uses the expression
#   H(X) = -trace( D_{op}(X||I) )
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

struct QuantumEntropy <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    m::Integer
    k::Integer

    function QuantumEntropy(X::AbstractExpr, m::Integer, k::Integer)
        children = (X,)
        n = size(X)[1]
        if size(X) != (n, n)
            throw(DimensionMismatch("X must be square"))
        end
        return new(:quantum_entropy, hash(children), children, (1, 1), m, k)
    end
end

function sign(atom::QuantumEntropy)
    return Positive()
end

function monotonicity(atom::QuantumEntropy)
    return (NoMonotonicity(),)
end

function curvature(atom::QuantumEntropy)
    return ConcaveVexity()
end

function evaluate(atom::QuantumEntropy)
    X = evaluate(atom.children[1])
    return quantum_entropy(X)
end

const MatrixOrConstant = Union{AbstractMatrix,Constant}

function quantum_entropy(X::AbstractExpr, m::Integer = 3, k::Integer = 3)
    #println("quantum_entropy general case")
    return QuantumEntropy(X, m, k)
end

function quantum_entropy(X::MatrixOrConstant, m::Integer = 0, k::Integer = 0)
    #println("quantum_entropy constant X")
    return -quantum_relative_entropy(X, Matrix(1.0 * I, size(X)))
end

function conic_form!(context::Context, atom::QuantumEntropy)
    X = atom.children[1]
    m = atom.m
    k = atom.k
    n = size(X)[1]
    eye = Matrix(1.0 * I, n, n)

    add_constraints_to_context(X ⪰ 0, context)

    is_complex = sign(X) == ComplexSign()
    if is_complex
        τ = ComplexVariable(n, n)
    else
        τ = Variable(n, n)
    end
    add_constraints_to_context(
        τ in RelativeEntropyEpiCone(X, eye, m, k),
        context,
    )

    # It's already a real mathematically, but need to make it a real type.
    τ = real(-tr(τ))
    return conic_form!(context, minimize(τ))
end
