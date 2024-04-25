# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
trace_logm(X, C) returns LinearAlgebra.tr(C*logm(X)) where X and C are a
positive definite matrices and C is constant.

trace_logm is concave in X.

This function implements the semidefinite programming approximation given
in the reference below.  Parameters m and k control the accuracy of this
approximation: m is the number of quadrature nodes to use and k the number
of square-roots to take. See reference for more details.

Implementation uses the expression
  LinearAlgebra.tr(C*logm(X)) = -LinearAlgebra.tr(C*D_{op}(I||X))
where D_{op} is the operator relative entropy:
  D_{op}(X||Y) = X^{1/2}*logm(X^{1/2} Y^{-1} X^{1/2})*X^{1/2}

All expressions and atoms are subtypes of AbstractExpr.
Please read expressions.jl first.

REFERENCE
  Ported from CVXQUAD which is based on the paper: "Lieb's concavity
  theorem, matrix geometric means and semidefinite optimization" by Hamza
  Fawzi and James Saunderson (arXiv:1512.03401)
"""
mutable struct TraceLogmAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    C::AbstractMatrix
    m::Integer
    k::Integer

    function TraceLogmAtom(
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
        if any(LinearAlgebra.eigvals(LinearAlgebra.Hermitian(C)) .< -1e-6)
            throw(DomainError(C, "C must be positive semidefinite"))
        end
        return new(children, (1, 1), C, m, k)
    end
end

head(io::IO, ::TraceLogmAtom) = print(io, "trace_logm")

Base.sign(::TraceLogmAtom) = NoSign()

monotonicity(::TraceLogmAtom) = (NoMonotonicity(),)

curvature(::TraceLogmAtom) = ConcaveVexity()

function evaluate(atom::TraceLogmAtom)
    return trace_logm(evaluate(atom.children[1]), atom.C)
end

function trace_logm(
    X::AbstractExpr,
    C::Union{AbstractMatrix,Constant},
    m::Integer = 3,
    k::Integer = 3,
)
    return TraceLogmAtom(X, evaluate(C), m, k)
end

function trace_logm(
    X::Union{AbstractMatrix,Constant},
    C::Union{AbstractMatrix,Constant},
    m::Integer = 3,
    k::Integer = 3,
)
    eye = Matrix(1.0 * LinearAlgebra.I, size(X))
    return -quantum_relative_entropy(C, X) + quantum_relative_entropy(C, eye)
end

function new_conic_form!(context::Context{T}, atom::TraceLogmAtom) where {T}
    X = atom.children[1]
    τ = if sign(X) == ComplexSign() || sign(constant(atom.C)) == ComplexSign()
        ComplexVariable(size(X))
    else
        Variable(size(X))
    end
    eye = Matrix(one(T) * LinearAlgebra.I, size(X))
    add_constraint!(
        context,
        τ in RelativeEntropyEpiCone(eye, X, atom.m, atom.k),
    )
    # It's already a real mathematically, but need to make it a real type.
    t = real(-LinearAlgebra.tr(atom.C * τ))
    return conic_form!(context, maximize(t))
end
