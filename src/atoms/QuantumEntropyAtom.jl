# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
quantum_entropy returns -LinearAlgebra.tr(X*log(X)) where X is a positive semidefinite.
Note this function uses logarithm base e, not base 2, so return value is in
units of nats, not bits.

Quantum entropy is concave. This function implements the semidefinite
programming approximation given in the reference below.  Parameters m and k
control the accuracy of this approximation: m is the number of quadrature
nodes to use and k the number of square-roots to take. See reference for
more details.

Implementation uses the expression
  H(X) = -trace( D_{op}(X||I) )
where D_{op} is the operator relative entropy:
  D_{op}(X||Y) = X^{1/2}*logm(X^{1/2} Y^{-1} X^{1/2})*X^{1/2}

All expressions and atoms are subtypes of AbstractExpr.
Please read expressions.jl first.

REFERENCE
  Ported from CVXQUAD which is based on the paper: "Lieb's concavity
  theorem, matrix geometric means and semidefinite optimization" by Hamza
  Fawzi and James Saunderson (arXiv:1512.03401)
"""
mutable struct QuantumEntropyAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    m::Integer
    k::Integer

    function QuantumEntropyAtom(X::AbstractExpr, m::Integer, k::Integer)
        children = (X,)
        n = size(X)[1]
        if size(X) != (n, n)
            throw(DimensionMismatch("X must be square"))
        end
        return new(children, (1, 1), m, k)
    end
end

head(io::IO, ::QuantumEntropyAtom) = print(io, "quantum_entropy")

Base.sign(::QuantumEntropyAtom) = Positive()

monotonicity(::QuantumEntropyAtom) = (NoMonotonicity(),)

curvature(::QuantumEntropyAtom) = ConcaveVexity()

function evaluate(atom::QuantumEntropyAtom)
    return quantum_entropy(evaluate(atom.children[1]))
end

function quantum_entropy(X::AbstractExpr, m::Integer = 3, k::Integer = 3)
    return QuantumEntropyAtom(X, m, k)
end

function quantum_entropy(
    X::Union{AbstractMatrix,Constant},
    m::Integer = 0,
    k::Integer = 0,
)
    return -quantum_relative_entropy(X, Matrix(1.0 * LinearAlgebra.I, size(X)))
end

function new_conic_form!(context::Context, atom::QuantumEntropyAtom)
    X = atom.children[1]
    n = size(X, 1)
    eye = Matrix(1.0 * LinearAlgebra.I, n, n)
    add_constraint!(context, X ⪰ 0)
    τ = if sign(X) == ComplexSign()
        ComplexVariable(n, n)
    else
        Variable(n, n)
    end
    add_constraint!(
        context,
        τ in RelativeEntropyEpiCone(X, eye, atom.m, atom.k),
    )
    # It's already a real mathematically, but need to make it a real type.
    τ = real(-LinearAlgebra.tr(τ))
    return conic_form!(context, minimize(τ))
end
