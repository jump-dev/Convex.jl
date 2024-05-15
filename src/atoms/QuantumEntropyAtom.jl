# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    QuantumEntropyAtom(X::AbstractExpr, m::Integer, k::Integer)

quantum_entropy returns -LinearAlgebra.tr(X*log(X)) where X is a positive
semidefinite.

Note this function uses logarithm base e, not base 2, so return value is in
units of nats, not bits.

Quantum entropy is concave. This function implements the semidefinite
programming approximation given in the reference below.  Parameters m and k
control the accuracy of this approximation: m is the number of quadrature nodes
to use and k the number of square-roots to take. See reference for more details.

Implementation uses the expression

    H(X) = -trace( D_{op}(X||I) )

where D_{op} is the operator relative entropy:

    D_{op}(X||Y) = X^{1/2}*logm(X^{1/2} Y^{-1} X^{1/2})*X^{1/2}

## Reference

Ported from CVXQUAD which is based on the paper: "Lieb's concavity theorem,
matrix geometric means and semidefinite optimization" by Hamza Fawzi and James
Saunderson (arXiv:1512.03401)
"""
mutable struct QuantumEntropyAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    m::Integer
    k::Integer

    function QuantumEntropyAtom(X::AbstractExpr, m::Integer, k::Integer)
        n = size(X, 1)
        if size(X) != (n, n)
            throw(DimensionMismatch("X must be square"))
        end
        return new((X,), (1, 1), m, k)
    end
end

head(io::IO, ::QuantumEntropyAtom) = print(io, "quantum_entropy")

Base.sign(::QuantumEntropyAtom) = Positive()

monotonicity(::QuantumEntropyAtom) = (NoMonotonicity(),)

curvature(::QuantumEntropyAtom) = ConcaveVexity()

evaluate(atom::QuantumEntropyAtom) = quantum_entropy(evaluate(atom.children[1]))

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

function new_conic_form!(
    context::Context{T},
    atom::QuantumEntropyAtom,
) where {T}
    X = only(atom.children)
    add_constraint!(context, X ⪰ 0)
    τ = if sign(X) == ComplexSign()
        ComplexVariable(size(X))
    else
        Variable(size(X))
    end
    I = Matrix(one(T) * LinearAlgebra.I(size(X, 1)))
    set = RelativeEntropyEpiConeSquare(size(X, 1), atom.m, atom.k)
    add_constraint!(context, GenericConstraint((τ, X, I), set))
    # It's already a real mathematically, but need to make it a real type.
    return conic_form!(context, real(-LinearAlgebra.tr(τ)))
end
