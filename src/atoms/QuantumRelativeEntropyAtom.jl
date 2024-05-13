# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    QuantumRelativeEntropy1Atom(
        A::AbstractExpr,
        B::AbstractExpr,
        m::Integer,
        k::Integer,
    )

quantum_relative_entropy returns LinearAlgebra.tr(A*(log(A)-log(B))) where A and
B are positive semidefinite matrices.

Note this function uses logarithm base e, not base 2, so return value is in
units of nats, not bits.

Quantum relative entropy is convex (jointly) in (A, B). This function implements
the semidefinite programming approximation given in the reference below.
Parameters m and k control the accuracy of this approximation: m is the number
of quadrature nodes to use and k the number of square-roots to take. See
reference for more details.

Implementation uses the expression

    D(A||B) = e'*D_{op} (A \\otimes I || I \\otimes B) )*e

where D_{op} is the operator relative entropy and e = vec(Matrix(I, n, n)).

## Reference

Ported from CVXQUAD which is based on the paper: "Lieb's concavity theorem,
matrix geometric means and semidefinite optimization" by Hamza Fawzi and James
Saunderson (arXiv:1512.03401)
"""
mutable struct QuantumRelativeEntropy1Atom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}
    m::Integer
    k::Integer

    function QuantumRelativeEntropy1Atom(
        A::AbstractExpr,
        B::AbstractExpr,
        m::Integer,
        k::Integer,
    )
        n = size(A, 1)
        if size(A) != size(B)
            throw(DimensionMismatch("A and B must be the same size"))
        elseif size(A) != (n, n)
            throw(DimensionMismatch("A and B must be square"))
        end
        return new((A, B), (1, 1), m, k)
    end
end

function head(io::IO, ::QuantumRelativeEntropy1Atom)
    return print(io, "quantum_relative_entropy")
end

Base.sign(::QuantumRelativeEntropy1Atom) = Positive()

function monotonicity(::QuantumRelativeEntropy1Atom)
    return (NoMonotonicity(), NoMonotonicity())
end

curvature(::QuantumRelativeEntropy1Atom) = ConvexVexity()

function evaluate(atom::QuantumRelativeEntropy1Atom)
    A = evaluate(atom.children[1])
    B = evaluate(atom.children[2])
    return quantum_relative_entropy(A, B)
end

function quantum_relative_entropy(
    A::AbstractExpr,
    B::AbstractExpr,
    m::Integer = 3,
    k::Integer = 3,
)
    return QuantumRelativeEntropy1Atom(A, B, m, k)
end

function new_conic_form!(
    context::Context{T},
    atom::QuantumRelativeEntropy1Atom,
) where {T}
    A, B = atom.children[1], atom.children[2]
    add_constraint!(context, A ⪰ 0)
    add_constraint!(context, B ⪰ 0)
    I = Matrix(one(T) * LinearAlgebra.I(size(A, 1)))
    m, k, e = atom.m, atom.k, vec(I)
    τ = Variable()
    add_constraint!(
        context,
        τ in RelativeEntropyEpiCone(kron(A, I), kron(I, conj(B)), m, k, e),
    )
    return conic_form!(context, τ)
end

mutable struct QuantumRelativeEntropy2Atom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    m::Integer
    k::Integer
    B::AbstractMatrix
    J::AbstractMatrix
    K::AbstractMatrix

    function QuantumRelativeEntropy2Atom(
        A::AbstractExpr,
        B::AbstractMatrix,
        m::Integer,
        k::Integer,
        nullspace_tol::Real,
    )
        n = size(A, 1)
        if size(A) != size(B)
            throw(DimensionMismatch("A and B must be the same size"))
        elseif size(A) != (n, n)
            throw(DimensionMismatch("A and B must be square"))
        elseif norm(B - B') > nullspace_tol
            throw(DomainError(B, "B must be Hermitian"))
        end
        # nullspace of A must contain nullspace of B
        v, U = LinearAlgebra.eigen(LinearAlgebra.Hermitian(B))
        if any(v .< -nullspace_tol)
            throw(DomainError(B, "B must be positive semidefinite"))
        end
        J = U'[v.>nullspace_tol, :]
        K = U'[v.<nullspace_tol, :]
        return new((A,), (1, 1), m, k, B, J, K)
    end
end

function head(io::IO, ::QuantumRelativeEntropy2Atom)
    return print(io, "quantum_relative_entropy")
end

Base.sign(::QuantumRelativeEntropy2Atom) = Positive()

monotonicity(::QuantumRelativeEntropy2Atom) = (NoMonotonicity(),)

curvature(::QuantumRelativeEntropy2Atom) = ConvexVexity()

function evaluate(atom::QuantumRelativeEntropy2Atom)
    A = evaluate(atom.children[1])
    return quantum_relative_entropy(A, atom.B)
end

function quantum_relative_entropy(
    A::AbstractExpr,
    B::Union{AbstractMatrix,Constant},
    m::Integer = 3,
    k::Integer = 3,
    nullspace_tol::Real = 1e-6,
)
    # This `evaluate` is safe since it is not a `fix!`ed variable
    # (it must be a constant or matrix)
    return QuantumRelativeEntropy2Atom(A, evaluate(B), m, k, nullspace_tol)
end

function quantum_relative_entropy(
    A::Union{AbstractMatrix,Constant},
    B::AbstractExpr,
    m::Integer = 3,
    k::Integer = 3,
)
    # This `evaluate` is safe since it is not a `fix!`ed variable
    # (it must be a constant or matrix)
    A = evaluate(A)
    return -quantum_entropy(A, m, k) - trace_logm(B, A, m, k)
end

function quantum_relative_entropy(
    A::Union{AbstractMatrix,Constant},
    B::Union{AbstractMatrix,Constant},
    m::Integer = 0,
    k::Integer = 0,
    nullspace_tol::Real = 1e-6,
)
    # These `evaluate`s are safe since it is not a `fix!`ed variable
    # (it must be a constant or matrix)
    A, B = evaluate(A), evaluate(B)
    if size(A) != size(B)
        throw(DimensionMismatch("A and B must be the same size"))
    elseif size(A) != (size(A)[1], size(A)[1])
        throw(DimensionMismatch("A and B must be square"))
    elseif norm(A - A') > nullspace_tol
        throw(DomainError(A, "A must be Hermitian"))
    elseif norm(B - B') > nullspace_tol
        throw(DomainError(B, "B must be Hermitian"))
    end
    v, U = LinearAlgebra.eigen(LinearAlgebra.Hermitian(A))
    if any(v .< -nullspace_tol)
        throw(DomainError(A, "A must be positive semidefinite"))
    end
    if any(LinearAlgebra.eigvals(LinearAlgebra.Hermitian(B)) .< -nullspace_tol)
        throw(DomainError(B, "B must be positive semidefinite"))
    end
    J = U'[v.>nullspace_tol, :]
    Ap = LinearAlgebra.Hermitian(J * A * J')
    Bp = LinearAlgebra.Hermitian(J * B * J')
    if any(LinearAlgebra.eigvals(Bp) .< nullspace_tol)
        return Inf
    end
    return real(LinearAlgebra.tr(Ap * (log(Ap) - log(Bp))))
end

function new_conic_form!(context::Context, atom::QuantumRelativeEntropy2Atom)
    A = only(atom.children)
    B, J, K, m, k = atom.B, atom.J, atom.K, atom.m, atom.k
    add_constraint!(context, A ⪰ 0)
    τ = if length(K) > 0
        add_constraint!(context, K * A * K' == 0)
        Ap = J * A * J'
        Bp = LinearAlgebra.Hermitian(J * B * J')
        -quantum_entropy(Ap, m, k) - real(LinearAlgebra.tr(Ap * log(Bp)))
    else
        -quantum_entropy(A, m, k) - real(LinearAlgebra.tr(A * log(B)))
    end
    return conic_form!(context, τ)
end
