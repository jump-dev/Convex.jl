# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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

function new_conic_form!(
    context::Context{T},
    atom::QuantumRelativeEntropy1Atom,
) where {T}
    A, B = atom.children[1], atom.children[2]
    add_constraint!(context, A ⪰ 0)
    add_constraint!(context, B ⪰ 0)
    I = Matrix(one(T) * LinearAlgebra.I(size(A, 1)))
    τ, X, Y = Variable(), kron(A, I), kron(I, conj(B))
    set = RelativeEntropyEpiConeSquare(size(X, 1), atom.m, atom.k, vec(I))
    add_constraint!(context, Constraint((τ, X, Y), set))
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
