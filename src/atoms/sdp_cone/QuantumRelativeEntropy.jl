"""
quantum_relative_entropy returns LinearAlgebra.tr(A*(log(A)-log(B))) where A and B
are positive semidefinite matrices.  Note this function uses logarithm
base e, not base 2, so return value is in units of nats, not bits.

Quantum relative entropy is convex (jointly) in (A,B). This function
implements the semidefinite programming approximation given in the reference
below.  Parameters m and k control the accuracy of this approximation: m is
the number of quadrature nodes to use and k the number of square-roots to
take. See reference for more details.

Implementation uses the expression
   D(A||B) = e'*D_{op} (A \\otimes I || I \\otimes B) )*e
where D_{op} is the operator relative entropy and e = vec(Matrix(I, n, n)).

All expressions and atoms are subtypes of AbstractExpr.
Please read expressions.jl first.

REFERENCE
  Ported from CVXQUAD which is based on the paper: "Lieb's concavity
  theorem, matrix geometric means and semidefinite optimization" by Hamza
  Fawzi and James Saunderson (arXiv:1512.03401)
"""
mutable struct QuantumRelativeEntropy1 <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}
    m::Integer
    k::Integer

    function QuantumRelativeEntropy1(
        A::AbstractExpr,
        B::AbstractExpr,
        m::Integer,
        k::Integer,
    )
        children = (A, B)
        if size(A) != size(B)
            throw(DimensionMismatch("A and B must be the same size"))
        end
        n = size(A)[1]
        if size(A) != (n, n)
            throw(DimensionMismatch("A and B must be square"))
        end
        return new(children, (1, 1), m, k)
    end
end

head(io::IO, ::QuantumRelativeEntropy1) = print(io, "quantum_relative_entropy")

mutable struct QuantumRelativeEntropy2 <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    m::Integer
    k::Integer
    B::AbstractMatrix
    J::AbstractMatrix
    K::AbstractMatrix

    function QuantumRelativeEntropy2(
        A::AbstractExpr,
        B::AbstractMatrix,
        m::Integer,
        k::Integer,
        nullspace_tol::Real,
    )
        children = (A,)

        if size(A) != size(B)
            throw(DimensionMismatch("A and B must be the same size"))
        end
        n = size(A)[1]
        if size(A) != (n, n)
            throw(DimensionMismatch("A and B must be square"))
        end
        if norm(B - B') > nullspace_tol
            throw(DomainError(B, "B must be Hermitian"))
        end
        # nullspace of A must contain nullspace of B
        v, U = LinearAlgebra.eigen(LinearAlgebra.Hermitian(B))
        if any(v .< -nullspace_tol)
            throw(DomainError(B, "B must be positive semidefinite"))
        end
        J = U'[v.>nullspace_tol, :]
        K = U'[v.<nullspace_tol, :]
        return new(children, (1, 1), m, k, B, J, K)
    end
end

head(io::IO, ::QuantumRelativeEntropy2) = print(io, "quantum_relative_entropy")

function Base.sign(::Union{QuantumRelativeEntropy1,QuantumRelativeEntropy2})
    return Positive()
end

function monotonicity(::QuantumRelativeEntropy1)
    return (NoMonotonicity(), NoMonotonicity())
end

monotonicity(::QuantumRelativeEntropy2) = (NoMonotonicity(),)

function curvature(::Union{QuantumRelativeEntropy1,QuantumRelativeEntropy2})
    return ConvexVexity()
end

function evaluate(atom::QuantumRelativeEntropy1)
    A = evaluate(atom.children[1])
    B = evaluate(atom.children[2])
    return quantum_relative_entropy(A, B)
end

function evaluate(atom::QuantumRelativeEntropy2)
    A = evaluate(atom.children[1])
    return quantum_relative_entropy(A, atom.B)
end

function quantum_relative_entropy(
    A::AbstractExpr,
    B::AbstractExpr,
    m::Integer = 3,
    k::Integer = 3,
)
    return QuantumRelativeEntropy1(A, B, m, k)
end

function quantum_relative_entropy(
    A::AbstractExpr,
    B::Union{AbstractMatrix,Constant},
    m::Integer = 3,
    k::Integer = 3,
    nullspace_tol::Real = 1e-6,
)
    return QuantumRelativeEntropy2(A, evaluate(B), m, k, nullspace_tol)
end

function quantum_relative_entropy(
    A::Union{AbstractMatrix,Constant},
    B::AbstractExpr,
    m::Integer = 3,
    k::Integer = 3,
)
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
    A = evaluate(A)
    B = evaluate(B)
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

function new_conic_form!(context::Context, atom::QuantumRelativeEntropy1)
    A = atom.children[1]
    B = atom.children[2]
    m = atom.m
    k = atom.k
    n = size(A)[1]
    eye = Matrix(1.0 * LinearAlgebra.I, n, n)
    e = vec(eye)
    add_constraint!(context, A ⪰ 0)
    add_constraint!(context, B ⪰ 0)
    τ = Variable()
    add_constraint!(
        context,
        τ in RelativeEntropyEpiCone(kron(A, eye), kron(eye, conj(B)), m, k, e),
    )
    return conic_form!(context, minimize(τ))
end

function new_conic_form!(context::Context, atom::QuantumRelativeEntropy2)
    A = atom.children[1]
    B = atom.B
    J = atom.J
    K = atom.K
    m = atom.m
    k = atom.k
    add_constraint!(context, A ⪰ 0)
    τ = if length(K) > 0
        add_constraint!(context, K * A * K' == 0)
        Ap = J * A * J'
        Bp = LinearAlgebra.Hermitian(J * B * J')
        -quantum_entropy(Ap, m, k) - real(LinearAlgebra.tr(Ap * log(Bp)))
    else
        -quantum_entropy(A, m, k) - real(LinearAlgebra.tr(A * log(B)))
    end
    return conic_form!(context, minimize(τ))
end
