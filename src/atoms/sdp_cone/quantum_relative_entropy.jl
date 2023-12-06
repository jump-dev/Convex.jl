#############################################################################
# quantum_relative_entropy returns tr(A*(log(A)-log(B))) where A and B
# are positive semidefinite matrices.  Note this function uses logarithm
# base e, not base 2, so return value is in units of nats, not bits.
#
# Quantum relative entropy is convex (jointly) in (A,B). This function
# implements the semidefinite programming approximation given in the reference
# below.  Parameters m and k control the accuracy of this approximation: m is
# the number of quadrature nodes to use and k the number of square-roots to
# take. See reference for more details.
#
# Implementation uses the expression
#    D(A||B) = e'*D_{op} (A \otimes I || I \otimes B) )*e
# where D_{op} is the operator relative entropy and e = vec(Matrix(I, n, n)).
#
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#
#REFERENCE
#   Ported from CVXQUAD which is based on the paper: "Lieb's concavity
#   theorem, matrix geometric means and semidefinite optimization" by Hamza
#   Fawzi and James Saunderson (arXiv:1512.03401)
#############################################################################

const MatrixOrConstant = Union{AbstractMatrix,Constant}

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
        v, U = eigen(Hermitian(B))
        if any(v .< -nullspace_tol)
            throw(DomainError(B, "B must be positive semidefinite"))
        end
        J = U'[v.>nullspace_tol, :]
        K = U'[v.<nullspace_tol, :]

        return new(children, (1, 1), m, k, B, J, K)
    end
end

head(io::IO, ::QuantumRelativeEntropy2) = print(io, "quantum_relative_entropy")

sign(atom::Union{QuantumRelativeEntropy1,QuantumRelativeEntropy2}) = Positive()

function monotonicity(atom::QuantumRelativeEntropy1)
    return (NoMonotonicity(), NoMonotonicity())
end
monotonicity(atom::QuantumRelativeEntropy2) = (NoMonotonicity(),)

function curvature(atom::Union{QuantumRelativeEntropy1,QuantumRelativeEntropy2})
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
    #println("quantum_relative_entropy general case")
    return QuantumRelativeEntropy1(A, B, m, k)
end

function quantum_relative_entropy(
    A::AbstractExpr,
    B::MatrixOrConstant,
    m::Integer = 3,
    k::Integer = 3,
    nullspace_tol::Real = 1e-6,
)
    #println("quantum_relative_entropy general case")
    return QuantumRelativeEntropy2(A, evaluate(B), m, k, nullspace_tol)
end

function quantum_relative_entropy(
    A::MatrixOrConstant,
    B::AbstractExpr,
    m::Integer = 3,
    k::Integer = 3,
)
    #println("quantum_relative_entropy constant A")
    A = evaluate(A)
    return -quantum_entropy(A, m, k) - trace_logm(B, A, m, k)
end

function quantum_relative_entropy(
    A::MatrixOrConstant,
    B::MatrixOrConstant,
    m::Integer = 0,
    k::Integer = 0,
    nullspace_tol::Real = 1e-6,
)
    #println("quantum_relative_entropy constant A and B")
    A = evaluate(A)
    B = evaluate(B)

    if size(A) != size(B)
        throw(DimensionMismatch("A and B must be the same size"))
    end
    if size(A) != (size(A)[1], size(A)[1])
        throw(DimensionMismatch("A and B must be square"))
    end
    if norm(A - A') > nullspace_tol
        throw(DomainError(A, "A must be Hermitian"))
    end
    if norm(B - B') > nullspace_tol
        throw(DomainError(B, "B must be Hermitian"))
    end

    # need to project down to support of A
    v, U = eigen(Hermitian(A))
    if any(v .< -nullspace_tol)
        throw(DomainError(A, "A must be positive semidefinite"))
    end
    if any(eigvals(Hermitian(B)) .< -nullspace_tol)
        throw(DomainError(B, "B must be positive semidefinite"))
    end
    J = U'[v.>nullspace_tol, :]
    Ap = Hermitian(J * A * J')
    Bp = Hermitian(J * B * J')

    if any(eigvals(Bp) .< nullspace_tol)
        return Inf
    end

    return real(tr(Ap * (log(Ap) - log(Bp))))
end

function _conic_form!(context::Context, atom::QuantumRelativeEntropy1)
    A = atom.children[1]
    B = atom.children[2]
    m = atom.m
    k = atom.k
    n = size(A)[1]
    eye = Matrix(1.0 * I, n, n)
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

function _conic_form!(context::Context, atom::QuantumRelativeEntropy2)
    A = atom.children[1]
    B = atom.B
    J = atom.J
    K = atom.K
    m = atom.m
    k = atom.k

    add_constraint!(context, A ⪰ 0)

    if length(K) > 0
        add_constraint!(context, K * A * K' == 0)
        Ap = J * A * J'
        Bp = Hermitian(J * B * J')
        τ = -quantum_entropy(Ap, m, k) - real(tr(Ap * log(Bp)))
    else
        τ = -quantum_entropy(A, m, k) - real(tr(A * log(B)))
    end

    return conic_form!(context, minimize(τ))
end
