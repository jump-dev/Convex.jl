# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
Returns LinearAlgebra.tr(K' * A^{1-t} * K * B^t) where A and B are positive semidefinite
matrices and K is an arbitrary matrix (possibly rectangular).

Disciplined convex programming information:
   lieb_ando(A,B,K,t) is concave in (A,B) for t in [0,1], and convex
   in (A,B) for t in [-1,0] or [1,2]. K is a fixed matrix.

Seems numerically unstable when t is on the endpoints of these ranges.

All expressions and atoms are subtypes of AbstractExpr.
Please read expressions.jl first.

REFERENCE
  Ported from CVXQUAD which is based on the paper: "Lieb's concavity
  theorem, matrix geometric means and semidefinite optimization" by Hamza
  Fawzi and James Saunderson (arXiv:1512.03401)
"""
function lieb_ando(
    A::Union{AbstractMatrix,Constant},
    B::Union{AbstractMatrix,Constant},
    K::Union{AbstractMatrix,Constant},
    t::Rational,
)
    if t < -1 || t > 2
        throw(DomainError(t, "t must be between -1 and 2"))
    end
    return real(LinearAlgebra.tr(K' * A^(1 - t) * K * B^t))
end

function lieb_ando(
    A::Union{AbstractMatrix,Constant},
    B::AbstractExpr,
    K::Union{AbstractMatrix,Constant},
    t::Rational,
)
    if t < -1 || t > 2
        throw(DomainError(t, "t must be between -1 and 2"))
    end
    KAK = K' * A^(1 - t) * K
    KAK = (KAK + KAK') / 2
    return trace_mpower(B, t, KAK)
end

function lieb_ando(
    A::AbstractExpr,
    B::Union{AbstractMatrix,Constant},
    K::Union{AbstractMatrix,Constant},
    t::Rational,
)
    if t < -1 || t > 2
        throw(DomainError(t, "t must be between -1 and 2"))
    end
    KBK = K * B^t * K'
    KBK = (KBK + KBK') / 2
    return trace_mpower(A, 1 - t, KBK)
end

function lieb_ando(
    A::AbstractExpr,
    B::AbstractExpr,
    K::Union{AbstractMatrix,Constant},
    t::Rational,
)
    n = size(A, 1)
    m = size(B, 1)
    Kvec = reshape(K', n * m, 1)
    KvKv = Kvec * Kvec'
    KvKv = (KvKv + KvKv') / 2
    Im = Matrix(1.0 * LinearAlgebra.I, m, m)
    In = Matrix(1.0 * LinearAlgebra.I, n, n)

    is_complex =
        sign(A) == ComplexSign() ||
        sign(B) == ComplexSign() ||
        sign(constant(K)) == ComplexSign()
    if is_complex
        T = HermitianSemidefinite(n * m)
    else
        T = Semidefinite(n * m)
    end

    if t >= 0 && t <= 1
        # Concave function
        add_constraint!(
            T,
            T in
            GeometricMeanHypoCone(kron(A, Im), kron(In, conj(B)), t, false),
        )
        return real(LinearAlgebra.tr(KvKv * T))
    elseif (t >= -1 && t <= 0) || (t >= 1 && t <= 2)
        # Convex function
        add_constraint!(
            T,
            Convex.GenericConstraint(
                (T, kron(A, Im), kron(In, conj(B))),
                GeometricMeanEpiConeSquare(t, size(T, 1)),
            ),
        )
        return real(LinearAlgebra.tr(KvKv * T))
    else
        throw(DomainError(t, "t must be between -1 and 2"))
    end
end
