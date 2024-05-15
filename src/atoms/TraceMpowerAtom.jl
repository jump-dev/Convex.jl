# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    TraceMpowerAtom(A::AbstractExpr, t::Rational, C::AbstractMatrix)

trace_mpower(A, t, C) returns LinearAlgebra.tr(C*A^t) where A and C are positive
definite matrices and C is constant and t ∈ [-1, 2].

When t ∈ [0,1], trace_mpower(A, t, C) is concave in A (for fixed
positive semidefinite matrix C) and convex for t ∈ [-1, 0) or (1, 2].

## Reference

Ported from CVXQUAD which is based on the paper: "Lieb's concavity theorem,
matrix geometric means and semidefinite optimization" by Hamza Fawzi and James
Saunderson (arXiv:1512.03401)
"""
mutable struct TraceMpowerAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    C::AbstractMatrix
    t::Rational

    function TraceMpowerAtom(A::AbstractExpr, t::Rational, C::AbstractMatrix)
        n = size(A, 1)
        if size(A) != size(C)
            throw(DimensionMismatch("A and C must be the same size"))
        elseif size(A) != (n, n)
            throw(DimensionMismatch("A and C must be square"))
        elseif norm(C - C') > 1e-6
            throw(DomainError(C, "C must be Hermitian"))
        elseif any(LinearAlgebra.eigvals(LinearAlgebra.Hermitian(C)) .< -1e-6)
            throw(DomainError(C, "C must be positive semidefinite"))
        elseif t < -1 || 2 < t
            throw(DomainError(t, "t must be in the range [-1, 2]"))
        end
        return new((A,), (1, 1), C, t)
    end
end

head(io::IO, ::TraceMpowerAtom) = print(io, "trace_mpower")

Base.sign(::TraceMpowerAtom) = NoSign()

monotonicity(::TraceMpowerAtom) = (NoMonotonicity(),)

function curvature(atom::TraceMpowerAtom)
    if 0 <= atom.t <= 1
        return ConcaveVexity()
    end
    return ConvexVexity()
end

function evaluate(atom::TraceMpowerAtom)
    return trace_mpower(evaluate(atom.children[1]), atom.t, atom.C)
end

function trace_mpower(
    A::AbstractExpr,
    t::Rational,
    C::Union{AbstractMatrix,Constant},
)
    # This `evaluate` is safe since it is not a `fix!`ed variable
    # (it must be a constant or matrix)
    return TraceMpowerAtom(A, t, evaluate(C))
end

function trace_mpower(
    A::Union{AbstractMatrix,Constant},
    t::Rational,
    C::Union{AbstractMatrix,Constant},
)
    return LinearAlgebra.tr(C * A^t)
end

function trace_mpower(
    A::Union{AbstractExpr,Value},
    t::Integer,
    C::Union{AbstractMatrix,Constant},
)
    return trace_mpower(A, t // 1, C)
end

function new_conic_form!(context::Context{T}, atom::TraceMpowerAtom) where {T}
    A = atom.children[1]
    n = size(A, 1)
    tmp = if sign(A) == ComplexSign()
        HermitianSemidefinite(n)
    else
        Semidefinite(n)
    end
    I = Matrix(one(T) * LinearAlgebra.I(n))
    if 0 <= atom.t <= 1
        add_constraint!(
            context,
            Constraint(
                (tmp, I, A),
                GeometricMeanHypoConeSquare(atom.t, n, false),
            ),
        )
    else
        add_constraint!(
            context,
            Convex.Constraint(
                (tmp, I, A),
                GeometricMeanEpiConeSquare(atom.t, size(A, 1)),
            ),
        )
    end
    # It's already a real mathematically, but Convex doesn't know it
    return conic_form!(context, real(LinearAlgebra.tr(atom.C * tmp)))
end
