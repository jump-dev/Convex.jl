# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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
        n = size(X, 1)
        if size(X) != size(C)
            throw(DimensionMismatch("X and C must be the same size"))
        elseif size(X) != (n, n)
            throw(DimensionMismatch("X and C must be square"))
        elseif norm(C - C') > 1e-6
            throw(DomainError(C, "C must be Hermitian"))
        elseif any(LinearAlgebra.eigvals(LinearAlgebra.Hermitian(C)) .< -1e-6)
            throw(DomainError(C, "C must be positive semidefinite"))
        end
        return new((X,), (1, 1), C, m, k)
    end
end

head(io::IO, ::TraceLogmAtom) = print(io, "trace_logm")

Base.sign(::TraceLogmAtom) = NoSign()

monotonicity(::TraceLogmAtom) = (NoMonotonicity(),)

curvature(::TraceLogmAtom) = ConcaveVexity()

evaluate(atom::TraceLogmAtom) = trace_logm(evaluate(atom.children[1]), atom.C)

function new_conic_form!(context::Context{T}, atom::TraceLogmAtom) where {T}
    X = atom.children[1]
    τ = if sign(X) == ComplexSign() || sign(constant(atom.C)) == ComplexSign()
        ComplexVariable(size(X))
    else
        Variable(size(X))
    end
    I = Matrix(one(T) * LinearAlgebra.I(size(X, 1)))
    set = RelativeEntropyEpiConeSquare(size(X, 1), atom.m, atom.k)
    add_constraint!(context, Constraint((τ, I, X), set))
    # It's already a real mathematically, but need to make it a real type.
    return conic_form!(context, real(-LinearAlgebra.tr(atom.C * τ)))
end
