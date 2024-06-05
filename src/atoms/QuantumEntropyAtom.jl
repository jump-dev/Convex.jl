# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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
    add_constraint!(context, Constraint((τ, X, I), set))
    # It's already a real mathematically, but need to make it a real type.
    return conic_form!(context, real(-LinearAlgebra.tr(τ)))
end
