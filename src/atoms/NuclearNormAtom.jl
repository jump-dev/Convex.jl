# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct NuclearNormAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    NuclearNormAtom(x::AbstractExpr) = new((x,), (1, 1))
end

head(io::IO, ::NuclearNormAtom) = print(io, "nuclearnorm")

Base.sign(::NuclearNormAtom) = Positive()

monotonicity(::NuclearNormAtom) = (NoMonotonicity(),)

curvature(::NuclearNormAtom) = ConvexVexity()

function evaluate(x::NuclearNormAtom)
    return sum(LinearAlgebra.svdvals(evaluate(x.children[1])))
end

nuclearnorm(x::AbstractExpr) = NuclearNormAtom(x)

# The complex case is example 1.20 of Watrous' "The Theory of Quantum Information"
# (the operator A is negated but this doesn't affect the norm)
# https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf
function new_conic_form!(context::Context{T}, x::NuclearNormAtom) where {T}
    A = only(AbstractTrees.children(x))
    if iscomplex(sign(A))
        # I'm not sure how to use MOI's `NormNuclearCone` in this case, so we'll
        # just do the extended formulation as an SDP ourselves:
        #   minimize (LinearAlgebra.tr(U) + LinearAlgebra.tr(V))/2
        #   subject to
        #            [U A; A' V] ⪰ 0
        # see eg Recht, Fazel, Parillo 2008 "Guaranteed Minimum-Rank Solutions
        # of Linear Matrix Equations via Nuclear Norm Minimization"
        # http://arxiv.org/pdf/0706.4138v1.pdf
        m, n = size(A)
        U = ComplexVariable(m, m)
        V = ComplexVariable(n, n)
        p = minimize(
            real(LinearAlgebra.tr(U) + LinearAlgebra.tr(V)) / 2,
            [U A; A' V] ⪰ 0,
        )
        return conic_form!(context, p)
    end
    t = conic_form!(context, Variable())
    f = operate(vcat, T, sign(x), t, conic_form!(context, A))
    m, n = size(A)
    MOI_add_constraint(context.model, f, MOI.NormNuclearCone(m, n))
    return t
end
