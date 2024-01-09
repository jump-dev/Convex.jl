#############################################################################
# nuclearnorm.jl
# Handles nuclear norm (the sum of the singular values of a matrix),
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

### Nuclear norm

mutable struct NuclearNormAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function NuclearNormAtom(x::AbstractExpr)
        children = (x,)
        return new(children, (1, 1))
    end
end

head(io::IO, ::NuclearNormAtom) = print(io, "nuclearnorm")

function Base.sign(x::NuclearNormAtom)
    return Positive()
end

# The monotonicity
function monotonicity(x::NuclearNormAtom)
    return (NoMonotonicity(),)
end

function curvature(x::NuclearNormAtom)
    return ConvexVexity()
end

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
        # I'm not sure how to use MOI's `NormNuclearCone` in this case, so we'll just do the extended formulation as an SDP ourselves:
        #   minimize (LinearAlgebra.tr(U) + LinearAlgebra.tr(V))/2
        #   subject to
        #            [U A; A' V] ⪰ 0
        # see eg Recht, Fazel, Parillo 2008 "Guaranteed Minimum-Rank Solutions of Linear Matrix Equations via Nuclear Norm Minimization"
        # http://arxiv.org/pdf/0706.4138v1.pdf
        m, n = size(A)
        U = ComplexVariable(m, m)
        V = ComplexVariable(n, n)
        p = minimize(
            real(LinearAlgebra.tr(U) + LinearAlgebra.tr(V)) / 2,
            [U A; A' V] ⪰ 0,
        )
        return conic_form!(context, p)
    else
        t = conic_form!(context, Variable())
        f = operate(vcat, T, sign(x), t, conic_form!(context, A))
        m, n = size(A)
        MOI_add_constraint(context.model, f, MOI.NormNuclearCone(m, n))
        return t
    end
end
