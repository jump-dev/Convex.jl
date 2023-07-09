#############################################################################
# rationalnorm.jl
#
# Handles the k-norms for k > 1 where k is rational. Reduces the
# k-norm constraint to at most 2d ceil(log2(n + m)) + O(d) second
# order cone constraints. Here d is the dimension of the problem and k
# is converted to a rational with value k = n / m. The procedure for
# this is based on the paper "Second-order cone programming" by
# F. Alizadeh and D. Goldfarb, Mathematical Programming, Series B,
# 95:3-51, 2001, and is documented in the pdf available at
# https://github.com/jump-dev/Convex.jl/raw/master/docs/supplementary/rational_to_socp.pdf
#############################################################################

### k-norm for rational k

struct RationalNormAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    k::Rational{Int64}

    function RationalNormAtom(x::AbstractExpr, k::Rational{Int})
        children = (x,)
        k >= 1 || error("p-norms not defined for p < 1")
        return new(:rationalnorm, hash(children), children, (1, 1), k)
    end
end

function sign(x::RationalNormAtom)
    return Positive()
end

# The monotonicity
function monotonicity(x::RationalNormAtom)
    return (sign(x.children[1]) * Nondecreasing(),)
end

function curvature(x::RationalNormAtom)
    return ConvexVexity()
end

function evaluate(x::RationalNormAtom)
    return sum(abs.(evaluate(x.children[1])) .^ x.k)^(1 / x.k)
end

function conic_form!(context::Context{T}, x::RationalNormAtom) where {T}
    v = conic_form!(context, only(x.children))
    d = length(only(x.children)) + 1
    t = Variable()
    t_obj = conic_form!(context, t)
    f = operate(vcat, T, t_obj, v)
    MOI_add_constraint(context.model, f, MOI.NormCone(Float64(x.k), d))
    return t_obj
end

function rationalnorm(x::AbstractExpr, k::Rational{Int})
    if sign(x) == ComplexSign()
        row, col = size(x)
        if row == 1 || col == 1
            return RationalNormAtom(abs(x), k)
        end
    else
        return RationalNormAtom(x, k)
    end
end
