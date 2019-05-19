#############################################################################
# relative_entropy.jl
# relative entropy (ie, sum_i( x_i log (x_i/y_i) ) of expressions x and y
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export relative_entropy, log_perspective
export sign, curvature, monotonicity, evaluate

# TODO: make this work for a *list* of inputs, rather than just for scalar/vector/matrix inputs

struct RelativeEntropyAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int, Int}

    function RelativeEntropyAtom(x::AbstractExpr, y::AbstractExpr)
        if sign(x) == ComplexSign() || sign(y) == ComplexSign()
            error("both arguments should be real but they are instead $(sign(x)) and $(sign(y))")
        else
            children = (x, y)
            return new(:entropy, hash(children), children, size(x))
        end
    end
end

function sign(x::RelativeEntropyAtom)
    return NoSign()
end

function monotonicity(x::RelativeEntropyAtom)
    return (NoMonotonicity(),NoMonotonicity())
end

function curvature(x::RelativeEntropyAtom)
    return ConvexVexity()
end

function evaluate(e::RelativeEntropyAtom)
    x = evaluate(e.children[1])
    y = evaluate(e.children[2])
    if any(isnan, y) return Inf end

    out = x.*log.(x./y)
    # fix value when x=0:
    # out will only be NaN if x=0, in which case the correct value is 0
    out[isnan.(out)] = 0
    return out
end

function conic_form!(e::RelativeEntropyAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, e)
        # transform to conic form:
        # x log x/y <= z
        # x log y/x >= -z
        # log y/x >= -z/x
        # y/x >= exp(-z/x)
        # y >= x exp(-z/x)
        # and cf the standard form for the exponential cone {(x,y,z): y*exp(x/y) <= z}
        z = Variable(e.size)
        x = e.children[1]
        y = e.children[2]
        objective = conic_form!(z, unique_conic_forms)
        for i=1:size(x,1)
            for j=1:size(x,2)
                conic_form!(ExpConstraint(-z[i,j], x[i,j], y[i,j]), unique_conic_forms)
            end
        end
        # need to constrain x>=0 and y>0.
        # x>=0 we get for free from the form of the exponential cone, so just add
        conic_form!(y>=0, unique_conic_forms) # nb we don't know how to ask for strict inequality
        cache_conic_form!(unique_conic_forms, e, objective)
    end
    return get_conic_form(unique_conic_forms, e)
end

relative_entropy(x::AbstractExpr, y::AbstractExpr) = sum(RelativeEntropyAtom(x, y))
# y*log(x/y)
log_perspective(x::AbstractExpr, y::AbstractExpr) = -relative_entropy(y, x)
