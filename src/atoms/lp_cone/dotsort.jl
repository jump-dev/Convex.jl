#############################################################################
# dotsort.jl
# dotsort(a,b) computes dot(sort(a), sort(b))
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export dotsort
export sign, curvature, monotonicity, evaluate

# This atom computes dot(sort(x), sort(w)), where w is constant
# for example, if w = [1 1 1 0 0 0 ... 0], it computes the sum of the 3 largest elements of x
struct DotSortAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    w::Value

    function DotSortAtom(x::AbstractExpr, w::Value)
        if sign(x) == ComplexSign()
            error("argument should be real instead it is $(sign(x))")
        else
            if !(length(w) == length(x))
                error("x and w must be the same size")
            end
            children = (x,)
            vecw = reshape(w, length(x))
            return new(:dotsort, hash((children, vecw)), children, (1,1), vecw)
        end
    end
end

function sign(x::DotSortAtom)
    if all(x.w .>= 0)
        return sign(x.children[1])
    elseif all(x.w .<= 0)
        return sign(x.children[1])
    else
        return NoSign()
    end
end

function monotonicity(x::DotSortAtom)
    if all(x.w .>= 0)
        return (Nondecreasing(), )
    else
        return (NoMonotonicity(), )
    end
end

function curvature(x::DotSortAtom)
    return ConvexVexity()
end

function evaluate(x::DotSortAtom)
    return sum(sort(vec(evaluate(x.children[1])), rev=true) .* sort(vec(x.w), rev=true))
end

function conic_form!(x::DotSortAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        y = x.children[1]
        w = x.w
        sy = size(y)
        if sy[1] > 1 && sy[2] > 1
            y = vec(y)
        end
        mu = Variable(size(y))
        nu = Variable(size(y))
        onesvec = ones(size(y))
        # given by the solution to
        # minimize sum(mu) + sum(nu)
        # subject to y*w' <= onesvec*nu' + mu*onesvec'
        objective = conic_form!(sum(mu) + sum(nu), unique_conic_forms)
        conic_form!(y*w' <= onesvec*nu' + mu*onesvec', unique_conic_forms)
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

dotsort(a::AbstractExpr, b::Value) = DotSortAtom(a, b)
dotsort(b::Value, a::AbstractExpr) = DotSortAtom(a, b)
