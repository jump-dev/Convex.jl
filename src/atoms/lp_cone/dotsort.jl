#############################################################################
# dotsort.jl
# dotsort(a,b) computes dot(sort(a), sort(b))
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

# This atom computes dot(sort(x), sort(w)), where w is constant
# for example, if w = [1 1 1 0 0 0 ... 0], it computes the sum of the 3 largest elements of x
mutable struct DotSortAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    w::Value

    function DotSortAtom(x::AbstractExpr, w::Value)
        if sign(x) == ComplexSign()
            error("Argument should be real instead it is $(sign(x))")
        else
            if !(length(w) == length(x))
                error("x and w must be the same size")
            end
            children = (x,)
            vecw = reshape(w, length(x))
            return new(children, (1, 1), vecw)
        end
    end
end

head(io::IO, ::DotSortAtom) = print(io, "dotsort")

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
        return (Nondecreasing(),)
    else
        return (NoMonotonicity(),)
    end
end

function curvature(x::DotSortAtom)
    return ConvexVexity()
end

function evaluate(x::DotSortAtom)
    return sum(
        sort(vec(evaluate(x.children[1])), rev = true) .*
        sort(vec(x.w), rev = true),
    )
end

function _conic_form!(context::Context{T}, x::DotSortAtom) where {T}
    y = only(x.children)
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
    add_constraint!(context, y * w' <= onesvec * nu' + mu * onesvec')
    objective = conic_form!(context, sum(mu) + sum(nu))
    return objective
end

dotsort(a::AbstractExpr, b::Value) = DotSortAtom(a, b)
dotsort(b::Value, a::AbstractExpr) = DotSortAtom(a, b)
