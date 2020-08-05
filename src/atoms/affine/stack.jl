import Base.vcat, Base.hcat, Base.hvcat
struct HcatAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple
    size::Tuple{Int, Int}

    function HcatAtom(args::AbstractExpr...)
        num_rows = args[1].size[1]
        num_cols = 0
        for arg in args
            if arg.size[1] != num_rows
                error("Cannot horizontally stack expressions of varying number of rows")
            end
            num_cols += arg.size[2]
        end
        children = tuple(args...)
        return new(:hcat, hash(children), children, (num_rows, num_cols))
    end
end

function sign(x::HcatAtom)
    return sum(map(sign, x.children))
end

function monotonicity(x::HcatAtom)
    return [Nondecreasing() for c in x.children]
end

function curvature(x::HcatAtom)
    return ConstVexity()
end

function evaluate(x::HcatAtom)
    return hcat(map(evaluate, x.children)...)
end

function template(x::HcatAtom, context::Context{T}) where {T}
    objectives = template.(children(x), Ref(context))
    operate(vcat, T, objectives...)
end

hcat(args::AbstractExpr...) = HcatAtom(args...)
hcat(args::AbstractExprOrValue...) = HcatAtom(map(arg -> convert(AbstractExpr, arg), args)...)
hcat(args::Value...) = Base.cat(args..., dims=Val(2))


# TODO: implement vertical concatenation in a more efficient way
vcat(args::AbstractExpr...) = transpose(HcatAtom(map(transpose, args)...))
vcat(args::AbstractExprOrValue...) = transpose(HcatAtom(map(arg -> transpose(convert(AbstractExpr, arg)), args)...))
vcat(args::Value...) = Base.cat(args..., dims=Val(1)) # Note: this makes general vcat slower for anyone using Convex...


function hvcat(rows::Tuple{Vararg{Int}}, args::AbstractExprOrValue...)
    nbr = length(rows)
    rs = Vector{Any}(undef, nbr)
    a = 1
    for i = 1:nbr
        rs[i] = hcat(args[a:a-1+rows[i]]...)
        a += rows[i]
    end
    return vcat(rs...)
end
