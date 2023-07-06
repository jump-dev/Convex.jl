import Base.vcat, Base.hcat, Base.hvcat
struct HcatAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple
    size::Tuple{Int,Int}

    function HcatAtom(args::AbstractExpr...)
        num_rows = args[1].size[1]
        num_cols = 0
        for arg in args
            if arg.size[1] != num_rows
                error(
                    "Cannot horizontally stack expressions of varying number of rows",
                )
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
    # Is this right? Should it be `hcat`?
    return operate(vcat, T, objectives...)
end
# TODO: fix piracy!

# * `Value` is not owned by Convex.jl
# * splatting creates zero-argument functions, which again are not owned by Convex.jl
Base.hcat(args::AbstractExpr...) = HcatAtom(args...)

function Base.hcat(args::AbstractExprOrValue...)
    if all(Base.Fix2(isa, Value), args)
        return Base.cat(args..., dims = Val(2))
    end
    return HcatAtom(map(arg -> convert(AbstractExpr, arg), args)...)
end

# TODO: implement vertical concatenation in a more efficient way
vcat(args::AbstractExpr...) = transpose(HcatAtom(map(transpose, args)...))
function vcat(args::AbstractExprOrValue...)
    if all(Base.Fix2(isa, Value), args)
        return Base.cat(args..., dims = Val(1))
    end
    return transpose(
        HcatAtom(map(arg -> transpose(convert(AbstractExpr, arg)), args)...),
    )
end

function hvcat(rows::Tuple{Vararg{Int}}, args::AbstractExprOrValue...)
    nbr = length(rows)
    rs = Vector{Any}(undef, nbr)
    a = 1
    for i in 1:nbr
        rs[i] = hcat(args[a:a-1+rows[i]]...)
        a += rows[i]
    end
    return vcat(rs...)
end
