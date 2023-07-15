import Base.vcat, Base.hcat, Base.hvcat
mutable struct HcatAtom <: AbstractExpr
    children::Tuple
    size::Tuple{Int,Int}

    function HcatAtom(args...)
        args = map(arg -> convert(AbstractExpr, arg), args)
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
        return new(children, (num_rows, num_cols))
    end
end

head(io::IO, ::HcatAtom) = print(io, "hcat")

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

function _conic_form!(context::Context{T}, x::HcatAtom) where {T}
    objectives = map(c -> conic_form!(context, c), children(x))
    # Suppose the child objectives for two children e1 (2 x 1) and e2 (2 x 2) look something like
    #  e1: x => 1 2 3
    #           4 5 6
    #      y => 2 4
    #           7 8
    #  e2: x => 1 1 1
    #           2 2 2
    #           3 3 3
    #           4 4 4
    # The objective of [e1 e2] will look like
    #            x => 1 2 3
    #                 4 5 6
    #                 1 1 1
    #                 2 2 2
    #                 3 3 3
    #                 4 4 4
    #            y => 2 4
    #                 7 8
    #                 0 0
    #                 0 0
    #                 0 0
    #                 0 0
    # builds the objective by aggregating a list of coefficients for each variable
    # from each child objective, and then vertically concatenating them
    return operate(vcat, T, sign(x), objectives...)
end
# TODO: fix piracy!

# * `Value` is not owned by Convex.jl
# * splatting creates zero-argument functions, which again are not owned by Convex.jl

function hcat(args::AbstractExprOrValue...)
    return HcatAtom(args...)
end
hcat(args::Value...) = Base.cat(args..., dims = Val(2))

# TODO: implement vertical concatenation in a more efficient way
function vcat(args::AbstractExprOrValue...)
    return transpose(HcatAtom(map(transpose, args)...))
end

vcat(args::Value...) = Base.cat(args..., dims = Val(1)) # Note: this makes general vcat slower for anyone using Convex...

function hvcat(rows::Tuple{Vararg{Int}}, args::AbstractExprOrValue...)
    nbr = length(rows)
    rs = Vector{Any}(undef, nbr)
    a = 1
    for i in 1:nbr
        rs[i] = HcatAtom(args[a:a-1+rows[i]]...)
        a += rows[i]
    end
    return vcat(rs...)
end
