import Base.getindex

const ArrayOrNothing = Union{AbstractArray,Nothing}

mutable struct IndexAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    rows::ArrayOrNothing
    cols::ArrayOrNothing
    inds::ArrayOrNothing

    function IndexAtom(
        x::AbstractExpr,
        rows::AbstractArray,
        cols::AbstractArray,
    )
        sz = (length(rows), length(cols))
        children = (x,)
        return new(children, sz, rows, cols, nothing)
    end

    function IndexAtom(x::AbstractExpr, inds::AbstractArray)
        sz = (length(inds), 1)
        children = (x,)
        return new(children, sz, nothing, nothing, inds)
    end
end

head(io::IO, ::IndexAtom) = print(io, "index")

## Type definition ends here

function sign(x::IndexAtom)
    return sign(x.children[1])
end

function monotonicity(x::IndexAtom)
    return (Nondecreasing(),)
end

function curvature(x::IndexAtom)
    return ConstVexity()
end

function evaluate(x::IndexAtom)
    result = if x.inds === nothing
        getindex(evaluate(x.children[1]), x.rows, x.cols)
    else
        getindex(evaluate(x.children[1]), x.inds)
    end
    # `output` needed to convert a 1-element array to a scalar (https://github.com/jump-dev/Convex.jl/issues/447)
    return output(result)
end

function _conic_form!(context::Context{T}, x::IndexAtom) where {T}
    obj = conic_form!(context, only(children(x)))

    m = length(x)
    n = length(x.children[1])

    if x.inds === nothing
        sz = length(x.cols) * length(x.rows)
        J = Array{Int}(undef, sz)
        k = 1

        num_rows = x.children[1].size[1]
        for c in x.cols
            for r in x.rows
                J[k] = num_rows * (convert(Int, c) - 1) + convert(Int, r)
                k += 1
            end
        end

        index_matrix = GBMatrix{T,T}(collect(1:sz), J, one(T), m, n)
    else
        index_matrix = GBMatrix{T,T}(collect(1:length(x.inds)), collect(x.inds), one(T), m, n)
    end

    return operate(add_operation, T, sign(x), index_matrix, obj)
end

## API Definition begins

function getindex(
    x::AbstractExpr,
    rows::AbstractVector{T},
    cols::AbstractVector{T},
) where {T<:Real}
    return IndexAtom(x, rows, cols)
end
getindex(x::AbstractExpr, inds::AbstractVector{<:Real}) = IndexAtom(x, inds)
getindex(x::AbstractExpr, ind::Real) = getindex(x, ind:ind)
getindex(x::AbstractExpr, row::Real, col::Real) = getindex(x, row:row, col:col)
function getindex(x::AbstractExpr, row::Real, cols::AbstractVector{<:Real})
    return getindex(x, row:row, cols)
end
function getindex(x::AbstractExpr, rows::AbstractVector{<:Real}, col::Real)
    return getindex(x, rows, col:col)
end
# XXX todo: speed test; there are lots of possible solutions for this
function getindex(x::AbstractExpr, I::AbstractMatrix{Bool})
    return [xi for (xi, ii) in zip(x, I) if ii]
end
function getindex(x::AbstractExpr, I::AbstractVector{Bool})
    return [xi for (xi, ii) in zip(x, I) if ii]
end
# Colon methods
# All rows and columns
function getindex(x::AbstractExpr, cln_r::Colon, cln_c::Colon)
    rows, cols = size(x)
    return getindex(x, 1:rows, 1:cols)
end
# All rows for this column(s)
getindex(x::AbstractExpr, cln_r::Colon, col) = getindex(x, 1:size(x)[1], col)
# All columns for this row(s)
getindex(x::AbstractExpr, row, cln_c::Colon) = getindex(x, row, 1:size(x)[2])

# Cartesian Index
getindex(x::AbstractExpr, c::CartesianIndex{N}) where {N} = x[Tuple(c)...]

## API Definition ends
