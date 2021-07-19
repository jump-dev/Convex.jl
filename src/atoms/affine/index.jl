import Base.getindex

const ArrayOrNothing = Union{AbstractArray, Nothing}

struct IndexAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    rows::ArrayOrNothing
    cols::ArrayOrNothing
    inds::ArrayOrNothing

    function IndexAtom(x::AbstractExpr, rows::AbstractArray, cols::AbstractArray)
        sz = (length(rows), length(cols))
        children = (x,)
        return new(:index, hash((children, rows, cols, nothing)), children, sz, rows, cols, nothing)
    end

    function IndexAtom(x::AbstractExpr, inds::AbstractArray)
        sz = (length(inds), 1)
        children = (x,)
        return new(:index, hash((children, nothing, nothing, inds)), children, sz, nothing, nothing, inds)
    end
end

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
    if x.inds === nothing
        return output(getindex(evaluate(x.children[1]), x.rows, x.cols)) # output() needed to convert a 1-element vector to a scalar, see the link to the Issue below.
    else
        return output(getindex(evaluate(x.children[1]), x.inds)) # output() needed to convert a 1-element vector to a scalar. In response to https://github.com/jump-dev/Convex.jl/issues/447.
    end
end

function conic_form!(x::IndexAtom, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, x)
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

            index_matrix = sparse(1:sz, J, 1.0, m, n)
        else
            index_matrix = sparse(1:length(x.inds), x.inds, 1.0, m, n)
        end
        objective = conic_form!(x.children[1], unique_conic_forms)
        objective = index_matrix * objective
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

## API Definition begins

getindex(x::AbstractExpr, rows::AbstractVector{T}, cols::AbstractVector{T}) where {T<:Real} = IndexAtom(x, rows, cols)
getindex(x::AbstractExpr, inds::AbstractVector{<:Real}) = IndexAtom(x, inds)
getindex(x::AbstractExpr, ind::Real) = getindex(x, ind:ind)
getindex(x::AbstractExpr, row::Real, col::Real) = getindex(x, row:row, col:col)
getindex(x::AbstractExpr, row::Real, cols::AbstractVector{<:Real}) = getindex(x, row:row, cols)
getindex(x::AbstractExpr, rows::AbstractVector{<:Real}, col::Real) = getindex(x, rows, col:col)
# XXX todo: speed test; there are lots of possible solutions for this
function getindex(x::AbstractExpr, I::AbstractMatrix{Bool})
    return [xi for (xi,ii) in zip(x,I) if ii]
end
function getindex(x::AbstractExpr, I::AbstractVector{Bool})
    return [xi for (xi,ii) in zip(x,I) if ii]
end
# Colon methods
# All rows and columns
function getindex(x::AbstractExpr, cln_r::Colon, cln_c::Colon)
    rows, cols = size(x)
    getindex(x, 1:rows, 1:cols)
end
# All rows for this column(s)
getindex(x::AbstractExpr, cln_r::Colon, col) = getindex(x, 1:size(x)[1], col)
# All columns for this row(s)
getindex(x::AbstractExpr, row, cln_c::Colon) = getindex(x, row, 1:size(x)[2])


## API Definition ends
