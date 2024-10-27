# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct IndexAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    rows::Union{AbstractArray,Nothing}
    cols::Union{AbstractArray,Nothing}
    inds::Union{AbstractArray,Nothing}
end

function IndexAtom(x::AbstractExpr, rows::AbstractArray, cols::AbstractArray)
    return IndexAtom((x,), (length(rows), length(cols)), rows, cols, nothing)
end

function IndexAtom(x::AbstractExpr, inds::AbstractArray)
    return IndexAtom((x,), (length(inds), 1), nothing, nothing, inds)
end

head(io::IO, ::IndexAtom) = print(io, "index")

Base.sign(x::IndexAtom) = sign(x.children[1])

monotonicity(::IndexAtom) = (Nondecreasing(),)

curvature(::IndexAtom) = ConstVexity()

function evaluate(x::IndexAtom)
    result = if x.inds === nothing
        # reshape to ensure we are respecting that a scalar row index
        # creates a column vector. We can't just check `length(x.rows)`
        # since that doesn't distinguish between `i` and `i:i`. But
        # we had that info when we set the size, so we will use it now.
        reshape(getindex(evaluate(x.children[1]), x.rows, x.cols), x.size)
    else
        getindex(evaluate(x.children[1]), x.inds)
    end
    # `output` needed to convert a 1-element array to a scalar, see
    # https://github.com/jump-dev/Convex.jl/issues/447)
    return output(result)
end

function _index(tape::SparseTape{T}, keep_rows::Vector{Int}) where {T}
    A = tape.operation.matrix
    indexed = A[keep_rows, :]
    af = SparseAffineOperation{T}(indexed, tape.operation.vector[keep_rows])
    return SparseTape{T}(af, tape.variables)
end

_index(tape::Vector, keep_rows::Vector{Int}) = tape[keep_rows]

function _index_real(
    obj_size::Tuple,
    obj_tape::Union{SparseTape,SPARSE_VECTOR},
    x::IndexAtom,
)
    if x.inds === nothing
        linear_indices = LinearIndices(CartesianIndices(obj_size))
        return _index(obj_tape, vec(linear_indices[x.rows, x.cols]))
    end
    return _index(obj_tape, vec(collect(x.inds)))
end

function new_conic_form!(context::Context{T}, x::IndexAtom) where {T}
    input = x.children[1]
    if !iscomplex(x) # real case
        input_tape = conic_form!(context, input)
        return _index_real(size(input), input_tape, x)
    end
    input_tape = conic_form!(context, input)
    re = _index_real(size(input), real(input_tape), x)
    im = _index_real(size(input), imag(input_tape), x)
    if re isa SPARSE_VECTOR
        @assert im isa SPARSE_VECTOR
        return ComplexStructOfVec(re, im)
    end
    return ComplexTape(re, im)
end

function Base.getindex(
    x::AbstractExpr,
    rows::AbstractVector{T},
    cols::AbstractVector{T},
) where {T<:Real}
    return IndexAtom(x, rows, cols)
end

function Base.getindex(x::AbstractExpr, inds::AbstractVector{<:Real})
    return IndexAtom(x, inds)
end

Base.getindex(x::AbstractExpr, ind::Real) = getindex(x, ind:ind)

function Base.getindex(x::AbstractExpr, row::Real, col::Real)
    return getindex(x, row:row, col:col)
end

function Base.getindex(x::AbstractExpr, row::Real, cols::AbstractVector{<:Real})
    # In this case, we must construct a column vector
    # https://github.com/jump-dev/Convex.jl/issues/509
    # Here we construct `getindex(x, row:row, cols)`
    # except with the size set to be that of a column vector
    return IndexAtom((x,), (length(cols), 1), row:row, cols, nothing)
end

function Base.getindex(x::AbstractExpr, rows::AbstractVector{<:Real}, col::Real)
    return getindex(x, rows, col:col)
end

function Base.getindex(x::AbstractExpr, I::Union{AbstractMatrix{Bool}, <:BitMatrix})
    return [xi for (xi, ii) in zip(x, I) if ii]
end

function Base.getindex(x::AbstractExpr, I::Union{<:AbstractVector{Bool}, <:BitVector})
    return [xi for (xi, ii) in zip(x, I) if ii]
end

# All rows and columns
Base.getindex(x::AbstractExpr, ::Colon, ::Colon) = x

# All rows for this column(s)
function Base.getindex(x::AbstractExpr, ::Colon, col)
    return getindex(x, 1:size(x, 1), col)
end

# All columns for this row(s)
function Base.getindex(x::AbstractExpr, row, ::Colon)
    return getindex(x, row, 1:size(x, 2))
end

# Cartesian Index
Base.getindex(x::AbstractExpr, c::CartesianIndex{N}) where {N} = x[Tuple(c)...]
