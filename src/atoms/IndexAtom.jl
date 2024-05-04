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

    function IndexAtom(
        x::AbstractExpr,
        rows::AbstractArray,
        cols::AbstractArray,
    )
        return new((x,), (length(rows), length(cols)), rows, cols, nothing)
    end

    function IndexAtom(x::AbstractExpr, inds::AbstractArray)
        return new((x,), (length(inds), 1), nothing, nothing, inds)
    end
end

head(io::IO, ::IndexAtom) = print(io, "index")

Base.sign(x::IndexAtom) = sign(x.children[1])

monotonicity(::IndexAtom) = (Nondecreasing(),)

curvature(::IndexAtom) = ConstVexity()

function evaluate(x::IndexAtom)
    result = if x.inds === nothing
        getindex(evaluate(x.children[1]), x.rows, x.cols)
    else
        getindex(evaluate(x.children[1]), x.inds)
    end
    # `output` needed to convert a 1-element array to a scalar, see
    # https://github.com/jump-dev/Convex.jl/issues/447)
    return output(result)
end

function _index_real!(
    context::Context{T},
    obj_size::Tuple,
    obj_tape::Union{SparseTape{T},SPARSE_VECTOR{T}},
    x::IndexAtom,
) where {T}
    if x.inds === nothing
        linear_indices =
            LinearIndices(CartesianIndices(obj_size))[x.rows, x.cols]
    else
        linear_indices = collect(x.inds)
    end
    sz = length(linear_indices)

    # Here, we are in the real case, so `obj_tape` is either a `SparseTape{T}`, or a `SPARSE_VECTOR`.
    # In the latter case, we can handle it directly
    if obj_tape isa SPARSE_VECTOR
        return obj_tape[vec(linear_indices)]
    end
    # Ok, in this case we have actual work to do. We will construct an auxiliary variable `out`,
    # which we will return, and we will constrain it to the values we want.
    # This speeds up formulation since we reduce the problem size, and send what we have over to MOI already.
    out = Variable(sz)
    out_tape = conic_form!(context, out)
    for (i, I) in enumerate(linear_indices)
        # For each index, we constrain an element of `out` via ScalarAffineFunction to the indexed value.
        saf = to_saf(obj_tape, I)
        push!(saf.terms, MOI.ScalarAffineTerm(T(-1), out_tape.variables[i]))
        MOI.Utilities.normalize_and_add_constraint(
            context.model,
            saf,
            MOI.EqualTo(T(0)),
        )
    end
    return out_tape
end

function new_conic_form!(context::Context{T}, x::IndexAtom) where {T}
    input = x.children[1]
    if !iscomplex(x) # real case
        input_tape = conic_form!(context, input)
        return _index_real!(context, size(input), input_tape, x)
    else # complex case
        input_tape = conic_form!(context, input)
        re = _index_real!(context, size(input), real(input_tape), x)
        im = _index_real!(context, size(input), imag(input_tape), x)
        if re isa SPARSE_VECTOR
            @assert im isa SPARSE_VECTOR
            return ComplexStructOfVec(re, im)
        else
            return ComplexTape(re, im)
        end
    end
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
    return getindex(x, row:row, cols)
end

function Base.getindex(x::AbstractExpr, rows::AbstractVector{<:Real}, col::Real)
    return getindex(x, rows, col:col)
end

function Base.getindex(x::AbstractExpr, I::AbstractMatrix{Bool})
    return [xi for (xi, ii) in zip(x, I) if ii]
end

function Base.getindex(x::AbstractExpr, I::AbstractVector{Bool})
    return [xi for (xi, ii) in zip(x, I) if ii]
end

# All rows and columns
function Base.getindex(x::AbstractExpr, ::Colon, ::Colon)
    rows, cols = size(x)
    return getindex(x, 1:rows, 1:cols)
end

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
