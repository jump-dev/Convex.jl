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

function new_conic_form!(context::Context{T}, x::IndexAtom) where {T}
    obj = conic_form!(context, only(AbstractTrees.children(x)))
    m = length(x)
    n = length(x.children[1])
    if x.inds === nothing
        sz = length(x.cols) * length(x.rows)
        if !iscomplex(x) # only real case handled here, for now
            obj = x.children[1]
            obj_tape = conic_form!(context, obj)
            linear_indices =
                LinearIndices(CartesianIndices(size(obj)))[x.rows, x.cols]
            # Here, we are in the real case, so `obj_tape` is either a `SparseTape{T}`, or a `Vector{T}`.
            # In the latter case, we can handle it directly
            if obj_tape isa Vector{T}
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
                push!(
                    saf.terms,
                    MOI.ScalarAffineTerm(T(-1), out_tape.variables[i]),
                )
                MOI.add_constraint(context.model, saf, MOI.EqualTo(T(0)))
            end
            return out_tape
        else
            J = Vector{Int}(undef, sz)
            k = 1
            num_rows = x.children[1].size[1]
            for c in x.cols
                for r in x.rows
                    J[k] = num_rows * (convert(Int, c) - 1) + convert(Int, r)
                    k += 1
                end
            end
            index_matrix = create_sparse(T, collect(1:sz), J, one(T), m, n)
        end
    else
        index_matrix = create_sparse(
            T,
            collect(1:length(x.inds)),
            collect(x.inds),
            one(T),
            m,
            n,
        )
    end
    return operate(add_operation, T, sign(x), index_matrix, obj)
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
