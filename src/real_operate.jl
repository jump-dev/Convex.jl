# Here we cover both real -> real and complex -> real
# First we cover real -> real, then complex -> real

# Only two types allowed here for real -> real
const AllAllowedReal{T} = Union{SparseTape{T},Vector{T}}

## Vararg

# `+`

# (Vector, Vector)
function real_operate(
    ::typeof(+),
    ::Type{T},
    v1::Vector{T},
    v2::Vector{T},
) where {T<:Real}
    return v1 + v2
end

# (Vector, SparseTape)
function real_operate(
    ::typeof(+),
    ::Type{T},
    v::Vector{T},
    tape::SparseTape{T},
) where {T<:Real}
    return real_operate(+, T, tape, v)
end

# (SparseTape, Vector)
function real_operate(
    ::typeof(+),
    ::Type{T},
    tape::SparseTape{T},
    v::Vector{T},
) where {T<:Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation(gbidentity(T, d), v))
end

# (SparseTape, SparseTape)
function real_operate(
    ::typeof(+),
    ::Type{T},
    tape1::SparseTape{T},
    tape2::SparseTape{T},
) where {T<:Real}
    @assert MOI.output_dimension(tape1) == MOI.output_dimension(tape2)
    op1 = SparseAffineOperation(tape1)
    op2 = SparseAffineOperation(tape2)

    if tape1.variables == tape2.variables
        op = SparseAffineOperation(
            op1.matrix + op2.matrix,
            op1.vector + op2.vector,
        )
        return SparseTape([op], tape1.variables)
    else
        mat = hcat(op1.matrix, op2.matrix)
        vec = op1.vector + op2.vector
        op = SparseAffineOperation(mat, vec)
        return SparseTape([op], vcat(tape1.variables, tape2.variables))
    end
end

# Reduce to 2-arg
function real_operate(
    ::typeof(+),
    ::Type{T},
    arg1::AllAllowedReal{T},
    arg2::AllAllowedReal{T},
    arg3::AllAllowedReal{T},
    args::AllAllowedReal{T}...,
) where {T<:Real}
    all_args = (arg1, arg2, arg3, args...)
    vec_args = (a for a in all_args if a isa Vector)
    tape_args = (a for a in all_args if a isa SparseTape)
    if isempty(tape_args)
        return sum(vec_args)
    else
        tape = foldl((a, b) -> real_operate(+, T, a, b), tape_args)
    end
    if isempty(vec_args)
        return tape
    else
        v = sum(vec_args)
        return real_operate(+, T, tape, v)
    end
end

# `-`

# Unary `-`
function real_operate(
    ::typeof(-),
    ::Type{T},
    tape::SparseTape{T},
) where {T<:Real}
    d = MOI.output_dimension(tape)

    return add_operation(
        tape,
        SparseAffineOperation(-gbidentity(T, d), GBVector{T,T}(d)),
    )
end

function real_operate(::typeof(-), ::Type{T}, v::Vector{T}) where {T<:Real}
    return -v
end

# 2+ args: reduce to unary - and +
function real_operate(
    ::typeof(-),
    ::Type{T},
    x::AllAllowedReal{T},
    ys::AllAllowedReal{T}...,
) where {T<:Real}
    mys = (real_operate(-, T, y) for y in ys)
    return real_operate(+, T, x, mys...)
end

# `vcat`

# 1-arg does nothing
function real_operate(
    ::typeof(vcat),
    ::Type{T},
    v::AllAllowedReal{T},
) where {T<:Real}
    return v
end

# https://github.com/JuliaSparse/SuiteSparseGraphBLAS.jl/blob/ce89efca856aa1afe65c5913423eb87ac4239f87/src/operations/concat.jl#L14C1-L35C1
using SuiteSparseGraphBLAS: AbstractGBVector, AbstractGBMatrix, cat!
function _vcat(
    tiles::Union{AbstractGBMatrix{T,F},AbstractGBVector{T,F}}...,
) where {T,F}
    return cat(collect(tiles))

    fills = getproperty.(tiles, :fill)
    if F <: Union{Nothing,Missing}
        fill = F()
    elseif all(y -> y == fills[1], fills)
        fill = fills[1]
    else
        fill = zero(F)
    end
    ncols = size(tiles[1], 2)
    nrows = sum(size.(tiles, 1))
    types = eltype.(tiles)
    t = types[1]
    for type in types[2:end]
        t = promote_type(t, type)
    end
    # sz = (tiles isa AbstractArray && ncols == 1) ? (nrows,) : (nrows, ncols)

    if ndims(tiles[1]) == 1
        @assert ncols == 1
        @assert all(==(1), ndims.(tiles))
        sz = (nrows,)
        vec = true
    else
        sz = (nrows, ncols)
        vec = false
    end
    C = similar(tiles[1], t, sz; fill) # TODO: FIXME, we want this to use promotion, but it's complicated.

    # if !vec
    nr = 1
    for x in tiles
        # @show size(x)
        C[nr:(nr+size(x, 1)-1)] = x
        nr += size(x, 1)
    end
    # else
    #     nr = 1
    #     for x in tiles
    #         @show size(x)
    #         @show nr:(nr+size(x, 1))
    #         @show size(C)
    #         C[nr:(nr+size(x, 1)-1)] = x
    #         nr += size(x, 1)
    #     end
    # end
    return C
end

# function Base.setindex!(
#     u::AbstractGBVector,
#     x,
#     I::Union{Vector,UnitRange,StepRange,Colon},
#     ::Colon;
#     mask = nothing,
#     accum = nothing,
#     desc = nothing,
# )
#     return SuiteSparseGraphBLAS.subassign!(u, x, I; mask, accum, desc)
# end

# we do all pairs of `SparseTape` and `Vector`, and then do 3+ arguments by iterating
function real_operate(
    ::typeof(vcat),
    ::Type{T},
    tape1::SparseTape{T},
    tape2::SparseTape{T},
) where {T<:Real}
    op1 = SparseAffineOperation(tape1)
    op2 = SparseAffineOperation(tape2)
    A = blockdiag(op1.matrix, op2.matrix)
    b = vcat(op1.vector, op2.vector)
    x = vcat(tape1.variables, tape2.variables)
    return SparseTape([SparseAffineOperation(A, b)], x)
end

blockdiag(x, y) = Base.blockdiag(x, y)

function blockdiag(xs::GBMatrix{T}...) where {T}
    # nrows = sum(size.(xs, 1))
    # ncols = sum(size.(xs, 2))
    # z = GBMatrix{T}(nrows, ncols)
    # n = m = 1
    # for x in xs
    #     z[n:(n+size(x, 1)-1), m:(m+size(x,2)-1)] = x
    # end
    # return z

    N = length(xs)

    entries = Matrix{GBMatrix{T,T}}(undef, N, N)

    heights = size.(xs, 1)
    for (i, x) in enumerate(xs)
        entries[i, i] = x
        n, m = size(x)
        # @assert n == m
        for j in 1:(i-1)
            entries[j, i] = GBMatrix{T,T}(heights[j], m)
        end
        for j in (i+1):lastindex(entries, 1)
            entries[j, i] = GBMatrix{T,T}(heights[j], m)
        end

    end

    return cat(entries)
end

function real_operate(
    ::typeof(vcat),
    ::Type{T},
    tape::SparseTape{T},
    v::Vector{T},
) where {T<:Real}
    op = SparseAffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling
    Z = GBMatrix{T,T}(n, m)
    A = vcat(op.matrix, Z)
    b = vcat(op.vector, GBVector{T,T}(v))
    return SparseTape([SparseAffineOperation(A, b)], tape.variables)
end

function real_operate(
    ::typeof(vcat),
    ::Type{T},
    v::Vector{T},
    tape::SparseTape{T},
) where {T<:Real}
    op = SparseAffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling
    Z = GBMatrix{T,T}(n, m)
    A = vcat(Z, op.matrix)
    b = _vcat(GBVector{T,T}(v), op.vector)
    return SparseTape([SparseAffineOperation(A, b)], tape.variables)
end

function real_operate(
    ::typeof(vcat),
    ::Type{T},
    v1::Vector{T},
    v2::Vector{T},
) where {T<:Real}
    return vcat(v1, v2)
end

function real_operate(
    ::typeof(vcat),
    ::Type{T},
    arg1::AllAllowedReal{T},
    arg2::AllAllowedReal{T},
    arg3::AllAllowedReal{T},
    args::AllAllowedReal{T}...,
) where {T<:Real}
    all_args = (arg1, arg2, arg3, args...)
    return foldl((a, b) -> real_operate(vcat, T, a, b), all_args)::SparseTape{T}
end

function real_operate(
    ::typeof(vcat),
    ::Type{T},
    arg1::SparseTape{T},
    arg2::SparseTape{T},
    arg3::SparseTape{T},
    args::SparseTape{T}...,
) where {T<:Real}
    all_args = (arg1, arg2, arg3, args...)
    ops = SparseAffineOperation.(all_args)
    A = blockdiag((op.matrix for op in ops)...)
    b = vcat((op.vector for op in ops)...)
    x = vcat((arg.variables for arg in all_args)...)
    return SparseTape([SparseAffineOperation(A, b)], x)
end

## Unary

# `sum`
function real_operate(
    ::typeof(sum),
    ::Type{T},
    tape::SparseTape{T},
) where {T<:Real}
    d = MOI.output_dimension(tape)
    # doesn't seem ideal for a sparse representation...
    A = ones(T, 1, d)
    return add_operation(tape, SparseAffineOperation(A, zeros(T, size(A, 1))))
end

function real_operate(::typeof(sum), ::Type{T}, v::Vector{T}) where {T<:Real}
    return [sum(v)]
end

## Binary

# `add_operation`
# Here the left-side argument may be a `SparseMatrixCSC{T}` or a `T`
# and the right-argument is either a `SparseTape{T}` or `Vector{T}`

function real_operate(
    ::typeof(add_operation),
    ::Type{T},
    A::AbstractMatrix,
    tape::SparseTape{T},
) where {T<:Real}
    return add_operation(tape, SparseAffineOperation(A, zeros(T, size(A, 1))))
end

function real_operate(
    ::typeof(add_operation),
    ::Type{T},
    A::AbstractMatrix,
    v::Vector{T},
) where {T<:Real}
    return Vector{T}(A * v)
end

gbidentity(T, d) = GBMatrix{T,T}(Diagonal(ones(T, d)))

function real_operate(
    ::typeof(add_operation),
    ::Type{T},
    x::Real,
    tape::SparseTape{T},
) where {T<:Real}
    d = MOI.output_dimension(tape)
    return add_operation(
        tape,
        SparseAffineOperation(gbidentity(T, d), GBVector{T,T}(d)),
    )
end

function real_operate(
    ::typeof(add_operation),
    ::Type{T},
    x::Real,
    v::Vector{T},
) where {T<:Real}
    return Vector{T}(real_convert(T, x) * v)
end

# Here we have our two complex -> real functions
# These are allowed these inputs:
const ComplexToRealInputs{T} =
    Union{ComplexTape{T},SparseTape{T},ComplexStructOfVec{T}}

# `real`
function real_operate(
    ::typeof(real),
    ::Type{T},
    c::ComplexToRealInputs{T},
) where {T}
    return real(c)
end

# `imag`
function real_operate(
    ::typeof(imag),
    ::Type{T},
    c::ComplexToRealInputs{T},
) where {T}
    return imag(c)
end
