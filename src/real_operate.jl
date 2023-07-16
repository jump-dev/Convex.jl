# Here we cover both real -> real and complex -> real
# First we cover real -> real, then complex -> real

# Only two types allowed here for real -> real
const AllAllowedReal{T} = Union{SparseTape{T},SPARSE_VECTOR{T}}

## Vararg

# `+`

# (Vector, Vector)
function real_operate(
    ::typeof(+),
    ::Type{T},
    v1::SPARSE_VECTOR{T},
    v2::SPARSE_VECTOR{T},
) where {T<:Real}
    return v1 + v2
end

# (Vector, SparseTape)
function real_operate(
    ::typeof(+),
    ::Type{T},
    v::SPARSE_VECTOR{T},
    tape::SparseTape{T},
) where {T<:Real}
    return real_operate(+, T, tape, v)
end

# (SparseTape, Vector)
function real_operate(
    ::typeof(+),
    ::Type{T},
    tape::SparseTape{T},
    v::SPARSE_VECTOR{T},
) where {T<:Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation(spidentity(T, d), v))
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
        SparseAffineOperation(-spidentity(T, d), spzeros(T, d)),
    )
end

function real_operate(
    ::typeof(-),
    ::Type{T},
    v::SPARSE_VECTOR{T},
) where {T<:Real}
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

function blockdiag(xs::GBMatrix{T,T}...) where {T}
    N = length(xs)
    entries = Matrix{GBMatrix{T,T}}(undef, N, N)
    heights = size.(xs, 1)
    for (i, x) in enumerate(xs)
        entries[i, i] = x
        m = size(x,2)
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
    v::SPARSE_VECTOR{T},
) where {T<:Real}
    op = SparseAffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling
    b = vcat(op.vector, v)
    # Workaround SparseSuiteGraphBLAS bug with vcat
    # where vcatting two (1,1) GBMatrix yields GBVector
    if op.matrix isa GBMatrix
        A = GBMatrix{T,T}(size(op.matrix, 1) + n, m)
        A[1:size(op.matrix, 1), :] = op.matrix
    else
        Z = spzeros(T, n, m)
        A = vcat(op.matrix, Z)
    end
    return SparseTape([SparseAffineOperation(A, b)], tape.variables)
end

function real_operate(
    ::typeof(vcat),
    ::Type{T},
    v::SPARSE_VECTOR{T},
    tape::SparseTape{T},
) where {T<:Real}
    op = SparseAffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling
    b = vcat(v, op.vector)

    # Workaround SparseSuiteGraphBLAS bug with vcat
    # where vcatting two (1,1) GBMatrix yields GBVector
    if op.matrix isa GBMatrix
        A = GBMatrix{T,T}(n + size(op.matrix, 1), m)
        A[(n+1):end, :] = op.matrix
    else
        Z = spzeros(T, n, m)
        A = vcat(Z, op.matrix)
    end
    return SparseTape([SparseAffineOperation(A, b)], tape.variables)
end

function real_operate(
    ::typeof(vcat),
    ::Type{T},
    v1::SPARSE_VECTOR{T},
    v2::SPARSE_VECTOR{T},
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

function real_operate(
    ::typeof(sum),
    ::Type{T},
    v::SPARSE_VECTOR{T},
) where {T<:Real}
    return [sum(v)]
end

## Binary

# `add_operation`
# Here the left-side argument may be a `SparseMatrixCSC{T}` or a `T`
# and the right-argument is either a `SparseTape{T}` or `SPARSE_VECTOR{T}`

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
    v::SPARSE_VECTOR{T},
) where {T<:Real}
    return SPARSE_VECTOR{T}(A * v)
end

function real_operate(
    ::typeof(add_operation),
    ::Type{T},
    x::Real,
    tape::SparseTape{T},
) where {T<:Real}
    d = MOI.output_dimension(tape)
    return add_operation(
        tape,
        SparseAffineOperation(spidentity(T, d), spzeros(T, d)),
    )
end

function real_operate(
    ::typeof(add_operation),
    ::Type{T},
    x::Real,
    v::SPARSE_VECTOR{T},
) where {T<:Real}
    return SPARSE_VECTOR{T}(real_convert(T, x) * v)
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
