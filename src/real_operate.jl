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
    return add_operation(
        tape,
        SparseAffineOperation(sparse(one(T) * I, d, d), v),
    )
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
        SparseAffineOperation(sparse(-one(T) * I, d, d), zeros(T, d)),
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

function real_operate(
    ::typeof(vcat),
    ::Type{T},
    tape::SparseTape{T},
    v::Vector{T},
) where {T<:Real}
    op = SparseAffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling
    Z = spzeros(T, n, m)
    A = vcat(op.matrix, Z)
    b = vcat(op.vector, v)
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
    Z = spzeros(T, n, m)
    A = vcat(Z, op.matrix)
    b = vcat(v, op.vector)
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
    v::Vector{T},
) where {T<:Real}
    return [sum(v)]
end

## Binary

# `add_operation`

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
    A,
    v::Vector{T},
) where {T<:Real}
    return A * v
end

function real_operate(
    ::typeof(add_operation),
    ::Type{T},
    x::Real,
    tape::SparseTape,
) where {T<:Real}
    d = MOI.output_dimension(tape)
    return add_operation(
        tape,
        SparseAffineOperation(sparse(T(x) * I, d, d), zeros(T, d)),
    )
end

# Here we have our two complex -> real functions

# `real`
#1. ComplexTape
function real_operate(::typeof(real), ::Type{T}, c::ComplexTape{T}) where {T}
    return real(c)
end

#2. SparseTape
function real_operate(::typeof(real), ::Type{T}, tape::SparseTape{T}) where {T}
    return tape
end

#3. ComplexStructOfVec
function real_operate(
    ::typeof(real),
    ::Type{T},
    v::ComplexStructOfVec{T},
) where {T}
    return real(v)
end

# `imag`
#1. ComplexTape
function real_operate(::typeof(imag), ::Type{T}, c::ComplexTape{T}) where {T}
    return imag(c)
end

#2. SparseTape
function real_operate(::typeof(imag), ::Type{T}, c::SparseTape{T}) where {T}
    return imag(c)
end

#3. ComplexStructOfVec
function real_operate(
    ::typeof(imag),
    ::Type{T},
    c::ComplexStructOfVec{T},
) where {T}
    return imag(c)
end
