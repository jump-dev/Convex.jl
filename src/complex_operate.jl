# Here we deal with real -> complex and complex -> complex functions

# Earlier, we called `promote_to_complex` to ensure `Vector{T}` is promoted to `ComplexStructOfVec{T}`. Therefore we are left with three possible input types.

# Here we consume and produce `AllAllowedComplex` elements
const AllAllowedComplex{T} =
    Union{ComplexTape{T},ComplexStructOfVec{T},SparseTape{T}}

## Vararg

# `+`

function complex_operate(
    ::typeof(+),
    ::Type{T},
    args::AllAllowedComplex{T}...,
) where {T}
    re = real_operate(+, T, (real(a) for a in args)...)
    im = real_operate(+, T, (imag(a) for a in args)...)
    return ComplexTape(re, im)
end

# Same for `-`
function complex_operate(
    ::typeof(-),
    ::Type{T},
    args::AllAllowedComplex{T}...,
) where {T}
    re = real_operate(-, T, (real(a) for a in args)...)
    im = real_operate(-, T, (imag(a) for a in args)...)
    return ComplexTape(re, im)
end

# vcat

# 1-arg does nothing
function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    arg::AllAllowedComplex{T},
) where {T<:Real}
    return arg
end

# 2-arg
function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    x::AllAllowedComplex{T},
    y::AllAllowedComplex{T},
) where {T<:Real}
    re = real_operate(vcat, T, real(x), real(y))
    im = real_operate(vcat, T, imag(x), imag(y))
    return ComplexTape(re, im)
end

# 3-arg: reduce
function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    arg1::AllAllowedComplex{T},
    arg2::AllAllowedComplex{T},
    arg3::AllAllowedComplex{T},
    args::AllAllowedComplex{T}...,
) where {T<:Real}
    all_args = (arg1, arg2, arg3, args...)
    return foldl(
        (a, b) -> complex_operate(vcat, T, a, b),
        all_args,
    )::ComplexTape{T}
end

## Unary
# `conj`
#1. ComplexTape
function complex_operate(::typeof(conj), ::Type{T}, c::ComplexTape{T}) where {T}
    im = real_operate(-, T, imag(c))
    return ComplexTape(real(c), im)
end

#2. SparseTape - should not be possible, since it has real output
# function complex_operate(::typeof(conj), ::Type{T}, tape::SparseTape{T}) where {T}
# return tape
# end

#3. ComplexStructOfVec
function complex_operate(
    ::typeof(conj),
    ::Type{T},
    v::ComplexStructOfVec{T},
) where {T}
    return conj(v)
end

# `sum`

# 1. ComplexTape
function complex_operate(::typeof(sum), ::Type{T}, c::ComplexTape{T}) where {T}
    re = real_operate(sum, T, real(c))
    im = real_operate(sum, T, imag(c))
    return ComplexTape(re, im)
end

#2. SparseTape - should not be possible, since this has real output
# function complex_operate(::typeof(sum), ::Type{T}, c::SparseTape{T}) where {T}
# return real_operate(sum, T, c)
# end

# 3. ComplexStructOfVec
function complex_operate(
    ::typeof(sum),
    ::Type{T},
    c::ComplexStructOfVec{T},
) where {T}
    return ComplexStructOfVec([sum(real(c))], [sum(imag(c))])
end

## Binary

function complex_operate(::typeof(add_operation), ::Type{T}, A, c) where {T}
    re = real_operate(
        -,
        T,
        real_operate(add_operation, T, real(A), real(c)),
        real_operate(add_operation, T, imag(A), imag(c)),
    )

    im = real_operate(
        +,
        T,
        real_operate(add_operation, T, real(A), imag(c)),
        real_operate(add_operation, T, imag(A), real(c)),
    )
    if re isa Vector && im isa Vector
        return ComplexStructOfVec(re, im)
    else
        return ComplexTape(re, im)
    end
end
