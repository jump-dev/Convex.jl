# Here we consume and produce `AllAllowedComplex` elements
const AllAllowedComplex = Union{<:ComplexTape,<:ComplexStructOfVec,<:SparseTape}

## Vararg

# `+`

# We will define the general case (`AllAllowedComplex`), and then specific all-real case.
function complex_operate(
    ::typeof(+),
    ::Type{T},
    args::AllAllowedComplex...,
) where {T}
    re = real_operate(+, T, (real(a) for a in args)...)
    # `imag` not defined for SparseTape
    imag_parts =
        (imag(a) for a in args if a isa ComplexTape || a isa ComplexStructOfVec)
    if isempty(imag_parts)
        error("Should not be possible; wrong dispatch hit?")
    else
        im = real_operate(+, T, imag_parts...)
        return ComplexTape(re, im)
    end
end

# All real
function complex_operate(::typeof(+), ::Type{T}, args::SparseTape...) where {T}
    return real_operate(+, T, args...)
end

# Same strat for `-`
# `AllAllowedComplex`
function complex_operate(
    ::typeof(-),
    ::Type{T},
    args::AllAllowedComplex...,
) where {T}
    re = real_operate(-, T, (real(a) for a in args)...)
    # `imag` not defined for SparseTape
    imag_parts =
        (imag(a) for a in args if a isa ComplexTape || a isa ComplexStructOfVec)
    if isempty(imag_parts)
        error("Should not be possible; wrong dispatch hit?")
    else
        im = real_operate(-, T, imag_parts...)
        return ComplexTape(re, im)
    end
end

# All real
function complex_operate(::typeof(-), ::Type{T}, args::SparseTape...) where {T}
    return real_operate(-, T, args...)
end

# vcat

# 1-arg does nothing
function complex_operate(::typeof(vcat), ::Type{T}, arg) where {T<:Real}
    return arg
end

# 2-arg
function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    v::AllAllowedComplex,
    tape::AllAllowedComplex,
) where {T<:Real}
    re = real_operate(vcat, T, real(v), real(tape))
    im = real_operate(vcat, T, imag(v), imag(tape))
    return ComplexTape(re, im)
end

# 3-arg: reduce
function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    arg1::AllAllowedComplex,
    arg2::AllAllowedComplex,
    arg3::AllAllowedComplex,
    args::Vararg{<:AllAllowedComplex},
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
function complex_operate(::typeof(conj), ::Type{T}, c::ComplexTape) where {T}
    im = real_operate(-, T, imag(c))
    return ComplexTape(real(c), im)
end

#2. SparseTape
function complex_operate(::typeof(conj), ::Type{T}, tape::SparseTape) where {T}
    return tape
end

#3. ComplexStructOfVec
function complex_operate(
    ::typeof(conj),
    ::Type{T},
    v::ComplexStructOfVec,
) where {T}
    return conj(v)
end

# `sum`

# 1. ComplexTape
function complex_operate(::typeof(sum), ::Type{T}, c::ComplexTape) where {T}
    re = real_operate(sum, T, real(c))
    im = real_operate(sum, T, imag(c))
    return ComplexTape(re, im)
end

#2. SparseTape
function complex_operate(::typeof(sum), ::Type{T}, c::SparseTape) where {T}
    return real_operate(sum, T, c)
end

# 3. ComplexStructOfVec
function complex_operate(
    ::typeof(sum),
    ::Type{T},
    c::ComplexStructOfVec,
) where {T}
    return ComplexStructOfVec(sum(real(c)), sum(imag(c)))
end

## Binary

# Here we treat `*` as binary, and we don't do quadratic stuff
# so either one is ComplexStructOfVec or the other argument is, or both.
function complex_operate(::typeof(*), ::Type{T}, A, c) where {T}
    re = real_operate(
        -,
        T,
        real_operate(*, T, real(A), real(c)),
        real_operate(*, T, imag(A), imag(c)),
    )

    im = real_operate(
        +,
        T,
        real_operate(*, T, real(A), imag(c)),
        real_operate(*, T, imag(A), real(c)),
    )
    return ComplexTape(re, im)
end
