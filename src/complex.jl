struct ComplexTape{T}
    real_tape::SparseVAFTape{T}
    imag_tape::SparseVAFTape{T}

    function ComplexTape(re::SparseVAFTape{T}, im::SparseVAFTape{T}) where {T}
        MOI.output_dimension(re) == MOI.output_dimension(im) ||
            DimensionMismatch()
        return new{T}(re, im)
    end
end

MOI.output_dimension(c::ComplexTape) = MOI.output_dimension(c.real_tape) # same value as for imag
Base.real(c::ComplexTape) = c.real_tape
Base.real(tape::SparseVAFTape) = tape
Base.imag(c::ComplexTape) = c.imag_tape

const ComplexTapeOrVec =
    Union{<:ComplexTape,<:ComplexStructOfVec,<:SparseVAFTape,<:AbstractVector}
const TapeOrVec = Union{<:SparseVAFTape,<:AbstractVector}

# this is pretty ugly...
function complex_operate(
    ::typeof(+),
    ::Type{T},
    args::ComplexTapeOrVec...,
) where {T}
    re = real_operate(+, T, (real(a) for a in args)...)
    # `imag` not defined for SparseVAFTape
    imag_parts =
        (imag(a) for a in args if a isa ComplexTape || a isa ComplexStructOfVec)
    if isempty(imag_parts)
        error("Wrong dispatch hit?")
    else
        im = real_operate(+, T, imag_parts...)
        return ComplexTape(re, im)
    end
end

function complex_operate(::typeof(+), ::Type{T}, args::TapeOrVec...) where {T}
    return real_operate(+, T, args...)
end

function complex_operate(::typeof(conj), ::Type{T}, c::ComplexTape) where {T}
    im = real_operate(-, T, imag(c))
    return ComplexTape(real(c), im)
end

function complex_operate(
    ::typeof(conj),
    ::Type{T},
    tape::SparseVAFTape,
) where {T}
    return tape
end

function complex_operate(::typeof(real), ::Type{T}, c::ComplexTape) where {T}
    return real(c)
end

function complex_operate(
    ::typeof(real),
    ::Type{T},
    tape::SparseVAFTape,
) where {T}
    return tape
end

function complex_operate(::typeof(imag), ::Type{T}, c::ComplexTape) where {T}
    return imag(c)
end

function complex_operate(::typeof(-), ::Type{T}, c::ComplexTape) where {T}
    re = real_operate(-, T, real(c))
    im = real_operate(-, T, imag(c))
    return ComplexTape(re, im)
end

function complex_operate(
    ::typeof(*),
    ::Type{T},
    A::AbstractArray{<:Real},
    c::ComplexTape,
) where {T}
    re = real_operate(*, T, A, real(c))
    im = real_operate(*, T, A, imag(c))
    return ComplexTape(re, im)
end

function complex_operate(
    ::typeof(*),
    ::Type{T},
    x::Real,
    c::ComplexTape,
) where {T}
    re = real_operate(*, T, x, real(c))
    im = real_operate(*, T, x, imag(c))
    return ComplexTape(re, im)
end

function complex_operate(::typeof(sum), ::Type{T}, c::ComplexTape) where {T}
    re = real_operate(sum, T, real(c))
    im = real_operate(sum, T, imag(c))
    return ComplexTape(re, im)
end

const ComplexValue = Union{
    <:Complex,
    AbstractArray{<:Complex},
    <:ComplexConstant,
    <:ComplexStructOfVec,
}

function complex_operate(
    ::typeof(*),
    ::Type{T},
    z::ComplexValue,
    tape::SparseVAFTape,
) where {T}
    tape = collapse(tape)
    re = real_operate(*, T, real(z), tape)
    im = real_operate(*, T, imag(z), tape)
    return ComplexTape(re, im)
end

function collapse(tape::ComplexTape)
    re = collapse(real(tape))
    im = collapse(imag(tape))
    return ComplexTape(re, im)
end

function complex_operate(
    ::typeof(*),
    ::Type{T},
    z::ComplexValue,
    tape::ComplexTape,
) where {T}
    tape = collapse(tape)
    re1 = real_operate(*, T, real(z), real(tape))
    re2 = real_operate(*, T, imag(z), imag(tape))
    re2neg = real_operate(-, T, re2)
    re = real_operate(+, T, re1, re2neg)

    im1 = real_operate(*, T, imag(z), real(tape))
    im2 = real_operate(*, T, real(z), imag(tape))
    im = real_operate(+, T, im1, im2)

    return ComplexTape(re, im)
end

function MOI_add_constraint(
    model,
    f::ComplexTape,
    set::Union{MOI.Zeros,MOI.Nonnegatives,MOI.Nonpositives},
)
    re_inds = MOI_add_constraint(model, to_vaf(real(f)), set)
    im_inds = MOI_add_constraint(model, to_vaf(imag(f)), set)
    return (re_inds, im_inds)
end

function MOI_add_constraint(model, f::ComplexTape, set::MOI.NormOneCone)
    re_inds = MOI_add_constraint(model, to_vaf(real(f)), set)
    im_inds = MOI_add_constraint(model, to_vaf(imag(f)), set)
    return (re_inds, im_inds)
end

# vcat

# 1-arg does nothing
function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    tape::ComplexTape,
) where {T<:Real}
    return tape
end

function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    tape1::ComplexTape,
    tape2::ComplexTape,
) where {T<:Real}
    re = real_operate(vcat, T, real(tape1), real(tape2))
    im = real_operate(vcat, T, imag(tape1), imag(tape2))
    return ComplexTape(re, im)
end

function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    tape::ComplexTape,
    v::Union{<:AbstractVector,<:ComplexStructOfVec},
) where {T<:Real}
    re = real_operate(vcat, T, real(tape), real(v))
    im = real_operate(vcat, T, imag(tape), imag(v))
    return ComplexTape(re, im)
end

function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    v::Union{<:AbstractVector,<:ComplexStructOfVec},
    tape::ComplexTape,
) where {T<:Real}
    re = real_operate(vcat, T, real(v), real(tape))
    im = real_operate(vcat, T, imag(v), imag(tape))
    return ComplexTape(re, im)
end

function complex_operate(
    ::typeof(vcat),
    ::Type{T},
    arg1::ComplexTapeOrVec,
    arg2::ComplexTapeOrVec,
    arg3::ComplexTapeOrVec,
    args::Vararg{<:ComplexTapeOrVec},
) where {T<:Real}
    all_args = (arg1, arg2, arg3, args...)
    return foldl(
        (a, b) -> complex_operate(vcat, T, a, b),
        all_args,
    )::ComplexTape{T}
end
