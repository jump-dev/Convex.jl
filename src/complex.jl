struct ComplexTape{T1 <: VAFTapes, T2 <: VAFTapes}
    real_tape::T1
    imag_tape::T2

    function ComplexTape(re::T1, im::T2) where {T1, T2}
        MOI.output_dimension(re) == MOI.output_dimension(im) || DimensionMismatch()
        return new{T1, T2}(re, im)
    end
end

MOI.output_dimension(c::ComplexTape) = MOI.output_dimension(c.real_tape) # same value as for imag
Base.real(c::ComplexTape) = c.real_tape
Base.real(tape::VAFTapes) = tape
Base.imag(c::ComplexTape) = c.imag_tape


const ComplexTapeOrVec = Union{<:ComplexTape, <:ComplexStructOfVec, <:VAFTapes, <:AbstractVector{<:Real}}
const TapeOrVec = Union{<:VAFTapes, <:AbstractVector{<:Real}}

# this is pretty ugly...
function operate(::typeof(+), ::Type{T}, args::ComplexTapeOrVec...) where {T}
    re = real_operate(+, T, (real(a) for a in args)...)
    # `imag` not defined for VAFTapes
    imag_parts =  (imag(a) for a in args if a isa ComplexTape || a isa ComplexStructOfVec)
    if isempty(imag_parts)
        error("Wrong dispatch hit?")
    else
        im = real_operate(+, T, imag_parts...)
        return ComplexTape(re, im)
    end
end

function operate(::typeof(+), ::Type{T}, args::TapeOrVec...) where {T}
    return real_operate(+, T, args...)
end

function operate(::typeof(conj), ::Type{T}, c::ComplexTape) where {T}
    im = operate(-, T, imag(c))
    return ComplexTape(real(c), im)
end

function operate(::typeof(conj), ::Type{T}, tape::VAFTapes) where {T}
    return tape
end

function operate(::typeof(real), ::Type{T}, c::ComplexTape) where {T}
    return real(c)
end

function operate(::typeof(real), ::Type{T}, tape::VAFTapes) where {T}
    return tape
end

function operate(::typeof(imag), ::Type{T}, c::ComplexTape) where {T}
    return imag(c)
end

function operate(::typeof(-), ::Type{T}, c::ComplexTape) where {T}
    re = operate(-, T, real(c))
    im = operate(-, T, imag(c))
    return ComplexTape(re, im)
end

function operate(::typeof(*), ::Type{T}, A::AbstractArray{<:Real}, c::ComplexTape) where {T}
    re = operate(*, T, A, real(c))
    im = operate(*, T, A, imag(c))
    return ComplexTape(re, im)
end

function operate(::typeof(*), ::Type{T}, x::Real, c::ComplexTape) where {T}
    re = operate(*, T, x, real(c))
    im = operate(*, T, x, imag(c))
    return ComplexTape(re, im)
end

function operate(::typeof(sum), ::Type{T}, c::ComplexTape) where {T}
    re = operate(sum, T, real(c))
    im = operate(sum, T, imag(c))
    return ComplexTape(re, im)
end


const ComplexValue = Union{<:Complex, AbstractArray{<:Complex}}

function operate(::typeof(*), ::Type{T}, z::ComplexValue, tape::VAFTapes) where {T}
    tape = collapse(tape)
    re = operate(*, T, real(z), tape)
    im = operate(*, T, imag(z), tape)
    return ComplexTape(re, im)
end

function collapse(tape::ComplexTape)
    re = collapse(real(tape))
    im = collapse(imag(tape))
    return ComplexTape(re, im)
end

function operate(::typeof(*), ::Type{T}, z::ComplexValue, tape::ComplexTape) where {T}
    tape = collapse(tape)
    re1 = operate(*, T, real(z), real(tape))
    re2 = operate(*, T, imag(z), imag(tape))
    re2neg = operate(-, T, re2)
    re = operate(+, T, re1, re2neg)
    
    im1 = operate(*, T, imag(z), real(tape))
    im2 = operate(*, T, real(z), imag(tape))
    im = operate(+, T, im1, im2)
    
    return ComplexTape(re, im)
end

function MOI_add_constraint(model, f::ComplexTape, set::Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives})
    re_inds = MOI_add_constraint(model, to_vaf(real(f)), set)
    im_inds = MOI_add_constraint(model, to_vaf(imag(f)), set)
    return (re_inds, im_inds)
end

function MOI_add_constraint(model, f::ComplexTape, set::MOI.NormOneCone)
    re_inds = MOI_add_constraint(model, to_vaf(real(f)), set)
    im_inds = MOI_add_constraint(model, to_vaf(imag(f)), set)
    return (re_inds, im_inds)
end
