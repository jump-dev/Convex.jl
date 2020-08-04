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
Base.imag(c::ComplexTape) = c.imag_tape

function MOIU.operate(::typeof(+), ::Type{T}, args::Union{<:ComplexTape, <:ComplexStructOfVec, <:VAFTapes, <:AbstractVector{<:Real}}...) where {T}
    re = MOIU.operate(+, T, (real(a) for a in args)...)
    im = MOIU.operate(+, T, (imag(a) for a in args)...)
    return ComplexTape(re, im)
end

function MOIU.operate(::typeof(conj), ::Type{T}, c::ComplexTape) where {T}
    im = MOIU.operate(-, T, imag(c))
    return ComplexTape(real(c), im)
end

function MOIU.operate(::typeof(conj), ::Type{T}, tape::VAFTapes) where {T}
    return tape
end

function MOIU.operate(::typeof(real), ::Type{T}, c::ComplexTape) where {T}
    return real(c)
end

function MOIU.operate(::typeof(real), ::Type{T}, tape::VAFTapes) where {T}
    return tape
end

function MOIU.operate(::typeof(imag), ::Type{T}, c::ComplexTape) where {T}
    return imag(c)
end

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractArray{<:Real}, c::ComplexTape) where {T}
    re = MOIU.operate(*, T, A, real(c))
    im = MOIU.operate(*, T, A, imag(c))
    return ComplexTape(re, im)
end

function MOIU.operate(::typeof(*), ::Type{T}, x::Real, c::ComplexTape) where {T}
    re = MOIU.operate(*, T, x, real(c))
    im = MOIU.operate(*, T, x, imag(c))
    return ComplexTape(re, im)
end

function MOIU.operate(::typeof(sum), ::Type{T}, c::ComplexTape) where {T}
    re = MOIU.operate(sum, T, real(c))
    im = MOIU.operate(sum, T, imag(c))
    return ComplexTape(re, im)
end


const ComplexValue = Union{<:Complex, AbstractArray{<:Complex}}

function MOIU.operate(::typeof(*), ::Type{T}, z::ComplexValue, tape::VAFTapes) where {T}
    tape = collapse(tape)
    re = MOIU.operate(*, T, real(z), tape)
    im = MOIU.operate(*, T, imag(z), tape)
    return ComplexTape(re, im)
end

function collapse(tape::ComplexTape)
    re = collapse(real(tape))
    im = collapse(imag(tape))
    return ComplexTape(re, im)
end

function MOIU.operate(::typeof(*), ::Type{T}, z::ComplexValue, tape::ComplexTape) where {T}
    tape = collapse(tape)
    re1 = MOIU.operate(*, T, real(z), real(tape))
    re2 = MOIU.operate(*, T, imag(z), imag(tape))
    re2neg = MOIU.operate(-, T, re2)
    re = MOIU.operate(+, T, re1, re2neg)
    
    im1 = MOIU.operate(*, T, imag(z), real(tape))
    im2 = MOIU.operate(*, T, real(z), imag(tape))
    im = MOIU.operate(+, T, im1, im2)
    
    return ComplexTape(re, im)
end

function MOI_add_constraint(model, f::ComplexTape, set::Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives})
    MOI_add_constraint(model, to_vaf(real(f)), set)
    MOI_add_constraint(model, to_vaf(imag(f)), set)
    return nothing
end
