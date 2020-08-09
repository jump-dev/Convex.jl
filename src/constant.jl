#############################################################################
# constant.jl
# Defines Constant, which is a subtype of AbstractExpr
#############################################################################

ispos(x::Real) = x >= 0
ispos(v::AbstractVecOrMat{<:Real}) = all(ispos, v)
isneg(x::Real) = x <= 0
isneg(v::AbstractVecOrMat{<:Real}) = all(isneg, v)

_size(x::Number) = (1, 1)
_size(x::AbstractVector) = (length(x), 1)
_size(x::AbstractMatrix) = size(x)

_sign(x::Union{Complex,AbstractVecOrMat{<:Complex}}) = ComplexSign()
function _sign(x::Value)
    if ispos(x)
        Positive()
    elseif isneg(x)
        Negative()
    else
        NoSign()
    end
end

struct Constant{T<:Value} <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    value::T
    size::Tuple{Int, Int}
    sign::Sign

    function Constant(x::Value, sign::Sign)
        x isa Complex && error("Real values expected")
        x isa AbstractArray && eltype(x) isa Complex && error("Real values expected")
        return new{typeof(x)}(:constant, objectid(x), x, _size(x), sign)
    end
    Constant(x::Value, check_sign::Bool=true) = Constant(x, check_sign ? _sign(x) : NoSign())
end

struct ComplexConstant{T<:Value} <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    real_constant::Constant{T}
    imag_constant::Constant{T}
    function ComplexConstant(re::Constant{T}, im::Constant{T}) where {T}
        size(re) == size(im) || error("size mismatch")
        new{T}(:complex_constant, rand(UInt64), size(re), re, im)
    end
end

AbstractTrees.children(c::ComplexConstant) = tuple()
vexity(::ComplexConstant) = ConstVexity()
sign(::ComplexConstant) = ComplexSign()

evaluate(c::ComplexConstant) = evaluate(c.real_constant) + im * evaluate(c.imag_constant)

struct ComplexStructOfVec{T <: AbstractVector{<:Real}}
    real_vec::T
    imag_vec::T
end

Base.real(c::ComplexStructOfVec) = c.real_vec
Base.imag(c::ComplexStructOfVec) = c.imag_vec

function template(C::ComplexConstant, context::Context)
    return ComplexStructOfVec(template(C.real_constant, context), template(C.imag_constant, context))
end

constant(x::Value) = Constant(x)
constant(x::AbstractArray{<:Complex}) = ComplexConstant(Constant(real(x)), Constant(imag(x)))
constant(x::Complex) = ComplexConstant(Constant(real(x)), Constant(imag(x)))

#### Constant Definition end     #####

vexity(::Constant) = ConstVexity()

# Lower (1,1)-matrices to scalars and (d,1)-matrices to vectors, for outputting to the user.
function output(x::Value)
    if size(x, 2) == 1
        if size(x, 1) == 1
            return x[]
        else
            return vec(x)
        end
    else
        return x
    end
end

evaluate(x::Constant) = output(x.value)

sign(x::Constant) = x.sign

# `real(::Real)` is a no-op and should be optimized out for `Constant{<:Real}`
real_conic_form(x::Constant{<:Number}) = [real(x.value)]
real_conic_form(x::Constant{<:AbstractVecOrMat}) = vec(real(x.value))

# `imag(::Real)` always returns 0, so we can avoid the implicit conversion to `Complex`
# by multiplication with `im` and just use an explicit call to `zeros` with the appropriate
# length
imag_conic_form(x::Constant{T}) where {T<:Real} = zeros(Complex{T}, 1)
imag_conic_form(x::Constant{<:AbstractVecOrMat{T}}) where {T<:Real} = zeros(Complex{T}, length(x))
imag_conic_form(x::Constant{<:Complex}) = [im * imag(x.value)]
imag_conic_form(x::Constant{<:AbstractVecOrMat{<:Complex}}) = im * vec(imag(x.value))

# We can more efficiently get the length of a constant by asking for the length of its
# value, which Julia can get via Core.arraylen for arrays and knows is 1 for scalars
length(x::Constant) = length(x.value)

function template(C::Constant, ::Context{T}) where T
    # this should happen at `Constant` creation?
    # No, we don't have access to `T` yet; that's problem-specific
    if eltype(C.value) != T
        C = Constant( T.(C.value) )
    end 
    return vectorize(C.value)
end
