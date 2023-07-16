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

mutable struct Constant{T<:Real} <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    value::Matrix{T}
    size::Tuple{Int,Int}
    sign::Sign

    function Constant(x::Value, sign::Sign)
        x isa Complex && error("Real values expected")
        x isa AbstractArray &&
            eltype(x) <: Complex &&
            error("Real values expected")

        # Convert to matrix
        mat = [x;;]
        return new{eltype(x)}(:constant, objectid(x), mat, _size(x), sign)
    end
    function Constant(x::Value, check_sign::Bool = true)
        return Constant(x, check_sign ? _sign(x) : NoSign())
    end
end
# Constant(x::Constant) = x

mutable struct ComplexConstant{T<:Real} <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int,Int}
    real_constant::Constant{T}
    imag_constant::Constant{T}
    function ComplexConstant(re::Constant{T}, im::Constant{T}) where {T}
        size(re) == size(im) || error("size mismatch")
        return new{T}(:complex_constant, rand(UInt64), size(re), re, im)
    end

    # function ComplexConstant(re::Constant{S1}, im::Constant{S2}) where {S1,S2}
    #     size(re) == size(im) || error("size mismatch")
    #     re, im = promote(re.value, im.value)
    #     re = Constant(re)
    #     im = Constant(im)
    #     return new{T}(:complex_constant, rand(UInt64), size(re), re, im)
    # end
end

AbstractTrees.children(c::ComplexConstant) = tuple()
vexity(::ComplexConstant) = ConstVexity()
sign(::ComplexConstant) = ComplexSign()

function evaluate(c::ComplexConstant)
    return evaluate(c.real_constant) + im * evaluate(c.imag_constant)
end

mutable struct ComplexStructOfVec{T<:Real}
    real_vec::Vector{T}
    imag_vec::Vector{T}
end

Base.real(c::ComplexStructOfVec) = c.real_vec
Base.imag(c::ComplexStructOfVec) = c.imag_vec
Base.conj(c::ComplexStructOfVec) = ComplexStructOfVec(real(c), -imag(c))

function _conic_form!(context::Context, C::ComplexConstant)
    return ComplexStructOfVec(
        conic_form!(context, C.real_constant),
        conic_form!(context, C.imag_constant),
    )
end

function constant(x)
    # Convert to matrix
    x = [x;;]
    if eltype(x) <: Real
        return Constant(x)
    else
        return ComplexConstant(Constant(real(x)), Constant(imag(x)))
    end
end
# constant(x::Complex) = ComplexConstant(Constant(real(x)), Constant(imag(x)))

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

# We can more efficiently get the length of a constant by asking for the length of its
# value, which Julia can get via Core.arraylen for arrays and knows is 1 for scalars
length(x::Constant) = length(x.value)

function _conic_form!(::Context{T}, C::Constant) where {T}
    # this should happen at `Constant` creation?
    # No, we don't have access to `T` yet; that's problem-specific
    if eltype(C.value) != T
        C = Constant(convert(Matrix{T}, C.value))
    end
    return vec(C.value)
end
