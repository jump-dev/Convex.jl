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
    size::Tuple{Int,Int}
    sign::Sign

    function Constant(x::Value, sign::Sign)
        return new{typeof(x)}(:constant, objectid(x), x, _size(x), sign)
    end
    function Constant(x::Value, check_sign::Bool = true)
        return Constant(x, check_sign ? _sign(x) : NoSign())
    end
end

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

evaluate(x::Value) = output(x)

sign(x::Constant) = x.sign

# `real(::Real)` is a no-op and should be optimized out for `Constant{<:Real}`
real_conic_form(x::Constant{<:Number}) = [real(x.value)]
real_conic_form(x::Constant{<:AbstractVecOrMat}) = vec(real(x.value))

# `imag(::Real)` always returns 0, so we can avoid the implicit conversion to `Complex`
# by multiplication with `im` and just use an explicit call to `zeros` with the appropriate
# length
imag_conic_form(x::Constant{T}) where {T<:Real} = zeros(Complex{T}, 1)
function imag_conic_form(x::Constant{<:AbstractVecOrMat{T}}) where {T<:Real}
    return zeros(Complex{T}, length(x))
end
imag_conic_form(x::Constant{<:Complex}) = [im * imag(x.value)]
function imag_conic_form(x::Constant{<:AbstractVecOrMat{<:Complex}})
    return im * vec(imag(x.value))
end

# We can more efficiently get the length of a constant by asking for the length of its
# value, which Julia can get via Core.arraylen for arrays and knows is 1 for scalars
length(x::Constant) = length(x.value)

function conic_form!(x::Constant, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, x)
        #real_Value = real_conic_form(x)
        #imag_Value = imag_conic_form(x)
        objective = ConicObj()
        objective[objectid(:constant)] =
            (real_conic_form(x), imag_conic_form(x))
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end
