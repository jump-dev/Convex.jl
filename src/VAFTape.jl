
function operate(op, ::Type{T}, sign::Sign, args...) where {T}
    if iscomplex(sign)
        return complex_operate(op, T, args...)
    else
        return real_operate(op, T, args...)
    end
end

# fallbacks
real_operate(op::F, ::Type{T}, args...) where {F,T} = op(args...)
complex_operate(op::F, ::Type{T}, args...) where {F,T} = op(args...)

function real_operate(::typeof(vcat), ::Type{T}, args...) where {T}
    @warn "real vcat fallback hit" typeof.(args) args
    return vcat(args...)
end

function complex_operate(::typeof(vcat), ::Type{T}, args...) where {T}
    @warn "complex vcat fallback hit" typeof.(args) args
    return vcat(args...)
end

# operate(op::F, ::Type{T}, args...) where {F,T} = (@show(args); op(args...))

# function real_operate(
#     ::typeof(+),
#     ::Type{T},
#     args::AbstractVector...,
# ) where {T<:Real}
#     return sum(args)
# end
