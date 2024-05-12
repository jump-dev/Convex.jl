# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

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

_matrix(x::AbstractArray) = Matrix(x)
_matrix(x::AbstractVector) = reshape(Vector(x), length(x), 1)
_matrix(x::Number) = _matrix([x])
_matrix(x::SparseArrays.AbstractSparseMatrix) = SparseArrays.sparse(x)
function _matrix(x::SparseArrays.AbstractSparseVector)
    return SparseArrays.sparse(reshape(x, length(x), 1))
end

mutable struct Constant{T<:Real} <: AbstractExpr
    head::Symbol
    value::Union{Matrix{T},SPARSE_MATRIX{T}}
    size::Tuple{Int,Int}
    sign::Sign

    function Constant(x::Value, sign::Sign)
        if x isa Complex || x isa AbstractArray{<:Complex}
            throw(DomainError(x, "Constant expects real values"))
        end
        return new{eltype(x)}(:constant, _matrix(x), _size(x), sign)
    end
end
function Constant(x::Value, check_sign::Bool = true)
    return Constant(x, check_sign ? _sign(x) : NoSign())
end
# Constant(x::Constant) = x

mutable struct ComplexConstant{T<:Real} <: AbstractExpr
    head::Symbol
    size::Tuple{Int,Int}
    real_constant::Constant{T}
    imag_constant::Constant{T}
    function ComplexConstant(re::Constant{T}, im::Constant{T}) where {T}
        size(re) == size(im) || error("size mismatch")
        return new{T}(:complex_constant, size(re), re, im)
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
Base.sign(::ComplexConstant) = ComplexSign()

function evaluate(c::ComplexConstant)
    return evaluate(c.real_constant) + im * evaluate(c.imag_constant)
end

mutable struct ComplexStructOfVec{T<:Real}
    real_vec::SPARSE_VECTOR{T}
    imag_vec::SPARSE_VECTOR{T}
end

Base.real(c::ComplexStructOfVec) = c.real_vec
Base.imag(c::ComplexStructOfVec) = c.imag_vec
Base.conj(c::ComplexStructOfVec) = ComplexStructOfVec(real(c), -imag(c))

function new_conic_form!(context::Context, C::ComplexConstant)
    return ComplexStructOfVec(
        conic_form!(context, C.real_constant),
        conic_form!(context, C.imag_constant),
    )
end

constant(x::Constant) = x
constant(x::ComplexConstant) = x
function constant(x)
    # Convert to matrix
    x = _matrix(x)
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

Base.sign(x::Constant) = x.sign

# We can more efficiently get the length of a constant by asking for the length of its
# value, which Julia can get via Core.arraylen for arrays and knows is 1 for scalars
Base.length(x::Constant) = length(x.value)

function new_conic_form!(::Context{T}, C::Constant) where {T}
    # this should happen at `Constant` creation?
    # No, we don't have access to `T` yet; that's problem-specific
    x = SPARSE_VECTOR{T}(vec(C.value))

    return x
end

# Only handle the real case for now
# basically a duplicate of `Constant`
# TODO- neither the sign nor the size must never change once we `formulate`
# even if the values do
mutable struct Parameter{T<:Real} <: AbstractExpr
    head::Symbol
    value::Union{Matrix{T},SPARSE_MATRIX{T}}
    size::Tuple{Int,Int}
    sign::Sign

    function Parameter(x::Value, sign::Sign)
        if x isa Complex || x isa AbstractArray{<:Complex}
            throw(DomainError(x, "Parameter expects real values"))
        end
        return new{eltype(x)}(:parameter, _matrix(x), _size(x), sign)
    end
end

function update_parameters!(context::Context{T}) where {T}
    for (p, inds) in pairs(context.constr_to_moi_inds)
        p isa Parameter || continue
        @assert length(p.value) == length(inds)
        for (elt, ci) in zip(p.value, inds)
            @info "Setting $ci to $elt"
            MOI.set(
                context.model,
                MOI.ConstraintSet(),
                ci,
                MOI.Parameter{T}(elt),
            )
        end
    end
end

vexity(::Parameter) = ConstVexity()

function Parameter(x::Value, check_sign::Bool = true)
    return Parameter(x, check_sign ? _sign(x) : NoSign())
end

function evaluate(x::Parameter)
    return output(x.value)
end

Base.sign(x::Parameter) = x.sign

Base.length(x::Parameter) = length(x.value)
AbstractTrees.children(::Parameter) = tuple()

function new_conic_form!(context::Context{T}, c::Parameter) where {T}
    @info "hi2"
    model = context.model
    variable_indices = MOI.VariableIndex[]
    constraint_indices = MOI.ConstraintIndex{
        MathOptInterface.VariableIndex,
        MathOptInterface.Parameter{T},
    }()
    for elt in c.value
        # TODO- convert?
        v, ci = MOI.add_constrained_variable(model, MOI.Parameter{T}(elt))
        push!(v, variable_indices)
        push!(constraint_indices, ci)
    end

    context.var_to_moi_indices[c] = variable_indices
    context.constr_to_moi_inds[c] = constraint_indices

    return _to_tape(MOI.VectorOfVariables(variable_indices))
end
