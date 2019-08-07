#############################################################################
# variable.jl
# Defines Variable, which is a subtype of AbstractExpr
#############################################################################

"""
    VarType

Describe the type of a [`AbstractVariable`](@ref): either continuous (`ContVar`),
integer-valued (`IntVar`), or binary (`BinVar`).
"""
@enum VarType BinVar IntVar ContVar

@doc "Indicates a [`AbstractVariable`](@ref) is a binary variable." BinVar
@doc "Indicates a [`AbstractVariable`](@ref) is integer-valued." IntVar
@doc "Indicates a [`AbstractVariable`](@ref) is continuous." ContVar

"""
    abstract type AbstractVariable{T <: Number} <: AbstractExpr end

An `AbstractVariable` should have `head` field, an `id_hash` field
and a `size` field to conform to the `AbstractExpr` interface, and
implement methods (or use the field-access fallbacks) for

* [`evaluate`](@ref), [`set_value!`](@ref): get or set the numeric value of the variable.
    `evaluate` should return `nothing` when no numeric value is set.
* [`vexity`](@ref), [`vexity!`](@ref): get or set the `vexity` of the variable. The
    `vexity` should be `AffineVexity()` unless the variable has been
    [`fix!`](@ref)'d, in which case it is `ConstVexity()`.
* [`sign`](@ref), [`vartype`](@ref), and [`constraints`](@ref): get the `Sign`, `VarType`,
    numeric type, and a (possibly empty) vector of constraints which are
    to be applied to any problem in which the variable is used.

Optionally, also implement [`sign!`](@ref), [`vartype!`](@ref), and [`add_constraint!`](@ref)
to allow users to modify those values or add a constraint.

The parameter `T` indicates the numeric type (e.g. `Float64`). If
the variable is a complex variable, `T` should be complex
(e.g. `Complex{Float64}`).

"""
abstract type AbstractVariable{T <: Number} <: AbstractExpr end
Base.eltype(::Type{<:AbstractVariable{T}}) where {T} = T

# Default implementation of `AbstractVariable` interface
"""
    constraints(x::AbstractVariable)

Returns the current constraints carried by `x`.
"""
constraints(x::AbstractVariable) = x.constraints

"""
    add_constraint!(x::AbstractVariable, C::Constraint)

Adds an constraint to those carried by `x`.
"""
add_constraint!(x::AbstractVariable, C::Constraint) = push!(x.constraints, C)

"""
    vartype(x::AbstractVariable)

Returns the current [`VarType`](@ref) of `x`.
"""
vartype(x::AbstractVariable) = x.vartype

"""
    vartype!(x::AbstractVariable, vt::VarType)

Sets the current [`VarType`](@ref) of `x` to `vt`.
"""
vartype!(x::AbstractVariable, vt::VarType) = x.vartype = vt

"""
    vexity(x::AbstractVariable)

Returns the current vexity of `x`.
"""
vexity(x::AbstractVariable) = x.vexity

"""
    vexity!(x::AbstractVariable, v::Vexity)

Sets the current vexity of `x` to `v`. Should only be called
by [`fix!`](@ref) and [`free!`](@ref).
"""
vexity!(x::AbstractVariable, v::Vexity) = x.vexity = v


"""
    sign(x::AbstractVariable)

Returns the current sign of `x`.
"""
sign(x::AbstractVariable) = x.sign

"""
    sign!(x::AbstractVariable, s::Sign)

Sets the current sign of `x` to `s`.
"""
sign!(x::AbstractVariable, s::Sign) = x.sign = s


"""
    evaluate(x::AbstractVariable)

Returns the current value of `x`.
"""
evaluate(x::AbstractVariable) = x.value === nothing ? nothing : output(x.value)

"""
    set_value!(x::AbstractVariable, v)

Sets the current value of `x` to `v`, appropriately converting `v` if necessary.
"""
function set_value! end

set_value!(x::AbstractVariable, ::Nothing) = x.value = nothing

function set_value!(x::AbstractVariable, v::AbstractArray)
    size(x) == size(v) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = convert(Array{eltype(x)}, v)
end

function set_value!(x::AbstractVariable, v::AbstractVector)
    size(x, 2) == 1 || throw(DimensionMismatch("Cannot set value of a variable of size $(size(x)) to a vector"))
    size(x, 1) == length(v) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = convert(Array{eltype(x)}, v)
end

function set_value!(x::AbstractVariable, v::Number)
    size(x) == (1,1) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = convert(eltype(x), v)
end
# End of `AbstractVariable` interface


# Some helpers for the constructors for `Variable`
iscomplex(::Type{T}) where {T <: Number} = T != real(T)
iscomplex(sign::Sign) = sign == ComplexSign()

defaultsign(::Type{T}) where {T <: Number} = iscomplex(T) ? ComplexSign() : NoSign()
defaultsign() = NoSign()
defaulttype(sign::Sign) = iscomplex(sign) ? Complex{Float64} : Float64
defaulttype() = Float64

# A concrete implementation of the `AbstractVariable` interface
mutable struct Variable{T <: Number} <: AbstractVariable{T}
    """
    Every `AbstractExpr` has a `head`; for a Variable it is set to `:variable`.
    """
    head::Symbol
    """
    A unique identifying hash used for caching.
    """
    id_hash::UInt64
    """
    The current value of the variable. Defaults to `nothing` until the
    variable has been [`fix!`](@ref)'d to a particular value, or the
    variable has been used in a problem which has been solved, at which
    point the optimal value is populated into this field.
    """
    value::ValueOrNothing
    """
    The size of the variable. Scalar variables have size `(1,1)`;
    `d`-dimensional vectors have size `(d, 1)`, and `n` by `m` matrices
    have size `(n,m)`.
    """
    size::Tuple{Int, Int}
    """
    `AffineVexity()` unless the variable is `fix!`'d, in which case it is
    `ConstVexity()`. Accessed by `vexity(v::Variable)`. To check if a
    `Variable` is fixed, use `vexity(v) == ConstVexity()`.
    """
    vexity::Vexity
    """
    The sign of the variable. Can be  `Positive()`, `Negative()`, `NoSign()`
    (i.e. real), or `ComplexSign()`. Accessed by `sign(v::Variable)`. 
    """
    sign::Sign
    """
    Vector of constraints which are enforced whenever the variable is used
    in a problem.
    """
    constraints::Vector{Constraint}
    """
    The `VarType` of the variable (binary, integer, or continuous).
    """
    vartype::VarType
    function Variable{T}(size::Tuple{Int, Int}, sign::Sign, vartype::VarType) where {T <: Number}
        if iscomplex(T) != iscomplex(sign)
            throw(ArgumentError("type parameter and sign must agree (both complex or both real). Got T=$T, sign = $sign"))
        end
        if iscomplex(T) && vartype != ContVar
            throw(ArgumentError("vartype must be `ContVar` for complex variables; got vartype = $vartype"))
        end
        this = new{T}(:variable, 0, nothing, size, AffineVexity(), sign, Constraint[], vartype)

        this.id_hash = objectid(this)
        return this
    end
end
Variable{T}(size::Tuple{Int, Int}, sign::Sign) where {T <: Number} = Variable{T}(size, sign, ContVar)
Variable{T}(size::Tuple{Int, Int}, vartype::VarType) where {T <: Number} = Variable{T}(size, defaultsign(T), vartype)
Variable{T}(size::Tuple{Int, Int}) where {T <: Number} = Variable{T}(size, defaultsign(T), ContVar)

Variable(size::Tuple{Int, Int}, sign::Sign, vartype::VarType) = Variable{defaulttype(sign)}(size, sign, vartype)
Variable(size::Tuple{Int, Int}, sign::Sign) = Variable{defaulttype(sign)}(size, sign, ContVar)
Variable(size::Tuple{Int, Int}, vartype::VarType) = Variable{defaulttype()}(size, defaultsign(), vartype)
Variable(size::Tuple{Int, Int}) = Variable{defaulttype()}(size, defaultsign(), ContVar)

Variable{T}(m::Int, n::Int, sign::Sign, vartype::VarType) where {T <: Number} = Variable{T}((m,n), sign, vartype)
Variable{T}(m::Int, n::Int, sign::Sign) where {T <: Number} = Variable{T}((m,n), sign, ContVar)
Variable{T}(m::Int, n::Int, vartype::VarType) where {T <: Number} = Variable{T}((m,n), defaultsign(T), vartype)
Variable{T}(m::Int, n::Int) where {T <: Number} = Variable{T}((m,n), defaultsign(T), ContVar)

Variable(m::Int, n::Int, sign::Sign, vartype::VarType) = Variable{defaulttype(sign)}((m,n), sign, vartype)
Variable(m::Int, n::Int, sign::Sign) = Variable{defaulttype(sign)}((m,n), sign, ContVar)
Variable(m::Int, n::Int, vartype::VarType) = Variable{defaulttype()}((m,n), defaultsign(), vartype)
Variable(m::Int, n::Int) = Variable{defaulttype()}((m,n), defaultsign(), ContVar)

Variable{T}(m::Int, sign::Sign, vartype::VarType) where {T <: Number} = Variable{T}((m,1), sign, vartype)
Variable{T}(m::Int, sign::Sign) where {T <: Number} = Variable{T}((m,1), sign, ContVar)
Variable{T}(m::Int, vartype::VarType) where {T <: Number} = Variable{T}((m,1), defaultsign(T), vartype)
Variable{T}(m::Int) where {T <: Number} = Variable{T}((m,1), defaultsign(T), ContVar)

Variable(m::Int, sign::Sign, vartype::VarType) = Variable{defaulttype(sign)}((m,1), sign, vartype)
Variable(m::Int, sign::Sign) = Variable{defaulttype(sign)}((m,1), sign, ContVar)
Variable(m::Int, vartype::VarType) = Variable{defaulttype()}((m,1), defaultsign(), vartype)
Variable(m::Int) = Variable{defaulttype()}((m,1), defaultsign(), ContVar)

Variable{T}(sign::Sign, vartype::VarType) where {T <: Number} = Variable{T}((1,1), sign, vartype)
Variable{T}(sign::Sign) where {T <: Number} = Variable{T}((1,1), sign, ContVar)
Variable{T}(vartype::VarType) where {T <: Number} = Variable{T}((1,1), defaultsign(T), vartype)
Variable{T}() where {T <: Number} = Variable{T}((1,1), defaultsign(T), ContVar)

Variable(sign::Sign, vartype::VarType) = Variable{defaulttype(sign)}((1,1), sign, vartype)
Variable(sign::Sign) = Variable{defaulttype(sign)}((1,1), sign, ContVar)
Variable(vartype::VarType) = Variable{defaulttype()}((1,1), defaultsign(), vartype)
Variable() = Variable{defaulttype()}((1,1), defaultsign(), ContVar)


###
# Compatability with old `sets` model
# Only dispatch to these methods when at least one set is given;
# otherwise dispatch to the inner constructors.
# There are choices of sign or no sign, type parameter or no type parameter, and various
# ways to enter the size for each. We will start with:
# Case 1: with sign, with type parameter.
function Variable{T}(size::Tuple{Int, Int}, sign::Sign, first_set::Symbol, more_sets::Symbol...)  where {T <: Number}
    sets = [first_set]
    append!(sets, more_sets)
    if :Bin in sets
        vartype = BinVar
    elseif :Int in sets
        vartype = IntVar
    else
        vartype = ContVar
    end
    x = Variable{T}(size, sign, vartype)
    if :Semidefinite in sets
        add_constraint!(x, x ⪰ 0)
    end
    return x
end
Variable{T}(m::Int, n::Int, sign::Sign, first_set::Symbol, more_sets::Symbol...)  where {T <: Number} = Variable{T}((m,n), sign, first_set, more_sets...)
Variable{T}(m::Int, sign::Sign, first_set::Symbol, more_sets::Symbol...)  where {T <: Number} = Variable{T}((m,1), sign, first_set, more_sets...)
Variable{T}(sign::Sign, first_set::Symbol, more_sets::Symbol...)  where {T <: Number} = Variable{T}((1,1), sign, first_set, more_sets...)

# Case 2: with sign, no type parameter.
Variable(size::Tuple{Int, Int}, sign::Sign, first_set::Symbol, more_sets::Symbol...) = Variable{defaulttype(sign)}(size, sign, first_set, more_sets...)
Variable(m::Int, n::Int, sign::Sign, first_set::Symbol, more_sets::Symbol...) = Variable{defaulttype(sign)}((m,n), sign, first_set, more_sets...)
Variable(m::Int, sign::Sign, first_set::Symbol, more_sets::Symbol...) = Variable{defaulttype(sign)}((m,1), sign, first_set, more_sets...)
Variable(sign::Sign, first_set::Symbol, more_sets::Symbol...) = Variable{defaulttype(sign)}((1,1), sign, first_set, more_sets...)

# Case 3: no sign, type parameter
Variable{T}(size::Tuple{Int, Int}, first_set::Symbol, more_sets::Symbol...) where {T <: Number} = Variable{T}(size, defaultsign(T), first_set, more_sets...)
Variable{T}(m::Int, n::Int, first_set::Symbol, more_sets::Symbol...) where {T <: Number} = Variable{T}((m,n), defaultsign(T), first_set, more_sets...)
Variable{T}(m::Int, first_set::Symbol, more_sets::Symbol...) where {T <: Number} = Variable{T}((m,1), defaultsign(T), first_set, more_sets...)
Variable{T}(first_set::Symbol, more_sets::Symbol...) where {T <: Number} = Variable{T}((1,1), defaultsign(T), first_set, more_sets...)

# Case 4: no sign, no type parameter
Variable(size::Tuple{Int, Int}, first_set::Symbol, more_sets::Symbol...) = Variable{defaulttype()}(size, defaultsign(), first_set, more_sets...)
Variable(m::Int, n::Int, first_set::Symbol, more_sets::Symbol...) = Variable{defaulttype()}((m,n), defaultsign(), first_set, more_sets...)
Variable(m::Int, first_set::Symbol, more_sets::Symbol...) = Variable{defaulttype()}((m,1), defaultsign(), first_set, more_sets...)
Variable(first_set::Symbol, more_sets::Symbol...) = Variable{defaulttype()}((1,1), defaultsign(), first_set, more_sets...)


Semidefinite(m::Integer, vartype::VarType = ContVar; eltype::Type{T} = Float64) where {T <: Number} = Semidefinite(m, m, vartype; eltype=T)

function Semidefinite(m::Integer, n::Integer, vartype::VarType = ContVar; eltype::Type{T} = Float64) where {T <: Number}
    if m == n
        x = Variable{T}((m, m), NoSign(), vartype)
        add_constraint!(x, x ⪰ 0)
        return x
    else
        error("Semidefinite matrices must be square")
    end
end

# The `defaultsign(T)` calls in the constructors will take care of setting the `sign`
ComplexVariable(args...; eltype::Type{T}=ComplexF64) where {T <: Complex{<:Real}} = Variable{T}(args...)
 
function HermitianSemidefinite(m::Integer, vartype::VarType = ContVar; eltype::Type{T}=ComplexF64) where {T <: Complex{<:Real}}
    return HermitianSemidefinite(m,m, vartype; eltype = eltype)
end

function HermitianSemidefinite(m::Integer, n::Integer, vartype::VarType = ContVar; eltype::Type{T}=ComplexF64) where {T <: Complex{<:Real}}
    if m == n
        x = Variable{eltype}((m, n), ComplexSign(), vartype)
        add_constraint!(x, x ⪰ 0)
        return x
    else
        error("`HermitianSemidefinite` matrices must be square")
    end
end


function real_conic_form(x::AbstractVariable)
    vec_size = length(x)
    return sparse(1.0I, vec_size, vec_size)
end

function imag_conic_form(x::AbstractVariable)
    vec_size = length(x)
    if sign(x) == ComplexSign()
        return im*sparse(1.0I, vec_size, vec_size)
    else
        return spzeros(vec_size, vec_size)
    end
end

# fix variables to hold them at their current value, and free them afterwards
function fix!(x::AbstractVariable)
    evaluate(x) === nothing && error("This variable has no value yet; cannot fix value to nothing!")
    vexity!(x, ConstVexity())
    x
end

function fix!(x::AbstractVariable, v)
    set_value!(x, v)
    fix!(x)
end

function free!(x::AbstractVariable)
    vexity!(x, AffineVexity())
    x
end
