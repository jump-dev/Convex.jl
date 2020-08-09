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
    abstract type AbstractVariable <: AbstractExpr end

An `AbstractVariable` should have `head` field, an `id_hash` field
and a `size` field to conform to the `AbstractExpr` interface, and
implement methods (or use the field-access fallbacks) for

* [`_value`](@ref), [`set_value!`](@ref): get or set the numeric value of the variable.
    `_value` should return `nothing` when no numeric value is set. Note:
    [`evaluate`](@ref) is the user-facing method to access the value of `x`.
* [`vexity`](@ref), [`vexity!`](@ref): get or set the `vexity` of the variable. The
    `vexity` should be `AffineVexity()` unless the variable has been
    [`fix!`](@ref)'d, in which case it is `ConstVexity()`.
* [`sign`](@ref), [`vartype`](@ref), and [`constraints`](@ref): get the `Sign`, `VarType`,
    numeric type, and a (possibly empty) vector of constraints which are
    to be applied to any problem in which the variable is used.

Optionally, also implement [`sign!`](@ref), [`vartype!`](@ref), and [`add_constraint!`](@ref)
to allow users to modify those values or add a constraint.

"""
abstract type AbstractVariable <: AbstractExpr end

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

Returns the current value of `x` if assigned; errors otherwise.
"""
evaluate(x::AbstractVariable) = _value(x) === nothing ? error("Value of the variable is yet to be calculated") : output(_value(x))

"""
    _value(x::AbstractVariable)

Raw access to the current value of `x`; used internally by Convex.jl.
"""
_value(x::AbstractVariable) = x.value


"""
    set_value!(x::AbstractVariable, v)

Sets the current value of `x` to `v`.
"""
function set_value! end

set_value!(x::AbstractVariable, ::Nothing) = x.value = nothing

function set_value!(x::AbstractVariable, v::AbstractArray)
    size(x) == size(v) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = sign(x) == ComplexSign() ? convert(Array{ComplexF64}, v) : convert(Array{Float64}, v)
end

function set_value!(x::AbstractVariable, v::AbstractVector)
    size(x, 2) == 1 || throw(DimensionMismatch("Cannot set value of a variable of size $(size(x)) to a vector"))
    size(x, 1) == length(v) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = sign(x) == ComplexSign() ? convert(Array{ComplexF64}, v) : convert(Array{Float64}, v)
end

function set_value!(x::AbstractVariable, v::Number)
    size(x) == (1,1) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = sign(x) == ComplexSign() ? convert(ComplexF64, v) : convert(Float64, v)
end
# End of `AbstractVariable` interface


# Some helpers for the constructors for `Variable`
iscomplex(sign::Sign) = sign == ComplexSign()
iscomplex(x::AbstractExpr) = iscomplex(sign(x))
iscomplex(x::Complex) = true
iscomplex(x::Real) = false
iscomplex(x::AbstractArray{<:Complex}) = true
iscomplex(x::AbstractArray{<:Real}) = false

# A concrete implementation of the `AbstractVariable` interface
mutable struct Variable <: AbstractVariable
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
    function Variable(size::Tuple{Int, Int}, sign::Sign, vartype::VarType)
        if iscomplex(sign) && vartype != ContVar
            throw(ArgumentError("vartype must be `ContVar` for complex variables; got vartype = $vartype"))
        end
        this = new(:variable, 0, nothing, size, AffineVexity(), sign, Constraint[], vartype)
        this.id_hash = objectid(this)
        return this
    end
end


struct ComplexVariable{T1, T2} <: AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int,Int}
    real_var::T1
    imag_var::T2
    constraints::Vector{Constraint}
    function ComplexVariable(v1::AbstractVariable, v2::AbstractVariable)
        size(v1) == size(v2) || throw(ArgumentError("Real and imaginary parts must have the same size"))
        if sign(v1) == ComplexSign() || sign(v2) == ComplexSign()
            throw(ArgumentError("Real and imaginary parts must be real, not complex."))
        end
        new{typeof(v1), typeof(v2)}(:ComplexVariable, rand(UInt64), size(v1), v1, v2, Constraint[])
    end
end

vartype(::ComplexVariable) = ContVar()
sign(::ComplexVariable) = ComplexSign()
vexity(c::ComplexVariable) = vexity(c.real_var) + vexity(c.imag_var)
function vexity!(c::ComplexVariable, v::Vexity)
    vexity!(c.real_var, v)
    vexity!(c.imag_var, v)
    return nothing
end

function _value(c::ComplexVariable)
    if _value(c.real_var) === nothing || _value(c.imag_var) === nothing
        return nothing
    else
        return _value(c.real_var) + im*_value(c.imag_var)
    end
end

function set_value!(c::ComplexVariable, val::Number) 
    set_value!(c.real_var, real(val))
    set_value!(c.imag_var, imag(val))
    return nothing
end

function set_value!(c::ComplexVariable, val::AbstractVector) 
    set_value!(c.real_var, real(val))
    set_value!(c.imag_var, imag(val))
    return nothing
end


function set_value!(c::ComplexVariable, val::AbstractArray) 
    set_value!(c.real_var, real(val))
    set_value!(c.imag_var, imag(val))
    return nothing
end


Variable(size::Tuple{Int, Int}, sign::Sign) = Variable(size, sign, ContVar)
Variable(size::Tuple{Int, Int}, vartype::VarType) = Variable(size, NoSign(), vartype)
Variable(size::Tuple{Int, Int}) = Variable(size, NoSign(), ContVar)
ComplexVariable(size::Tuple{Int, Int}) = ComplexVariable(Variable(size), Variable(size))

Variable(m::Int, n::Int, sign::Sign, vartype::VarType) = Variable((m,n), sign, vartype)
Variable(m::Int, n::Int, sign::Sign) = Variable((m,n), sign, ContVar)
Variable(m::Int, n::Int, vartype::VarType) = Variable((m,n), NoSign(), vartype)
Variable(m::Int, n::Int) = Variable((m,n), NoSign(), ContVar)
ComplexVariable(m::Int, n::Int) =  ComplexVariable(Variable(m, n), Variable(m, n))

Variable(m::Int, sign::Sign, vartype::VarType) = Variable((m,1), sign, vartype)
Variable(m::Int, sign::Sign) = Variable((m,1), sign, ContVar)
Variable(m::Int, vartype::VarType) = Variable((m,1), NoSign(), vartype)
Variable(m::Int) = Variable((m,1), NoSign(), ContVar)
ComplexVariable(m::Int) =  ComplexVariable(Variable(m), Variable(m))


Variable(sign::Sign, vartype::VarType) = Variable((1,1), sign, vartype)
Variable(sign::Sign) = Variable((1,1), sign, ContVar)
Variable(vartype::VarType) = Variable((1,1), NoSign(), vartype)
Variable() = Variable((1,1), NoSign(), ContVar)
ComplexVariable() =  ComplexVariable(Variable(), Variable())


###
# Compatability with old `sets` model
# Only dispatch to these methods when at least one set is given;
# otherwise dispatch to the inner constructors.
# There are choices of sign, no sign, or `ComplexVariable` and various
# ways to enter the size for each. We will start with:
# Case 1: with sign
function Variable(size::Tuple{Int, Int}, sign::Sign, first_set::Symbol, more_sets::Symbol...)
    sets = [first_set]
    append!(sets, more_sets)
    if :Bin in sets
        vartype = BinVar
    elseif :Int in sets
        vartype = IntVar
    else
        vartype = ContVar
    end
    x = Variable(size, sign, vartype)
    if :Semidefinite in sets
        add_constraint!(x, x ⪰ 0)
    end
    return x
end
Variable(m::Int, n::Int, sign::Sign, first_set::Symbol, more_sets::Symbol...) = Variable((m,n), sign, first_set, more_sets...)
Variable(m::Int, sign::Sign, first_set::Symbol, more_sets::Symbol...) = Variable((m,1), sign, first_set, more_sets...)
Variable(sign::Sign, first_set::Symbol, more_sets::Symbol...) = Variable((1,1), sign, first_set, more_sets...)

# Case 2: no sign
Variable(size::Tuple{Int, Int}, first_set::Symbol, more_sets::Symbol...) = Variable(size, NoSign(), first_set, more_sets...)
Variable(m::Int, n::Int, first_set::Symbol, more_sets::Symbol...) = Variable((m,n), NoSign(), first_set, more_sets...)
Variable(m::Int, first_set::Symbol, more_sets::Symbol...) = Variable((m,1), NoSign(), first_set, more_sets...)
Variable(first_set::Symbol, more_sets::Symbol...) = Variable((1,1), NoSign(), first_set, more_sets...)

# Case 3: ComplexVariable
# ComplexVariable(size::Tuple{Int, Int}, first_set::Symbol, more_sets::Symbol...) = Variable(size, ComplexSign(), first_set, more_sets...)
# ComplexVariable(m::Int, n::Int, first_set::Symbol, more_sets::Symbol...) = Variable((m,n), ComplexSign(), first_set, more_sets...)
# ComplexVariable(m::Int, first_set::Symbol, more_sets::Symbol...) = Variable((m,1), ComplexSign(), first_set, more_sets...)
# ComplexVariable(first_set::Symbol, more_sets::Symbol...) = Variable((1,1), ComplexSign(), first_set, more_sets...)


Semidefinite(m::Integer) = Semidefinite(m, m)

function Semidefinite(m::Integer, n::Integer)
    if m == n
        x = Variable((m, m), NoSign(), ContVar)
        add_constraint!(x, x ⪰ 0)
        return x
    else
        error("Semidefinite matrices must be square")
    end
end

Semidefinite((m, n)::Tuple{Integer, Integer}) = Semidefinite(m, n)

function HermitianSemidefinite(m::Integer, n::Integer = m)
    if m == n
        x = ComplexVariable((m, n))
        add_constraint!(x, x ⪰ 0)
        return x
    else
        error("`HermitianSemidefinite` matrices must be square")
    end
end

HermitianSemidefinite((m, n)::Tuple{Integer, Integer}) = HermitianSemidefinite(m, n)

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

"""
    fix!(x::AbstractVariable, v = value(x))

Fixes `x` to `v`. It is subsequently treated as a constant in future
optimization problems. See also [`free!`](@ref).
"""
function fix! end

function fix!(x::AbstractVariable)
    _value(x) === nothing && error("This variable has no value yet; cannot fix value to nothing!")
    vexity!(x, ConstVexity())
    x
end

function fix!(x::AbstractVariable, v)
    set_value!(x, v)
    fix!(x)
end

"""
    free!(x::AbstractVariable)

Frees a previously [`fix!`](@ref)'d variable `x`, to treat it once again as a
variable to optimize over.
"""
function free! end

function free!(x::AbstractVariable)
    vexity!(x, AffineVexity())
    x
end
