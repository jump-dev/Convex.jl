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
* [`sign`](@ref), [`vartype`](@ref), and [`get_constraints`](@ref): get the `Sign`, `VarType`,
    numeric type, and a (possibly empty) vector of constraints which are
    to be applied to any problem in which the variable is used.

Optionally, also implement [`sign!`](@ref), [`vartype!`](@ref), and [`add_constraint!`](@ref)
to allow users to modify those values or add a constraint.

"""
abstract type AbstractVariable <: AbstractExpr end

# Default implementation of `AbstractVariable` interface
"""
    get_constraints(x::AbstractVariable)

Returns the current constraints carried by `x`.
"""
get_constraints(x::AbstractVariable) = x.constraints

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
    Base.sign(x::AbstractVariable)

Returns the current sign of `x`.
"""
Base.sign(x::AbstractVariable) = x.sign

"""
    sign!(x::AbstractVariable, s::Sign)

Sets the current sign of `x` to `s`.
"""
sign!(x::AbstractVariable, s::Sign) = x.sign = s

"""
    evaluate(x::AbstractVariable)

Returns the current value of `x` if assigned; errors otherwise.
"""
function evaluate(x::AbstractVariable)
    if _value(x) === nothing
        error("Value of the variable is yet to be calculated")
    end
    return output(copy(_value(x)))
end

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
    size(x) == size(v) ||
        throw(DimensionMismatch("Variable and value sizes do not match!"))
    if iscomplex(x) && !(eltype(v) <: Complex)
        v = complex.(v)
    end
    return x.value = v
end

function set_value!(x::AbstractVariable, v::AbstractVector)
    size(x, 2) == 1 || throw(
        DimensionMismatch(
            "Cannot set value of a variable of size $(size(x)) to a vector",
        ),
    )
    size(x, 1) == length(v) ||
        throw(DimensionMismatch("Variable and value sizes do not match!"))
    if iscomplex(x) && !(eltype(v) <: Complex)
        v = complex.(v)
    end
    return x.value = v
end

function set_value!(x::AbstractVariable, v::Number)
    size(x) == (1, 1) ||
        throw(DimensionMismatch("Variable and value sizes do not match!"))
    if sign(x) == ComplexSign() && !(v isa Complex)
        v = complex(v)
    end
    return x.value = v
end
# End of `AbstractVariable` interface

# Some helpers for the constructors for `Variable`
iscomplex(sign::Sign) = sign == ComplexSign()
iscomplex(x::AbstractVariable) = iscomplex(sign(x))
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
    value::Union{Value,Nothing}
    """
    The size of the variable. Scalar variables have size `(1,1)`;
    `d`-dimensional vectors have size `(d, 1)`, and `n` by `m` matrices
    have size `(n,m)`.
    """
    size::Tuple{Int,Int}
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
    function Variable(size::Tuple{Int,Int}, sign::Sign, vartype::VarType)
        if iscomplex(sign)
            if vartype == ContVar
                return ComplexVariable(size)
            else
                throw(
                    ArgumentError(
                        "vartype must be `ContVar` for complex variables; got vartype = $vartype",
                    ),
                )
            end
        end
        this = new(
            :variable,
            0,
            nothing,
            size,
            AffineVexity(),
            sign,
            Constraint[],
            vartype,
        )
        this.id_hash = objectid(this)
        return this
    end
end

function Base.isequal(x::AbstractVariable, y::AbstractVariable)
    return x.id_hash == y.id_hash
end
Base.hash(x::AbstractVariable, h::UInt) = hash(x.id_hash, h)

mutable struct ComplexVariable <: AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int,Int}
    value::Union{Value,Nothing}
    vexity::Vexity
    constraints::Vector{Constraint}
    function ComplexVariable(sz::Tuple{Int,Int})
        return new(
            :ComplexVariable,
            rand(UInt64),
            sz,
            nothing,
            AffineVexity(),
            Constraint[],
        )
    end
end

vartype(::ComplexVariable) = ContVar
Base.sign(::ComplexVariable) = ComplexSign()

Variable(size::Tuple{Int,Int}, sign::Sign) = Variable(size, sign, ContVar)
function Variable(size::Tuple{Int,Int}, vartype::VarType)
    return Variable(size, NoSign(), vartype)
end
Variable(size::Tuple{Int,Int}) = Variable(size, NoSign(), ContVar)

function Variable(m::Int, n::Int, sign::Sign, vartype::VarType)
    return Variable((m, n), sign, vartype)
end
Variable(m::Int, n::Int, sign::Sign) = Variable((m, n), sign, ContVar)
Variable(m::Int, n::Int, vartype::VarType) = Variable((m, n), NoSign(), vartype)
Variable(m::Int, n::Int) = Variable((m, n), NoSign(), ContVar)
ComplexVariable(m::Int, n::Int) = ComplexVariable((m, n))

Variable(m::Int, sign::Sign, vartype::VarType) = Variable((m, 1), sign, vartype)
Variable(m::Int, sign::Sign) = Variable((m, 1), sign, ContVar)
Variable(m::Int, vartype::VarType) = Variable((m, 1), NoSign(), vartype)
Variable(m::Int) = Variable((m, 1), NoSign(), ContVar)
ComplexVariable(m::Int) = ComplexVariable((m, 1))

Variable(sign::Sign, vartype::VarType) = Variable((1, 1), sign, vartype)
Variable(sign::Sign) = Variable((1, 1), sign, ContVar)
Variable(vartype::VarType) = Variable((1, 1), NoSign(), vartype)
Variable() = Variable((1, 1), NoSign(), ContVar)
ComplexVariable() = ComplexVariable((1, 1))

###
# Compatability with old `sets` model
# Only dispatch to these methods when at least one set is given;
# otherwise dispatch to the inner constructors.
# There are choices of sign, no sign, or `ComplexVariable` and various
# ways to enter the size for each. We will start with:
# Case 1: with sign
function Variable(
    size::Tuple{Int,Int},
    sign::Sign,
    first_set::Symbol,
    more_sets::Symbol...,
)
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
function Variable(
    m::Int,
    n::Int,
    sign::Sign,
    first_set::Symbol,
    more_sets::Symbol...,
)
    return Variable((m, n), sign, first_set, more_sets...)
end
function Variable(m::Int, sign::Sign, first_set::Symbol, more_sets::Symbol...)
    return Variable((m, 1), sign, first_set, more_sets...)
end
function Variable(sign::Sign, first_set::Symbol, more_sets::Symbol...)
    return Variable((1, 1), sign, first_set, more_sets...)
end

# Case 2: no sign
function Variable(size::Tuple{Int,Int}, first_set::Symbol, more_sets::Symbol...)
    return Variable(size, NoSign(), first_set, more_sets...)
end
function Variable(m::Int, n::Int, first_set::Symbol, more_sets::Symbol...)
    return Variable((m, n), NoSign(), first_set, more_sets...)
end
function Variable(m::Int, first_set::Symbol, more_sets::Symbol...)
    return Variable((m, 1), NoSign(), first_set, more_sets...)
end
function Variable(first_set::Symbol, more_sets::Symbol...)
    return Variable((1, 1), NoSign(), first_set, more_sets...)
end

# Case 3: ComplexVariable
function ComplexVariable(
    size::Tuple{Int,Int},
    first_set::Symbol,
    more_sets::Symbol...,
)
    sets = [first_set]
    append!(sets, more_sets)
    if :Bin in sets
        throw(ArgumentError("Complex variables cannot be restricted to binary"))
    elseif :Int in sets
        throw(
            ArgumentError("Complex variables cannot be restricted to integer"),
        )
    end
    x = ComplexVariable(size)
    if :Semidefinite in sets
        add_constraint!(x, x ⪰ 0)
    end
    return x
end
function ComplexVariable(
    m::Int,
    n::Int,
    first_set::Symbol,
    more_sets::Symbol...,
)
    return ComplexVariable((m, n), first_set, more_sets...)
end
function ComplexVariable(m::Int, first_set::Symbol, more_sets::Symbol...)
    return ComplexVariable((m, 1), first_set, more_sets...)
end
function ComplexVariable(first_set::Symbol, more_sets::Symbol...)
    return ComplexVariable((1, 1), first_set, more_sets...)
end

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

Semidefinite((m, n)::Tuple{Integer,Integer}) = Semidefinite(m, n)

function HermitianSemidefinite(m::Integer, n::Integer = m)
    if m == n
        x = ComplexVariable((m, n))
        add_constraint!(x, x ⪰ 0)
        return x
    else
        error("`HermitianSemidefinite` matrices must be square")
    end
end

function HermitianSemidefinite((m, n)::Tuple{Integer,Integer})
    return HermitianSemidefinite(m, n)
end

"""
    fix!(x::AbstractVariable, v = value(x))

Fixes `x` to `v`. It is subsequently treated as a constant in future
optimization problems. See also [`free!`](@ref).
"""
function fix! end

function fix!(x::AbstractVariable)
    _value(x) === nothing &&
        error("This variable has no value yet; cannot fix value to nothing!")
    vexity!(x, ConstVexity())
    return x
end

function fix!(x::AbstractVariable, v)
    set_value!(x, v)
    return fix!(x)
end

"""
    free!(x::AbstractVariable)

Frees a previously [`fix!`](@ref)'d variable `x`, to treat it once again as a
variable to optimize over.
"""
function free! end

function free!(x::AbstractVariable)
    vexity!(x, AffineVexity())
    return x
end
