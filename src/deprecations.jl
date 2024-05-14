# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

# TODO(odow): Remove code in a future release

function Base.:<(x::AbstractExpr, y::AbstractExpr)
    @warn(
        "Strict inequality operator `<` has never enforced the strict " *
        "inequality, and instead adds the non-strict inequality `<=`. This " *
        "operator will be removed in the next release.",
        maxlog = 1,
    )
    return x <= y
end

Base.:<(lhs::AbstractExpr, rhs::Value) = <(lhs, constant(rhs))

Base.:<(lhs::Value, rhs::AbstractExpr) = <(constant(lhs), rhs)

function Base.:>(x::AbstractExpr, y::AbstractExpr)
    @warn(
        "Strict inequality operator `>` has never enforced the strict " *
        "inequality, and instead adds the non-strict inequality `>=`. This " *
        "operator will be removed in the next release.",
        maxlog = 1,
    )
    return x >= y
end

Base.:>(lhs::AbstractExpr, rhs::Value) = >(lhs, constant(rhs))

Base.:>(lhs::Value, rhs::AbstractExpr) = >(constant(lhs), rhs)

@deprecate norm_inf(x::AbstractExpr) norm(x, Inf)

@deprecate norm_1(x::AbstractExpr) norm(x, 1)

@deprecate norm_fro(x::AbstractExpr) norm(x, 2)

@deprecate get_vectorized_size(x::AbstractExpr) length(x)

@deprecate operatornorm(x::AbstractExpr) LinearAlgebra.opnorm(x)

function Base.in(x::AbstractExpr, y::Symbol)
    if !(y in (:semidefinite, :SDP))
        error("Set $y not understood")
    end
    @warn(
        "Using `x in $y` to construct a semidefinite constraint is " *
        "deprecated. Use `isposdef(x)` or `x ⪰ 0` instead.",
        maxlog = 1,
    )
    return isposdef(x)
end

# Compatability with old `sets` model.
#
# Only dispatch to these methods when at least one set is given.
function Variable(
    size::Tuple{Int,Int},
    sign::Sign,
    set::Symbol,
    sets::Symbol...,
)
    @warn(
        "Using symbols in `Variable` constructor is deprecated. Use " *
        "`Convex.BinVar` or `Convex.IntVar` instead.",
        maxlog = 1,
    )
    sets = [set, sets...]
    x = if :Bin in sets
        Variable(size, sign, BinVar)
    elseif :Int in sets
        Variable(size, sign, IntVar)
    else
        Variable(size, sign, ContVar)
    end
    if :Semidefinite in sets
        add_constraint!(x, x ⪰ 0)
    end
    return x
end

function Variable(sign::Sign, set::Symbol, sets::Symbol...)
    return Variable((1, 1), sign, set, sets...)
end

function Variable(size::Tuple{Int,Int}, set::Symbol, sets::Symbol...)
    return Variable(size, NoSign(), set, sets...)
end

function Variable(m::Int, set::Symbol, sets::Symbol...)
    return Variable((m, 1), NoSign(), set, sets...)
end

function Variable(set::Symbol, sets::Symbol...)
    return Variable((1, 1), NoSign(), set, sets...)
end

function ComplexVariable(size::Tuple{Int,Int}, set::Symbol, sets::Symbol...)
    @warn(
        "Using symbols in `ComplexVariable` constructor is deprecated. Use " *
        "`isposdef(x)` to enforce semidefiniteness instead of `:Semidefinite.",
        maxlog = 1,
    )
    sets = [set, sets...]
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

function ComplexVariable(set::Symbol, sets::Symbol...)
    return ComplexVariable((1, 1), set, sets...)
end

# `+` on constraints
function warn_deprecated_constraint_concatenation()
    @warn(
        "Concatenating collections of constraints together with `+` or `+=` to produce a new list of constraints is deprecated. Instead, use `vcat` to concatenate collections of constraints.",
        maxlog = 1
    )
    return
end

function Base.:+(x::Array{<:Constraint}, y::Array{<:Constraint})
    warn_deprecated_constraint_concatenation()
    return vcat(x, y)
end

function Base.:+(x::Constraint, y::Constraint)
    @warn(
        "Adding constraints together (with `+` or `+=`) to produce a list of constraints is deprecated. Instead, construct a list of constraints via `[constraint1, constraint2]`",
        maxlog = 1
    )
    return [x, y]
end

function Base.:+(x::Constraint, y::Array{<:Constraint})
    warn_deprecated_constraint_concatenation()
    return vcat(x, y)
end

function Base.:+(x::Array{<:Constraint}, y::Constraint)
    warn_deprecated_constraint_concatenation()
    return vcat(x, y)
end
