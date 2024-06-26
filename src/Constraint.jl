# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct Constraint{S<:MOI.AbstractVectorSet}
    child::AbstractExpr
    set::S
    dual::Union{Value,Nothing}

    function Constraint(child::AbstractExpr, set::MOI.AbstractVectorSet)
        return new{typeof(set)}(child, set, nothing)
    end
end

function Constraint{S}(child::AbstractExpr) where {S<:MOI.AbstractVectorSet}
    return Constraint(child, set_with_size(S, size(child)))
end

iscomplex(c::Constraint) = iscomplex(c.child)

function set_with_size(
    ::Type{S},
    sz::Tuple{Int,Int},
) where {S<:MOI.AbstractVectorSet}
    if sz[2] != 1
        error(
            "Cannot constrain a matrix of size `$sz` to be long to the cone " *
            "`$S`, there should be only one column.",
        )
    end
    return MOI.Utilities.set_with_dimension(S, sz[1])
end

head(io::IO, c::Constraint) = head(io, c.set)

function head(io::IO, set::MOI.AbstractVectorSet)
    return print(io, replace("$(typeof(set))", "MathOptInterface" => "MOI"))
end

AbstractTrees.children(c::Constraint) = (c.child,)

# A default fallback which means that we are unsure.
is_feasible(x, set, tol) = missing

function is_feasible(
    x::Vector,
    set::Union{
        MOI.Nonnegatives,
        MOI.Nonpositives,
        MOI.Zeros,
        MOI.SecondOrderCone,
        MOI.RotatedSecondOrderCone,
        MOI.ExponentialCone,
        MOI.DualExponentialCone,
        MOI.PowerCone,
        MOI.DualPowerCone,
        MOI.GeometricMeanCone,
        MOI.NormCone,
        MOI.NormInfinityCone,
        MOI.NormOneCone,
    },
    tol,
)
    return MOI.Utilities.distance_to_set(x, set) <= tol
end

function is_feasible(x::AbstractMatrix, set::MOI.AbstractVectorSet, tol)
    return is_feasible(vec(x), set, tol)
end

function is_feasible(x::Number, set::MOI.AbstractVectorSet, tol)
    return is_feasible([x], set, tol)
end

vexity(c::Constraint) = vexity(vexity(c.child), c.set)

function vexity(
    vex,
    # An enumeration of sets that are most likely to be used by Convex.jl
    ::Union{
        MOI.SecondOrderCone,
        MOI.RotatedSecondOrderCone,
        MOI.ExponentialCone,
        MOI.DualExponentialCone,
        MOI.PowerCone,
        MOI.DualPowerCone,
        MOI.PositiveSemidefiniteConeSquare,
        MOI.GeometricMeanCone,
        MOI.NormCone,
        MOI.NormInfinityCone,
        MOI.NormOneCone,
        MOI.NormSpectralCone,
        MOI.NormNuclearCone,
        MOI.RelativeEntropyCone,
        MOI.LogDetConeSquare,
        MOI.RootDetConeSquare,
    },
)
    if !(vex == ConstVexity() || vex == AffineVexity())
        return NotDcp()
    end
    return ConvexVexity()
end

function vexity(::Any, set::MOI.AbstractVectorSet)
    return error(
        "`Convex.vexity(vex, ::$(typeof(set)))`: is not yet implemented. Please open an issue at https://github.com/jump-dev/Convex.jl",
    )
end

function _add_constraint!(context::Context, c::Constraint)
    if vexity(c.child) == ConstVexity()
        # This `evaluate` call is safe, since even if it refers to a `fix!`'d
        # variable, it happens when we are formulating the problem (not at
        # expression-time), so there is no time for the variable to be
        # re-`fix!`'d to a different value (or `free!`'d)
        feas = is_feasible(evaluate(c.child), c.set, CONSTANT_CONSTRAINT_TOL[])
        # There are three possible values of feas: true, false, and missing.
        if feas === true
            # Do nothing. The constraint is satisfied. Do not add it to the
            # solver.
            return
        elseif feas === false
            # We have proven the constraint is not satisfied. Set a flag and
            # bail. We don't need to add it to the solver.
            context.detected_infeasible_during_formulation = true
            return
        else
            # The feasibility check was unsure, likely because a method was
            # missing. Pass the constraint to the solver.
            @assert ismissing(feas)
        end
    end
    f = conic_form!(context, c.child)
    context.constr_to_moi_inds[c] = MOI_add_constraint(context.model, f, c.set)
    return
end

function populate_dual!(
    model::MOI.ModelLike,
    c::Constraint,
    indices::MOI.ConstraintIndex,
)
    ret = MOI.get(model, MOI.ConstraintDual(), indices)
    c.dual = output(reshape(ret, c.child.size))
    return
end

function populate_dual!(model::MOI.ModelLike, c::Constraint, indices::NTuple{2})
    re = MOI.get(model, MOI.ConstraintDual(), indices[1])
    imag = MOI.get(model, MOI.ConstraintDual(), indices[2])
    c.dual = output(reshape(re + im * imag, c.child.size))
    return
end

function _promote_size(lhs::AbstractExpr, rhs::AbstractExpr)
    if lhs.size == rhs.size || lhs.size == (1, 1)
        sz = rhs.size
        if lhs.size == (1, 1) && rhs.size != (1, 1)
            lhs = lhs * ones(rhs.size)
        end
    elseif rhs.size == (1, 1)
        sz = lhs.size
        if rhs.size == (1, 1) && lhs.size != (1, 1)
            rhs = rhs * ones(lhs.size)
        end
    else
        error(
            "Cannot create constraint between expressions of size " *
            "$(lhs.size) and $(rhs.size)",
        )
    end
    return lhs, rhs
end

# ==============================================================================
#     Nonnegatives
# ==============================================================================

function set_with_size(::Type{MOI.Nonnegatives}, sz::Tuple{Int,Int})
    return MOI.Nonnegatives(prod(sz))
end

head(io::IO, ::MOI.Nonnegatives) = print(io, "≥")

function vexity(vex, ::MOI.Nonnegatives)
    if vex == ConvexVexity()
        return NotDcp()
    elseif vex == ConcaveVexity()
        return ConvexVexity()
    end
    return vex
end

function Base.:>=(lhs::AbstractExpr, rhs::AbstractExpr)
    if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
        error(
            "Cannot create constraint between expressions of sign " *
            "$(sign(lhs)) and $(sign(rhs))",
        )
    end
    lhs, rhs = _promote_size(lhs, rhs)
    return Constraint{MOI.Nonnegatives}(lhs - rhs)
end

Base.:>=(lhs::AbstractExpr, rhs::Value) = >=(lhs, constant(rhs))

Base.:>=(lhs::Value, rhs::AbstractExpr) = >=(constant(lhs), rhs)

# ==============================================================================
#     Nonnpositives
# ==============================================================================

function set_with_size(::Type{MOI.Nonpositives}, sz::Tuple{Int,Int})
    return MOI.Nonpositives(prod(sz))
end

head(io::IO, ::MOI.Nonpositives) = print(io, "≤")

function vexity(vex, ::MOI.Nonpositives)
    if vex == ConcaveVexity()
        return NotDcp()
    end
    return vex
end

function Base.:<=(lhs::AbstractExpr, rhs::AbstractExpr)
    if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
        error(
            "Cannot create constraint between expressions of sign " *
            "$(sign(lhs)) and $(sign(rhs))",
        )
    end
    lhs, rhs = _promote_size(lhs, rhs)
    return Constraint{MOI.Nonpositives}(lhs - rhs)
end

Base.:<=(lhs::AbstractExpr, rhs::Value) = <=(lhs, constant(rhs))

Base.:<=(lhs::Value, rhs::AbstractExpr) = <=(constant(lhs), rhs)

# ==============================================================================
#     Zeros
# ==============================================================================

function set_with_size(::Type{MOI.Zeros}, sz::Tuple{Int,Int})
    return MOI.Zeros(prod(sz))
end

head(io::IO, ::MOI.Zeros) = print(io, "==")

function vexity(vex, ::MOI.Zeros)
    if vex == ConvexVexity() || vex == ConcaveVexity()
        return NotDcp()
    end
    return vex
end

function Base.:(==)(lhs::AbstractExpr, rhs::AbstractExpr)
    lhs, rhs = _promote_size(lhs, rhs)
    return Constraint{MOI.Zeros}(lhs - rhs)
end

Base.:(==)(lhs::AbstractExpr, rhs::Value) = ==(lhs, constant(rhs))

Base.:(==)(lhs::Value, rhs::AbstractExpr) = ==(constant(lhs), rhs)

# ==============================================================================
#     PositiveSemidefiniteConeSquare
# ==============================================================================

function set_with_size(
    ::Type{MOI.PositiveSemidefiniteConeSquare},
    sz::Tuple{Int,Int},
)
    if sz[1] != sz[2]
        error("Positive semidefinite expressions must be square")
    end
    return MOI.PositiveSemidefiniteConeSquare(sz[1])
end

head(io::IO, ::MOI.PositiveSemidefiniteConeSquare) = print(io, "PSD")

function is_feasible(
    x::AbstractMatrix,
    ::MOI.PositiveSemidefiniteConeSquare,
    tol,
)
    return x ≈ transpose(x) && LinearAlgebra.eigmin(x) >= -tol
end

function LinearAlgebra.isposdef(x::AbstractExpr)
    if iscomplex(x)
        return Constraint{MOI.PositiveSemidefiniteConeSquare}(
            [real(x) -imag(x); imag(x) real(x)],
        )
    end
    return Constraint{MOI.PositiveSemidefiniteConeSquare}(x)
end

⪰(x::AbstractExpr, y::AbstractExpr) = isposdef(x - y)

function ⪰(x::AbstractExpr, y::Value)
    if all(y .== 0)
        return isposdef(x)
    end
    return isposdef(x - constant(y))
end

function ⪰(x::Value, y::AbstractExpr)
    if all(x .== 0)
        return isposdef(-y)
    end
    return isposdef(constant(x) - y)
end

⪯(x::AbstractExpr, y::AbstractExpr) = ⪰(y, x)

⪯(x::Value, y::AbstractExpr) = ⪰(y, x)

⪯(x::AbstractExpr, y::Value) = ⪰(y, x)
