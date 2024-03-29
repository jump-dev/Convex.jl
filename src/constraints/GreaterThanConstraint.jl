mutable struct GreaterThanConstraint <: Constraint
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int,Int}
    dual::Union{Value,Nothing}

    function GreaterThanConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
        if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
            error(
                "Cannot create inequality constraint between expressions of sign $(sign(lhs)) and $(sign(rhs))",
            )
        end
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
                "Cannot create inequality constraint between expressions of size $(lhs.size) and $(rhs.size)",
            )
        end
        return new(lhs, rhs, sz, nothing)
    end
end

head(io::IO, ::GreaterThanConstraint) = print(io, "≥")

function vexity(c::GreaterThanConstraint)
    vex = -vexity(c.lhs) + (vexity(c.rhs))
    if vex == ConcaveVexity()
        return NotDcp()
    end
    return vex
end

function _add_constraint!(
    context::Context{T},
    c::GreaterThanConstraint,
) where {T}
    f = conic_form!(context, c.lhs - c.rhs)
    if f isa AbstractVector
        if !all(f .>= -CONSTANT_CONSTRAINT_TOL[])
            @warn "Constant constraint is violated"
            context.detected_infeasible_during_formulation[] = true
        end
        return
    end
    set = MOI.Nonnegatives(MOI.output_dimension(f))
    context.constr_to_moi_inds[c] = MOI_add_constraint(context.model, f, set)
    return
end

Base.:>=(lhs::AbstractExpr, rhs::AbstractExpr) = GreaterThanConstraint(lhs, rhs)

Base.:>=(lhs::AbstractExpr, rhs::Value) = >=(lhs, constant(rhs))

Base.:>=(lhs::Value, rhs::AbstractExpr) = >=(constant(lhs), rhs)

function populate_dual!(model::MOI.ModelLike, c::GreaterThanConstraint, indices)
    ret = MOI.get(model, MOI.ConstraintDual(), indices)
    c.dual = output(reshape(ret, c.size))
    return
end
