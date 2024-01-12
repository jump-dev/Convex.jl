mutable struct EqualToConstraint <: Constraint
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int,Int}
    dual::Union{Value,Nothing}

    function EqualToConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
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
                "Cannot create equality constraint between expressions of size $(lhs.size) and $(rhs.size)",
            )
        end
        return new(lhs, rhs, sz, nothing)
    end
end

head(io::IO, ::EqualToConstraint) = print(io, "==")

function vexity(c::EqualToConstraint)
    vex = vexity(c.lhs) + (-vexity(c.rhs))
    # You can't have equality constraints with concave/convex expressions
    if vex == ConvexVexity() || vex == ConcaveVexity()
        return NotDcp()
    end
    return vex
end

function _add_constraint!(context::Context{T}, eq::EqualToConstraint) where {T}
    f = conic_form!(context, eq.lhs - eq.rhs)
    if f isa AbstractVector
        # a trivial constraint without variables like `5 == 0`
        if !all(abs.(f) .<= CONSTANT_CONSTRAINT_TOL[])
            @warn "Constant constraint is violated"
            context.detected_infeasible_during_formulation[] = true
        end
        return
    end
    context.constr_to_moi_inds[eq] =
        MOI_add_constraint(context.model, f, MOI.Zeros(MOI.output_dimension(f)))
    return
end

Base.:(==)(lhs::AbstractExpr, rhs::AbstractExpr) = EqualToConstraint(lhs, rhs)

Base.:(==)(lhs::AbstractExpr, rhs::Value) = ==(lhs, constant(rhs))

Base.:(==)(lhs::Value, rhs::AbstractExpr) = ==(constant(lhs), rhs)

function populate_dual!(model::MOI.ModelLike, c::EqualToConstraint, indices)
    if iscomplex(c)
        re = MOI.get(model, MOI.ConstraintDual(), indices[1])
        imag = MOI.get(model, MOI.ConstraintDual(), indices[2])
        c.dual = output(reshape(re + im * imag, c.size))
    else
        ret = MOI.get(model, MOI.ConstraintDual(), indices)
        c.dual = output(reshape(ret, c.size))
    end
    return
end
