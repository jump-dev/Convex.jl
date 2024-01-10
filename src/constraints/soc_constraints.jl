mutable struct SOCConstraint <: Constraint
    children::Tuple
    dual::Union{Value,Nothing}

    SOCConstraint(args::AbstractExpr...) = new(args, nothing)
end

head(io::IO, ::SOCConstraint) = print(io, "soc")

function _add_constraint!(context::Context{T}, c::SOCConstraint) where {T}
    f = operate(
        vcat,
        T,
        sum(map(sign, c.children)),
        map(child -> conic_form!(context, child), c.children)...,
    )
    context.constr_to_moi_inds[c] = MOI_add_constraint(
        context.model,
        f,
        MOI.SecondOrderCone(MOI.output_dimension(f)),
    )
    return
end

function populate_dual!(model::MOI.ModelLike, c::SOCConstraint, indices)
    c.dual = output(MOI.get(model, MOI.ConstraintDual(), indices))
    return
end
