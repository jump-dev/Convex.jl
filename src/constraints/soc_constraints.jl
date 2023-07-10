# TODO: Document this
mutable struct SOCConstraint <: Constraint
    children::Tuple
    dual::ValueOrNothing

    function SOCConstraint(args::AbstractExpr...)
        children = tuple(args...)
        return new(children, nothing)
    end
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

    return nothing
end

function populate_dual!(
    model::MOI.ModelLike,
    constr::SOCConstraint,
    MOI_constr_indices,
)
    return constr.dual =
        output(MOI.get(model, MOI.ConstraintDual(), MOI_constr_indices))
end
