# TODO: Document this
mutable struct SOCConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    children::Tuple
    dual::ValueOrNothing

    function SOCConstraint(args::AbstractExpr...)
        children = tuple(args...)
        id_hash = hash((children, :soc))
        return new(:soc, id_hash, children, nothing)
    end
end

function _add_constraints_to_context(
    c::SOCConstraint,
    context::Context{T},
) where {T}
    f = operate(
        vcat,
        T,
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
