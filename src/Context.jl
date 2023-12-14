mutable struct Context{T,M}
    # MOI model
    model::M

    # Used for populating variable values after solving
    var_id_to_moi_indices::OrderedDict{
        UInt64,
        Union{
            Vector{MOI.VariableIndex},
            Tuple{Vector{MOI.VariableIndex},Vector{MOI.VariableIndex}},
        },
    }
    # `id_hash` -> `AbstractVariable`
    id_to_variables::OrderedDict{UInt64,Any}

    # Used for populating constraint duals
    constr_to_moi_inds::IdDict{Any,Any}

    detected_infeasible_during_formulation::Ref{Bool}

    # Cache
    # conic_form_cache::DataStructures.WeakKeyIdDict{Any, Any}
    conic_form_cache::IdDict{Any,Any}
end

function Context{T}(optimizer_factory) where {T}
    model = MOI.instantiate(optimizer_factory, with_bridge_type = T)
    return Context{T,typeof(model)}(
        model,
        OrderedDict{UInt64,Vector{MOI.VariableIndex}}(),
        OrderedDict{UInt64,Any}(),
        IdDict{Any,Any}(),
        false,
        IdDict{Any,Any}(),
    )
end

function Base.empty!(context::Context)
    MOI.empty!(context.model)
    empty!(context.var_id_to_moi_indices)
    empty!(context.id_to_variables)
    empty!(context.constr_to_moi_inds)
    context.detected_infeasible_during_formulation[] = false
    empty!(context.conic_form_cache)
    return
end
