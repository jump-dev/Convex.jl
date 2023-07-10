struct Context{T,M}
    model::M
    var_id_to_moi_indices::OrderedDict{
        UInt64,
        Union{
            Vector{MOI.VariableIndex},
            Tuple{Vector{MOI.VariableIndex},Vector{MOI.VariableIndex}},
        },
    }
    id_to_variables::OrderedDict{UInt64,Any}

    # Used for populating constraint duals
    constr_to_moi_inds::IdDict{Any,Any}

    detected_infeasible_during_formulation::Ref{Bool}

    # Unused
    # constraint_id_to_moi_index::OrderedDict{UInt64, Any}
    # constraint_id_to_constraint::OrderedDict{UInt64, Any}
end

function Context{T}(optimizer) where {T}
    model = MOIB.full_bridge_optimizer(
        MOIU.CachingOptimizer(
            MOIU.UniversalFallback(MOIU.Model{T}()),
            optimizer,
        ),
        T,
    )
    return Context{T,typeof(model)}(
        model,
        OrderedDict{UInt64,Vector{MOI.VariableIndex}}(),
        OrderedDict{UInt64,Any}(),
        IdDict{Any,Any}(),
        false,
    )
    # OrderedDict{UInt64,Any}(), OrderedDict{UInt64,Any}())
end
