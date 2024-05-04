# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct Context{T,M}
    # MOI model
    model::M

    # Used for populating variable values after solving
    var_id_to_moi_indices::OrderedCollections.OrderedDict{
        UInt64,
        Union{
            Vector{MOI.VariableIndex},
            Tuple{Vector{MOI.VariableIndex},Vector{MOI.VariableIndex}},
        },
    }
    # `id_hash` -> `AbstractVariable`
    id_to_variables::OrderedCollections.OrderedDict{UInt64,Any}

    # Used for populating constraint duals
    constr_to_moi_inds::IdDict{Any,Any}

    detected_infeasible_during_formulation::Ref{Bool}

    # Cache
    # conic_form_cache::DataStructures.WeakKeyIdDict{Any, Any}
    conic_form_cache::IdDict{Any,Any}
    vaf_cache::IdDict{Any, MOI.VectorAffineFunction{T}} # keys are SparseTape{T}
end

function get_vaf!(context::Context, obj)
    return get!(() -> to_vaf(obj), context.vaf_cache, obj)
end

function Context{T}(optimizer_factory; add_cache::Bool = false) where {T}
    model = if add_cache
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}()),
            MOI.instantiate(optimizer_factory; with_bridge_type = T),
        )
    else
        MOI.instantiate(optimizer_factory; with_bridge_type = T)
    end
    return Context{T,typeof(model)}(
        model,
        OrderedCollections.OrderedDict{UInt64,Vector{MOI.VariableIndex}}(),
        OrderedCollections.OrderedDict{UInt64,Any}(),
        IdDict{Any,Any}(),
        false,
        IdDict{Any,Any}(),
        IdDict{Any, MOI.VectorAffineFunction{T}}()
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
