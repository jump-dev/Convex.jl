struct Optimizer{T,M} <: MOI.AbstractOptimizer
    context::Context{T,M}
    moi_to_convex::OrderedDict{MOI.VariableIndex,UInt64}
    convex_to_moi::Dict{UInt64,Vector{MOI.VariableIndex}}
    constraint_map::Vector{MOI.ConstraintIndex}
    function Optimizer(context::Context{T,M}) where {T,M}
        return new{T,M}(
            context,
            OrderedDict{MOI.VariableIndex,UInt64}(),
            Dict{UInt64,MOI.VariableIndex}(),
            MOI.ConstraintIndex[],
        )
    end
    function Optimizer{T}(model::MOI.ModelLike) where {T}
        return Optimizer(Context{T}(model))
    end
end

Optimizer(model::MOI.ModelLike) = Optimizer{Float64}(model)

MOI.is_empty(model::Optimizer) = MOI.is_empty(model.context.model)

function MOI.empty!(model::Optimizer)
    empty!(model.context)
    empty!(model.moi_to_convex)
    empty!(model.convex_to_moi)
    empty!(model.constraint_map)
    return
end

function MOI.get(model::Optimizer, attr::MOI.SolverName)
    inner = MOI.get(model.context.model, attr)
    return "Convex with $inner"
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

function _add_variable(model::Optimizer, vi::MOI.VariableIndex)
    var = Variable()
    model.moi_to_convex[vi] = var.id_hash
    model.context.var_id_to_moi_indices[var.id_hash] = [vi]
    model.context.id_to_variables[var.id_hash] = var
    return
end

function MOI.supports_add_constrained_variables(
    model::Optimizer,
    S::Type{MOI.Reals},
)
    return MOI.supports_add_constrained_variables(model.context.model, S)
end

function MOI.supports_add_constrained_variables(
    model::Optimizer,
    S::Type{<:MOI.AbstractVectorSet},
)
    return MOI.supports_add_constrained_variables(model.context.model, S)
end

function MOI.add_constrained_variables(
    model::Optimizer,
    set::MOI.AbstractVectorSet,
)
    vis, ci = MOI.add_constrained_variables(model.context.model, set)
    for vi in vis
        _add_variable(model, vi)
    end
    return vis, ci
end

function MOI.supports_add_constrained_variable(
    model::Optimizer,
    S::Type{<:MOI.AbstractScalarSet},
)
    return MOI.supports_add_constrained_variable(model.context.model, S)
end

function MOI.add_constrained_variable(
    model::Optimizer,
    set::MOI.AbstractScalarSet,
)
    vi, ci = MOI.add_constrained_variable(model.context.model, set)
    _add_variable(model, vi)
    return vi, ci
end

function MOI.add_variable(model::Optimizer)
    vi = MOI.add_variable(model.context.model)
    _add_variable(model, vi)
    return vi
end

function MOI.supports_constraint(
    model::Optimizer,
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
)
    return MOI.supports_constraint(model.context.model, F, S)
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.AbstractFunction,
    set::MOI.AbstractSet,
)
    return MOI.add_constraint(model.context.model, func, set)
end

function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:Union{MOI.EqualTo{T},MOI.GreaterThan{T},MOI.LessThan{T}}},
) where {T}
    return true
end

function _expr(vi::MOI.VariableIndex, model)
    return model.context.id_to_variables[model.moi_to_convex[vi]]
end

function _expr(v::Value, _)
    return Constant(v)
end

Base.:+(e::AbstractExpr) = e

function _expr(func::MOI.ScalarNonlinearFunction, model)
    if func.head == :^
        if length(func.args) == 2 && func.args[2] == 2
            return square(_expr(func.args[1], model))
        end
        error(
            "Power with exponent different from 2 is not supported by Convex.jl",
        )
    end
    expr = Expr(:call, func.head, _expr.(func.args, model)...)
    return eval(expr)
end

function _constraint(expr::AbstractExpr, set::MOI.EqualTo)
    return expr == MOI.constant(set)
end

function _constraint(expr::AbstractExpr, set::MOI.LessThan)
    return expr <= MOI.constant(set)
end

function _constraint(expr::AbstractExpr, set::MOI.GreaterThan)
    return expr >= MOI.constant(set)
end

function MOI.add_constraint(
    model::Optimizer{T},
    func::MOI.ScalarNonlinearFunction,
    set::MOI.AbstractScalarSet,
) where {T}
    constraint = _constraint(_expr(func, model), set)
    add_constraint!(model.context, constraint)
    push!(model.constraint_map, model.context.constr_to_moi_inds[constraint])
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(
        length(model.constraint_map),
    )
end

function MOI.is_valid(model::Optimizer, i::MOI.Index)
    return MOI.is_valid(model.context.model, i)
end

MOI.optimize!(model::Optimizer) = MOI.optimize!(model.context.model)

function MOI.supports(
    model::Optimizer,
    attr::Union{MOI.AbstractModelAttribute,MOI.AbstractOptimizerAttribute},
)
    return MOI.supports(model.context.model, attr)
end

function MOI.set(
    model::Optimizer,
    attr::Union{MOI.AbstractModelAttribute,MOI.AbstractOptimizerAttribute},
    value,
)
    return MOI.set(model.context.model, attr, value)
end

function MOI.get(
    model::Optimizer,
    attr::Union{MOI.AbstractModelAttribute,MOI.AbstractOptimizerAttribute},
)
    return MOI.get(model.context.model, attr)
end

function MOI.supports(
    model::Optimizer,
    attr::MOI.AbstractVariableAttribute,
    VI::Type{MOI.VariableIndex},
)
    return MOI.supports(model.context.model, attr, VI)
end

function MOI.set(
    model::Optimizer,
    attr::MOI.AbstractVariableAttribute,
    vi::MOI.VariableIndex,
    value,
)
    return MOI.set(model.context.model, attr, vi, value)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.AbstractVariableAttribute,
    vi::MOI.VariableIndex,
)
    return MOI.get(model.context.model, attr, vi)
end

function MOI.supports(
    model::Optimizer,
    attr::MOI.AbstractConstraintAttribute,
    CI::Type{<:MOI.ConstraintIndex},
)
    return MOI.supports(model.context.model, attr, CI)
end

function MOI.set(
    model::Optimizer,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
    value,
)
    return MOI.set(model.context.model, attr, ci, value)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
)
    return MOI.get(model.context.model, attr, ci)
end

function MOI.get(model::Optimizer, I::Type{<:MOI.Index}, name::String)
    return MOI.get(model.context.model, I, name)
end
