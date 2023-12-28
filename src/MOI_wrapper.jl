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
    function Optimizer{T}(optimizer_constructor) where {T}
        return Optimizer(Context{T}(optimizer_constructor))
    end
end

Optimizer(optimizer_constructor) = Optimizer{Float64}(optimizer_constructor)

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

function _expr(::Optimizer, v::Value)
    return Constant(v)
end

function _expr(model::Optimizer, x::MOI.VariableIndex)
    return model.context.id_to_variables[model.moi_to_convex[x]]
end

function _expr(model::Optimizer, f::MOI.AbstractScalarFunction)
    return _expr(model, convert(MOI.ScalarNonlinearFunction, f))
end

function _expr(model::Optimizer, f::MOI.ScalarNonlinearFunction)
    args = _expr.(model, f.args)
    if f.head == :+
        if length(args) == 1
            return args[1]
        elseif length(args) == 2
            return args[1] + args[2]
        else
            return sum(args)
        end
    elseif f.head == :- && 1 <= length(args) <= 2
        return -(args...)
    elseif f.head == :*
        return *(args...)
    elseif f.head == :/ && length(args) == 2
        if !(f.args[2] isa Real)
            msg = "denominator must be a scalar constant. Got $(f.args[2])"
            throw(MOI.UnsupportedNonlinearOperator(:/, msg))
        end
        return args[1] / f.args[2]
    elseif f.head == :^ && length(args) == 2
        if f.args[2] != 2
            msg = "Power with exponent different from 2 is not supported by Convex.jl"
            throw(MOI.UnsupportedNonlinearOperator(:^, msg))
        end
        return square(args[1])
    elseif f.head == :min
        if length(f.args) == 1
            return args[1]
        elseif length(args) == 2
            return min(args[1], args[2])
        else
            return minimum(args)
        end
    elseif f.head == :max
        if length(f.args) == 1
            return args[1]
        elseif length(args) == 2
            return max(args[1], args[2])
        else
            return maximum(args)
        end
    elseif f.head == :abs
        return abs(args[1])
    elseif f.head == :sqrt
        return sqrt(args[1])
    elseif f.head == :exp
        return exp(args[1])
    elseif f.head == :log
        return log(args[1])
    end
    return throw(MOI.UnsupportedNonlinearOperator(f.head))
end

function MOI.get(::Optimizer, ::MOI.ListOfSupportedNonlinearOperators)
    return Symbol[:+, :-, :*, :/, :^, :min, :max, :abs, :sqrt, :exp, :log]
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
    constraint = _constraint(_expr(model, func), set)
    add_constraint!(model.context, constraint)
    push!(model.constraint_map, model.context.constr_to_moi_inds[constraint])
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(
        length(model.constraint_map),
    )
end

function MOI.is_valid(model::Optimizer, i::MOI.Index)
    return MOI.is_valid(model.context.model, i)
end

function MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction})
    return true
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}, func::MOI.ScalarNonlinearFunction)
    cfp = conic_form!(model.context, _expr(func, model))
    obj = scalar_fn(cfp)
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    return
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
