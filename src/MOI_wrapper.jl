# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

struct Optimizer{T,M} <: MOI.AbstractOptimizer
    context::Context{T,M}
    moi_to_convex::OrderedCollections.OrderedDict{MOI.VariableIndex,Any}
    convex_to_moi::Dict{UInt64,Vector{MOI.VariableIndex}}
    constraint_map::Vector{MOI.ConstraintIndex}

    function Optimizer(context::Context{T,M}) where {T,M}
        return new{T,M}(
            context,
            OrderedCollections.OrderedDict{MOI.VariableIndex,UInt64}(),
            Dict{UInt64,MOI.VariableIndex}(),
            MOI.ConstraintIndex[],
        )
    end
    function Optimizer{T}(optimizer_constructor) where {T}
        # See https://github.com/jump-dev/Convex.jl/issues/564 for details
        return Optimizer(Context{T}(optimizer_constructor; add_cache = true))
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
    model.moi_to_convex[vi] = var
    model.context.var_to_moi_indices[var] = [vi]
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
    ::Optimizer,
    ::Type{MOI.VectorNonlinearFunction},
    ::Type{<:MOI.AbstractVectorSet},
)
    # This can cause false positives because:
    #  1) some sets might not be supported by Convex.jl (e.g., `vexity` might
    #     be missing)
    #  2) whether we support the constraint can depend on the vexity of the
    #     function, which we currently don't know.
    # Rather than attempt an enumeration of supported sets here, let's just
    # pass things on and hope that there is a nice error message elsewhere in
    # the callchain.
    return true
end

function _expr(::Optimizer, v::Value)
    return Constant(v)
end

function _expr(model::Optimizer, x::MOI.VariableIndex)
    return model.moi_to_convex[x]
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

function _expr(model::Optimizer, f::MOI.VectorNonlinearFunction)
    return vcat(_expr.(model, f.rows)...)
end

function MOI.get(::Optimizer, ::MOI.ListOfSupportedNonlinearOperators)
    return Symbol[:+, :-, :*, :/, :^, :min, :max, :abs, :sqrt, :exp, :log]
end

function MOI.add_constraint(
    model::Optimizer{T},
    func::MOI.VectorNonlinearFunction,
    set::MOI.AbstractVectorSet,
) where {T}
    constraint = Constraint(_expr(model, func), set)
    if vexity(constraint) == Convex.NotDcp()
        msg = "\n\n[Convex.jl] The constraint is not convex according to the axioms of Disciplined Convex Programming.\n\n"
        throw(MOI.AddConstraintNotAllowed{typeof(func),typeof(set)}(msg))
    end
    add_constraint!(model.context, constraint)
    push!(model.constraint_map, model.context.constr_to_moi_inds[constraint])
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(
        length(model.constraint_map),
    )
end

function MOI.is_valid(model::Optimizer, i::MOI.Index)
    return MOI.is_valid(model.context.model, i)
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction},
)
    return true
end

function MOI.set(
    model::Optimizer{T},
    attr::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction},
    func::MOI.ScalarNonlinearFunction,
) where {T}
    obj_fn = _expr(model, func)
    vex = vexity(obj_fn)
    if MOI.get(model, MOI.ObjectiveSense()) == MOI.MAX_SENSE
        vex = -vex
    end
    if vex in (Convex.NotDcp(), Convex.ConcaveVexity())
        msg = "\n\n[Convex.jl] The objective is not convex according to the axioms of Disciplined Convex Programming.\n\n"
        throw(MOI.SetAttributeNotAllowed(attr, msg))
    end
    cfp = conic_form!(model.context, obj_fn)
    obj = _to_scalar_moi(T, cfp)
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

function MOI.get(
    model::Optimizer,
    attr::Union{MOI.ConstraintDual,MOI.ConstraintPrimal},
    ci::MOI.ConstraintIndex{MOI.VectorNonlinearFunction,S},
) where {S<:MOI.AbstractVectorSet}
    return MOI.get(model.context.model, attr, model.constraint_map[ci.value])
end

function MOI.get(model::Optimizer, I::Type{<:MOI.Index}, name::String)
    return MOI.get(model.context.model, I, name)
end
