# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

function add_variables!(model, var::AbstractVariable)
    if sign(var) == ComplexSign()
        L = MOI.add_variables(model, length(var))
        R = MOI.add_variables(model, length(var))
        return (L, R)
    end
    return MOI.add_variables(model, length(var))
end

function _primal_start(
    ::Type{T},
    x::Vector{MOI.VariableIndex},
    start::AbstractArray,
) where {T}
    return convert(Vector{T}, real.(reshape(start, length(x))))
end

function _primal_start(
    ::Type{T},
    x::Tuple{Vector{MOI.VariableIndex},Vector{MOI.VariableIndex}},
    start::AbstractArray,
) where {T}
    y = reshape(start, length(x[1]))
    return convert(Vector{T}, real.(y)), convert(Vector{T}, imag.(y))
end

function _primal_start(::Type{T}, x::Vector, start::Number) where {T}
    return _primal_start(T, x, fill(start, length(x)))
end

function _primal_start(::Type{T}, x::Tuple, start::Number) where {T}
    return _primal_start(T, x, fill(start, length(x[1])))
end

function _add_variable_primal_start(context::Convex.Context{T}) where {T}
    attr = MOI.VariablePrimalStart()
    for (x, moi_indices) in context.var_to_moi_indices
        if x.value === nothing
            continue
        elseif moi_indices isa Tuple  # Variable is complex
            start_re, start_im = _primal_start(T, moi_indices, x.value)
            MOI.set(context.model, attr, moi_indices[1], start_re)
            MOI.set(context.model, attr, moi_indices[2], start_im)
        else
            @assert moi_indices isa Vector{MOI.VariableIndex}
            start = _primal_start(T, moi_indices, x.value)
            MOI.set(context.model, attr, moi_indices, start)
        end
    end
    return
end

"""
    solve!(
        problem::Problem,
        optimizer_factory;
        silent_solver = false,
        warmstart::Bool = true,
    )

Solves the problem, populating `problem.optval` with the optimal value, as well
as the values of the variables (accessed by [`evaluate`](@ref)) and constraint
duals (accessed by `cons.dual`), where applicable.

Optional keyword arguments:

 * `silent_solver`: whether the solver should be silent (and not emit output or
   logs) during the solution process.
 * `warmstart` (default: `false`): whether the solver should start the
   optimization from a previous optimal value (according to the current primal
   value of the variables in the problem, which can be set by [`set_value!`](@ref).
"""
function solve!(
    p::Problem,
    optimizer_factory;
    silent_solver = false,
    warmstart::Bool = false,
)
    if problem_vexity(p) in (ConcaveVexity(), NotDcp())
        throw(DCPViolationError())
    end
    context, (stats...) = @timed Context(p, optimizer_factory)
    if !silent_solver
        s = round(stats.time; digits = 2), Base.format_bytes(stats.bytes)
        @info "[Convex.jl] Compilation finished: $(s[1]) seconds, $(s[2]) of memory allocated"
    end
    if silent_solver
        MOI.set(context.model, MOI.Silent(), true)
    end
    attr = MOI.VariablePrimalStart()
    if warmstart && MOI.supports(context.model, attr, MOI.VariableIndex)
        _add_variable_primal_start(context)
    end
    if context.detected_infeasible_during_formulation
        p.status = MOI.INFEASIBLE
    else
        MOI.optimize!(context.model)
        p.status = MOI.get(context.model, MOI.TerminationStatus())
    end
    p.model = context.model
    if p.status != MOI.OPTIMAL
        @warn "Problem wasn't solved optimally" status = p.status
    end
    primal_status = MOI.get(context.model, MOI.PrimalStatus())
    for (var, indices) in context.var_to_moi_indices
        if vexity(var) == ConstVexity()
            continue  # Fixed variable
        elseif primal_status == MOI.NO_SOLUTION
            set_value!(var, nothing)
        elseif indices isa Tuple  # Complex-valued variable
            primal_re = MOI.get(p.model, MOI.VariablePrimal(), indices[1])
            primal_im = MOI.get(p.model, MOI.VariablePrimal(), indices[2])
            set_value!(
                var,
                _unpack_vector(primal_re, size(var)) +
                im * _unpack_vector(primal_im, size(var)),
            )
        else
            @assert !iscomplex(sign(var))
            primal = MOI.get(p.model, MOI.VariablePrimal(), indices)
            set_value!(var, _unpack_vector(primal, size(var)))
        end
    end
    dual_status = MOI.get(context.model, MOI.DualStatus())
    for (c, indices) in context.constr_to_moi_inds
        if dual_status == MOI.NO_SOLUTION
            c.dual = nothing
        else
            populate_dual!(context.model, c, indices)
        end
    end
    return
end

function populate_dual!(::MOI.ModelLike, c::Constraint, ::Nothing)
    # If the Constraint does not support `populate_dual!` it must return
    # `nothing` from `_add_constraint!`.
    c.dual = nothing
    return
end

function _unpack_vector(v::AbstractVector, size::Tuple{Int,Int})
    if length(v) == 1
        return v[]
    end
    return reshape(v, size)
end
