function add_variables!(model, var::AbstractVariable)
    if sign(var) == ComplexSign()
        L = MOI.add_variables(model, length(var))
        R = MOI.add_variables(model, length(var))
        return (L, R)
    end
    return MOI.add_variables(model, length(var))
end

"""
    solve!(
        problem::Problem,
        optimizer_factory;
        silent_solver = false,
    )

Solves the problem, populating `problem.optval` with the optimal value, as well
as the values of the variables (accessed by [`evaluate`](@ref)) and constraint
duals (accessed by `cons.dual`), where applicable.

Optional keyword arguments:

 * `silent_solver`: whether the solver should be silent (and not emit output or
   logs) during the solution process.
"""
function solve!(p::Problem, optimizer_factory; silent_solver = false)
    if problem_vexity(p) in (ConcaveVexity(), NotDcp())
        throw(DCPViolationError())
    end
    context = Context(p, optimizer_factory)
    if silent_solver
        MOI.set(context.model, MOI.Silent(), true)
    end
    if context.detected_infeasible_during_formulation[]
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
    for (id, indices) in context.var_id_to_moi_indices
        var = context.id_to_variables[id]
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
