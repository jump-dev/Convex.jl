function add_variables!(model, var::AbstractVariable)
    var.id_hash == objectid(:constant) &&
        error("Internal error: constant used as variable")
    return if sign(var) == ComplexSign()
        L = MOI.add_variables(model, length(var))
        R = MOI.add_variables(model, length(var))
        return (L, R)
    else
        return MOI.add_variables(model, length(var))
    end
end

scalar_fn(x::Number) = x # for `satisfy` problems? Not sure...
scalar_fn(x) = only(MOI.Utilities.scalarize(x))
scalar_fn(x::SparseTape) = scalar_fn(to_vaf(x))
scalar_fn(v::MOI.AbstractScalarFunction) = v

"""
    latex_formulation(problem::Problem, optimizer=MOI.Utilities.Model{Float64}())

Prints a LaTeX formulation of the problem. Optionally, pass an `optimizer`
(like `SCS.Optimizer`) as the second argument, to see the formulation provided
for that specific optimizer (since MathOptInterface will reformulate
the problem based on what the optimizer can support).

Uses `MathOptInterface.Utilities.latex_formulation`.
"""
function latex_formulation(problem::Problem, optimizer = MOI.Utilities.Model{Float64}())
    context = Context(problem, optimizer)
    return MOI.Utilities.latex_formulation(context.model)
end

"""
    solve!(problem::Problem, optimizer_factory;
        silent_solver = false,
    )

Solves the problem, populating `problem.optval` with the optimal value,
as well as the values of the variables (accessed by [`evaluate`](@ref))
and constraint duals (accessed by `cons.dual`), where applicable.
Optional keyword arguments:
* `silent_solver`: whether the solver should be silent (and not emit output or logs) during the solution process.
"""
function solve!(p::Problem, optimizer_factory; silent_solver = false)
    context = Context(p, optimizer_factory)
    model = context.model

    if silent_solver
        MOI.set(model, MOI.Silent(), true)
    end

    if context.detected_infeasible_during_formulation[]
        p.status = MOI.INFEASIBLE
    else
        MOI.optimize!(model)
        p.status = MOI.get(model, MOI.TerminationStatus())
    end

    p.model = model

    if p.status != MOI.OPTIMAL
        @warn "Problem wasn't solved optimally" status = p.status
    end

    dual_status = MOI.get(model, MOI.DualStatus())
    primal_status = MOI.get(model, MOI.PrimalStatus())

    var_to_indices = context.var_id_to_moi_indices
    id_to_variables = context.id_to_variables

    if primal_status != MOI.NO_SOLUTION
        for (id, var_indices) in var_to_indices
            var = id_to_variables[id]
            vexity(var) == ConstVexity() && continue
            if var_indices isa Tuple
                vectorized_value_re =
                    MOI.get(model, MOI.VariablePrimal(), var_indices[1])
                vectorized_value_im =
                    MOI.get(model, MOI.VariablePrimal(), var_indices[2])
                set_value!(
                    var,
                    unpackvec(vectorized_value_re, size(var), false) +
                    im * unpackvec(vectorized_value_im, size(var), false),
                )
            else
                vectorized_value =
                    MOI.get(model, MOI.VariablePrimal(), var_indices)
                set_value!(
                    var,
                    unpackvec(
                        vectorized_value,
                        size(var),
                        iscomplex(sign(var)),
                    ),
                )
            end
        end
    else
        for (id, var_indices) in var_to_indices
            var = id_to_variables[id]
            vexity(var) == ConstVexity() && continue
            set_value!(var, nothing)
        end
    end

    if dual_status != MOI.NO_SOLUTION
        for (constr, MOI_constr_indices) in pairs(context.constr_to_moi_inds)
            populate_dual!(model, constr, MOI_constr_indices)
        end
    else
        # Empty duals
        for constr in keys(context.constr_to_moi_inds)
            constr.dual = nothing
        end
    end

    return nothing
end

function populate_dual!(model::MOI.ModelLike, constr::Constraint, ::Nothing)
    constr.dual = nothing
    return nothing
end

# This is type unstable!
function unpackvec(v::AbstractVector, size::Tuple{Int,Int}, iscomplex::Bool)
    if iscomplex && length(v) == 2
        return v[1] + im * v[2]
    elseif iscomplex
        l = length(v) ÷ 2
        # use views?
        return reshape(v[1:l], size) + im * reshape(v[l+1:end], size)
    elseif length(v) == 1
        return v[]
    else
        return reshape(v, size)
    end
end
