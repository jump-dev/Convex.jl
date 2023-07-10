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

# this is kind of awful; currently needed since sometimes Convex treats numbers as vectors or matrices
# need to think of a better way, possibly via breaking changes to Convex's syntax
function promote_size(values)
    d = only(
        unique(
            MOI.output_dimension(v) for
            v in values if v isa MOI.AbstractFunction ||
            v isa VectorAffineFunctionAsMatrix ||
            v isa VAFTape ||
            v isa SparseTape
        ),
    )
    return (v isa Number ? fill(v, d) : v for v in values)
end

scalar_fn(x::Number) = x # for `satisfy` problems? Not sure...
scalar_fn(x) = only(MOIU.scalarize(x))
scalar_fn(x::SparseTape) = scalar_fn(to_vaf(x))
scalar_fn(v::MOI.AbstractScalarFunction) = v

"""
    latex_formulation(problem::Problem, optimizer=MOIU.Model{Float64}())

Prints a LaTeX formulation of the problem. Optionally, pass an `optimizer`
(like `SCS.Optimizer`) as the second argument, to see the formulation provided
for that specific optimizer (since MathOptInterface will reformulate
the problem based on what the optimizer can support).

Uses `MathOptInterface.Utilities.latex_formulation`.
"""
function latex_formulation(problem::Problem, optimizer = MOIU.Model{Float64}())
    context = Context(problem, optimizer)
    return MOIU.latex_formulation(context.model)
end

function solve!(problem::Problem{T}, optimizer_factory; kwargs...) where {T}
    optimizer = MOI.instantiate(optimizer_factory)
    return solve!(problem, optimizer; kwargs...)
end

function solve!(
    p::Problem{T},
    optimizer::MOI.ModelLike;
    silent_solver = false,
) where {T}
    context = Context(p, optimizer)
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
            vectorized_value = MOI.get(model, MOI.VariablePrimal(), var_indices)
            set_value!(
                var,
                unpackvec(vectorized_value, size(var), iscomplex(sign(var))),
            )
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

# In `packvec`, we use `real` even in the `iscomplex == false` branch
# for type stability. Note here `iscomplex` refers to the `Sign` of the variable
# associated to the values here, not whether or not the variable actually holds
# a complex number at this point. Hence, one could have `iscomplex==true` and
# `isreal(value)==true` or even `value isa Real` could hold.
function packvec(value::Number, iscomplex::Bool)
    return iscomplex ? [real(value), imag(value)] : [real(value)]
end

function packvec(value::AbstractArray, iscomplex::Bool)
    value = reshape(value, length(value))
    if iscomplex
        return [real(value); imag(value)]
    else
        return real(value)
    end
end
