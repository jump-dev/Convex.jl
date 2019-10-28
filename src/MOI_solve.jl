using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

export solve!

# Convert from sets used within Convex to MOI sets
function get_MOI_set(cone, length_inds)
    if cone == :SDP
        set = MOI.PositiveSemidefiniteConeTriangle(Int(sqrt(.25 + 2 * length_inds) - .5))
    elseif cone == :Zero
        set = MOI.Zeros(length_inds)
    elseif cone == :Free
        set = MOI.Reals(length_inds)
    elseif cone == :NonNeg
        set = MOI.Nonnegatives(length_inds)
    elseif cone == :NonPos
        set = MOI.Nonpositives(length_inds)
    elseif cone == :SOC
        set = MOI.SecondOrderCone(length_inds)
    elseif cone == :SOCRotated
        set = MOI.RotatedSecondOrderCone(length_inds)
    elseif cone == :ExpPrimal
        set = MOI.ExponentialCone()
    elseif cone == :ExpDual
        set = MOI.DualExponentialCone()
    else
        error("Cone $cone not found somehow; please file an issue.")
    end
    return set
end

# add_terms! turns a matrix or vector `X` into a collection of `VectorAffineTerm`'s, which it adds to a list (the first argument).
# This replicates the behavior of the relation `A[row_indices, col_indices] = X`, where `X` could be a vector or matrix.
# Here, `row_indices` corresponds to the output index of a `VectorAffineTerm`
# Instead of `col_indices`, here we have their wrapped counterparts, `var_indices` (a vector of `MOI.VariableIndex`'s).
function add_terms!(terms::Vector{MOI.VectorAffineTerm{T}}, matrix::SparseMatrixCSC, row_indices, var_indices) where T
    CIs = CartesianIndices((row_indices, 1:length(var_indices)))
    I, J, V = findnz(matrix)
    for  k = eachindex(I, J, V)
        term = MOI.VectorAffineTerm{T}(CIs[I[k],J[k]][1], MOI.ScalarAffineTerm{T}(V[k], var_indices[CIs[I[k],J[k]][2]]))
        push!(terms, term)
    end
end

function add_terms!(terms::Vector{MOI.VectorAffineTerm{T}}, vector::SparseVector, row_indices, var_indices) where T
    CIs = CartesianIndices((row_indices, 1:length(var_indices)))
    I, V = findnz(vector)
    for k = eachindex(I, V)
        term =  MOI.VectorAffineTerm{T}(CIs[I[k]][1], MOI.ScalarAffineTerm{T}(V[k], var_indices[CIs[I[k]][2]]))
        push!(terms, term)
    end
end

function add_terms!(terms::Vector{MOI.VectorAffineTerm{T}}, matrix::Matrix, row_indices, var_indices) where T
    CIs = CartesianIndices((row_indices, 1:length(var_indices)))
    for ci = CartesianIndices(size(matrix))
        term = MOI.VectorAffineTerm{T}(CIs[ci][1], MOI.ScalarAffineTerm{T}(matrix[ci], var_indices[CIs[ci][2]]))
        push!(terms, term)
    end
end

function add_terms!(terms::Vector{MOI.VectorAffineTerm{T}}, vector::Vector, row_indices, var_indices) where T
    CIs = CartesianIndices((row_indices, 1:length(var_indices)))
    for k = 1:length(vector)
        term = MOI.VectorAffineTerm{T}(CIs[k][1], MOI.ScalarAffineTerm{T}(vector[k], var_indices[CIs[k][2]]))
        push!(terms, term)
    end
end


# This function generates a MOI set and VectorAffineFunction corresponding to a constraint that
# Convex represents with a `ConicConstr` object.
# This is part of what used to be `Convex.conic_problem` (under MPB instead of MOI).
# The constraint is of the form Ax + b ∈ set
# where `A` is a matrix, `x` a vector of variables, and `b` a constant vector.
# MOI represents this type of constraint as a VectorAffineFunction, where `A` is represented
# by a collection of `VectorAffineTerm`'s.
function make_MOI_constr(conic_constr::ConicConstr, var_to_indices, id_to_variables, T)
    total_constraint_size = sum(conic_constr.sizes)
    constr_index = 0
    # Under MPB, we would have
    # A = spzeros(T, total_constraint_size, var_size)
    # Instead, we hold a vector of `VectorAffineTerm`'s.
    terms =  MOI.VectorAffineTerm{T}[]
    b = spzeros(T, total_constraint_size)

    for i = 1:length(conic_constr.objs)
        sz = conic_constr.sizes[i]
        for (id, val) in conic_constr.objs[i]
            if id == objectid(:constant)
                for l in 1:sz
                    b[constr_index + l] = ifelse(val[1][l] == 0, val[2][l],  val[1][l])
                end
            else
                var_indices = var_to_indices[id]
                if id_to_variables[id].sign == ComplexSign()
                    l = length(var_indices) ÷ 2
                    # Under MPB, we would perform the following to update `A`
                    # A[constr_index + 1 : constr_index + sz, var_range[1] : var_range[1] + length(id_to_variables[id])-1] = -val[1]
                    # where `var_range` is the analog of `var_indices` (i.e. unwrapped indices).
                    # Here, instead we use `add_terms!` to reproduce this behavior with MOI's `VectorAffineTerm`s.
                    # We also differ by a minus sign due a difference in convention between MPB and MOI.
                    add_terms!(terms, val[1], constr_index + 1 : constr_index + sz, var_indices[1:l])

                    # MPB version:
                    # A[constr_index + 1 : constr_index + sz, var_range[1] + length(id_to_variables[id]) : var_range[2]] = -val[2]
                    add_terms!(terms, val[2], constr_index + 1 : constr_index + sz, var_indices[l+1 : end])

                else
                    # MPB version:
                    # A[constr_index + 1 : constr_index + sz, var_range[1] : var_range[2]] = -val[1]
                    add_terms!(terms, val[1], constr_index + 1 : constr_index + sz, var_indices)

                end
            end
        end
        constr_index += sz
    end
    set = get_MOI_set(conic_constr.cone, total_constraint_size)
    constr_fn =  MOI.VectorAffineFunction{T}(terms, b)
    return set, constr_fn
end


# Under MPB, this was the function `find_variable_ranges`.
# This function gets MOI variable indices for each variable used as a key
# inside a `ConicObj` used in the problem.
function get_variable_indices!(model, conic_constraints, id_to_variables)
    var_to_indices = Dict{UInt64,Vector{MOI.VariableIndex}}()
    for conic_constr in conic_constraints
        for i = 1:length(conic_constr.objs)
            for (id, val) in conic_constr.objs[i]
                if !haskey(var_to_indices, id) && id != objectid(:constant)
                    var = id_to_variables[id]
                    if var.sign == ComplexSign()
                        var_size = 2 * length(var)
                    else
                        var_size = length(var)
                    end
                    var_to_indices[id] = MOI.add_variables(model, var_size)
                end
            end
        end
    end
    return var_to_indices
end



# This function traverses the problem tree and loads
# an MOI model with the corresponding problem instance.
# The problem is written in a simple epigraph form.
# Under MPB, this was `Convex.conic_problem`.
function load_MOI_model!(model, problem::Problem{T}) where {T}
    if length(problem.objective) != 1
        error("Objective must be a scalar")
    end

    unique_conic_forms = UniqueConicForms()
    objective, objective_var_id = conic_form!(problem, unique_conic_forms)
    conic_constraints = unique_conic_forms.constr_list
    conic_constr_to_constr = unique_conic_forms.conic_constr_to_constr
    id_to_variables = unique_conic_forms.id_to_variables

    # var_to_indices maps from variable id to the MOI `VariableIndex`'s corresponding to the variable
    var_to_indices = get_variable_indices!(model, conic_constraints, id_to_variables)

    # the objective: maximize or minimize a scalar variable
    objective_index = var_to_indices[objective_var_id][] # get the `MOI.VariableIndex` corresponding to the objective
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(objective_index))
    MOI.set(model, MOI.ObjectiveSense(), problem.head == :maximize ? MOI.MAX_SENSE : MOI.MIN_SENSE)

    # Constraints: Generate a MOI function and a MOI sets for each `ConicConstr` object in the problem
    MOI_constr_fn = Union{MOI.VectorAffineFunction{T},MOI.SingleVariable}[]
    MOI_sets = Any[]
    for conic_constr in conic_constraints
        set, constr_fn = make_MOI_constr(conic_constr, var_to_indices, id_to_variables, T)
        push!(MOI_sets, set)
        push!(MOI_constr_fn, constr_fn)
    end

    # Add integral and boolean constraints
    for var_id in keys(var_to_indices)
        variable = id_to_variables[var_id]
        if :Int in variable.sets
            var_indices = var_to_indices[var_id]
            for idx = eachindex(var_indices)
                push!(MOI_constr_fn, MOI.SingleVariable(var_indices[idx]))
                push!(MOI_sets, MOI.Integer())
            end
        end
        if :Bin in variable.sets
            var_indices = var_to_indices[var_id]
            for idx in eachindex(var_indices)
                push!(MOI_constr_fn, MOI.SingleVariable(var_indices[idx]))
                push!(MOI_sets, MOI.ZeroOne())
            end
        end
    end

    # Add all the constraints to the model and collect the corresponding MOI indices
    constraint_indices = MOI.add_constraints(model, MOI_constr_fn, MOI_sets)

    return id_to_variables, conic_constr_to_constr, conic_constraints, var_to_indices, constraint_indices
end


function solve!(problem::Problem{T}, optimizer::MOI.ModelLike;
    check_vexity = true,
    verbose = true,
    warmstart = false) where {T}

    if check_vexity
        vex = vexity(problem)
    end

    model = MOIB.full_bridge_optimizer(
        MOIU.CachingOptimizer(
            MOIU.UniversalFallback(MOIU.Model{T}()),
            optimizer
        ),
        T
    )
    
    id_to_variables, conic_constr_to_constr, conic_constraints, var_to_indices, constraint_indices = load_MOI_model!(model, problem)

    if warmstart
        warmstart_variables!(model, var_to_indices, id_to_variables, T, verbose)
    end
 
    MOI.optimize!(model)
    problem.model = model

    # populate the status, primal variables, and dual variables (when possible)
    moi_populate_solution!(model, problem, id_to_variables, conic_constr_to_constr, conic_constraints, var_to_indices, constraint_indices)
    
    if problem.status != MOI.OPTIMAL && verbose
        @warn "Problem status $(problem.status); solution may be inaccurate."
    end
end


function warmstart_variables!(model, var_to_indices, id_to_variables, T, verbose)
    if MOI.supports(model, MOI.VariablePrimalStart(), MOI.VariableIndex)
        for (id, var_inds) in pairs(var_to_indices)
            x = id_to_variables[id]
            value = x.value
            value === nothing && continue
            value_vec = packvec(value, sign(x) == ComplexSign())
            value_vec = convert(Vector{T}, value_vec)
            MOI.set(model, MOI.VariablePrimalStart(), var_inds, value_vec)
        end
    elseif verbose
        @warn "Skipping variable warmstart; the solver does not support it."
    end
end

function moi_populate_solution!(model::MOI.ModelLike, problem, id_to_variables, conic_constr_to_constr, conic_constraints, var_to_indices, constraint_indices)
    status = MOI.get(model, MOI.TerminationStatus())
    problem.status = status

    # should check when this is allowed
    objective = MOI.get(model, MOI.ObjectiveValue())
    problem.optval = objective

    dual_status = MOI.get(model, MOI.DualStatus())
    primal_status = MOI.get(model, MOI.PrimalStatus())

    if primal_status != MOI.NO_SOLUTION
        for (id, var_indices) in var_to_indices
            var = id_to_variables[id]
            vectorized_value =  MOI.get(model, MOI.VariablePrimal(), var_indices)
            var.value = unpackvec(vectorized_value, size(var), var.sign == ComplexSign())
        end
    else
        for (id, var_indices) in var_to_indices
            var = id_to_variables[id]
            var.value = nothing
        end
    end

    if dual_status != MOI.NO_SOLUTION
        for (idx, conic_constr) in enumerate(conic_constraints)
            haskey(conic_constr_to_constr, conic_constr) || continue
            constr = conic_constr_to_constr[conic_constr]
            MOI_constr_indices = constraint_indices[idx]
            dual_value_vectorized = MOI.get(model, MOI.ConstraintDual(), MOI_constr_indices)
            iscomplex = sign(constr.lhs) == ComplexSign() || sign(constr.rhs) == ComplexSign()
            constr.dual = unpackvec(dual_value_vectorized, constr.size, iscomplex)
        end
    else
        for (idx, conic_constr) in enumerate(conic_constraints)
            haskey(conic_constr_to_constr, conic_constr) || continue
            constr = conic_constr_to_constr[conic_constr]
            constr.dual = nothing
        end
    end

end


# This is type unstable!
function unpackvec(v::AbstractVector, size::Tuple{Int,Int}, iscomplex::Bool)
    if iscomplex && length(v) == 2
        return v[1] + im * v[2]
    elseif iscomplex
        l = length(v) ÷ 2
        # use views?
        return reshape(v[1:l], size) + im * reshape(v[l+1 : end], size)
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
    iscomplex ? [real(value), imag(value)] : [real(value)]
end

function packvec(value::AbstractArray, iscomplex::Bool)
    value = reshape(value, length(value))
    if iscomplex
        return [real(value); imag(value)]
    else
        return real(value) 
    end
end
