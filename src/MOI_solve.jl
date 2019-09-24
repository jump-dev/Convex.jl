using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

export solve!

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
        error("Cone $cone not found somehow")
    end
    return set
end


function add_terms!(terms::Vector{MOI.VectorAffineTerm{T}}, matrix::SparseMatrixCSC, rowinds, colinds, vars) where T
    CIs = CartesianIndices((rowinds, colinds))
    I, J, V = findnz(matrix)
    for  k = eachindex(I, J, V)
        term = MOI.VectorAffineTerm{T}(CIs[I[k],J[k]][1] , MOI.ScalarAffineTerm{T}(V[k], vars[CIs[I[k],J[k]][2]]))
        push!(terms, term)
    end
end

function add_terms!(terms::Vector{MOI.VectorAffineTerm{T}}, vector::SparseVector, rowinds, colinds, vars) where T
    CIs = CartesianIndices((rowinds, colinds))
    I, V = findnz(vector)
    for k = eachindex(I, V)
        term =  MOI.VectorAffineTerm{T}(CIs[I[k]][1] , MOI.ScalarAffineTerm{T}(V[k], vars[CIs[I[k]][2]]))
        push!(terms, term)
    end
end

function add_terms!(terms::Vector{MOI.VectorAffineTerm{T}}, matrix::Matrix, rowinds, colinds, vars) where T
    CIs = CartesianIndices((rowinds, colinds))
    for ci = CartesianIndices(size(matrix))
        term = MOI.VectorAffineTerm{T}(CIs[ci][1] , MOI.ScalarAffineTerm{T}(matrix[ci], vars[CIs[ci][2]]))
        push!(terms, term)
    end
end

function add_terms!(terms::Vector{MOI.VectorAffineTerm{T}}, vector::Vector, rowinds, colinds, vars) where T
    CIs = CartesianIndices((rowinds, colinds))
    for k = 1:length(vector)
        term = MOI.VectorAffineTerm{T}(CIs[k][1] , MOI.ScalarAffineTerm{T}(vector[k], vars[CIs[k][2]]))
        push!(terms, term)
    end
end



function process_constr!(constr_fns, sets, T, constraint, var_to_ranges, vars, id_to_variables)
        total_constraint_size = sum(constraint.sizes)
        constr_index = 0
        # A = spzeros(T, total_constraint_size, var_size)
        b = spzeros(T, total_constraint_size)
        terms =  MOI.VectorAffineTerm{T}[]
        for i = 1:length(constraint.objs)
            sz = constraint.sizes[i]
            for (id, val) in constraint.objs[i]
                if id == objectid(:constant)
                    for l in 1:sz
                        # b[constr_index + l] = val[1][l] == 0 ? val[2][l] : val[1][l]
                        b[constr_index + l] = ifelse(val[1][l] == 0, val[2][l],  val[1][l])
                    end
                else
                    var_range = var_to_ranges[id]
                    if id_to_variables[id].sign == ComplexSign()
                        # A[constr_index + 1 : constr_index + sz, var_range[1] : var_range[1] + length(id_to_variables[id])-1] = -val[1]
                        add_terms!(terms, val[1], constr_index + 1 : constr_index + sz, var_range[1] : var_range[1] + length(id_to_variables[id])-1, vars)

                        # A[constr_index + 1 : constr_index + sz, var_range[1] + length(id_to_variables[id]) : var_range[2]] = -val[2]
                        add_terms!(terms, val[2], constr_index + 1 : constr_index + sz, var_range[1] + length(id_to_variables[id]) : var_range[2], vars)

                    else
                        # A[constr_index + 1 : constr_index + sz, var_range[1] : var_range[2]] = -val[1]
                        add_terms!(terms, val[1], constr_index + 1 : constr_index + sz,var_range[1] : var_range[2], vars)

                    end
                end
            end
            constr_index += sz
        end
        push!(sets, get_MOI_set(constraint.cone, total_constraint_size))
        push!(constr_fns, MOI.VectorAffineFunction{T}(terms, b))

end

function load_MOI_model!(model, problem::Problem{T}) where {T}
    if length(problem.objective) != 1
        error("Objective must be a scalar")
    end

    # conic problems have the form
    # minimize c'*x
    # st       b - Ax \in cones
    # our job is to take the conic forms of the objective and constraints
    # and convert them into vectors b and c and a matrix A
    # one chunk of rows in b and in A corresponds to each constraint,
    # and one chunk of columns in b and A corresponds to each variable,
    # with the size of the chunk determined by the size of the constraint or of the variable

    unique_conic_forms = UniqueConicForms()
    objective, objective_var_id = conic_form!(problem, unique_conic_forms)
    constraints = unique_conic_forms.constr_list
    conic_constr_to_constr = unique_conic_forms.conic_constr_to_constr
    id_to_variables = unique_conic_forms.id_to_variables

    # var_to_ranges maps from variable id to the (start_index, stop_index) pairs of the columns of A corresponding to that variable
    # var_size is the sum of the lengths of all variables in the problem
    # constr_size is the sum of the lengths of all constraints in the problem
    var_size, constr_size, var_to_ranges = find_variable_ranges(constraints, id_to_variables)
    
    # c = spzeros(T, var_size, 1)
    objective_range = var_to_ranges[objective_var_id]

    vars = MOI.add_variables(model, var_size)

    # the objective
    val = problem.head == :maximize ? -one(T) : one(T)
    obj_terms = [  MOI.ScalarAffineTerm(val, vars[k]) for k in objective_range[1]:objective_range[2] ]
    obj = MOI.ScalarAffineFunction{T}(obj_terms, zero(T))
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), obj)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    MOI_constr_fn = Union{MOI.VectorAffineFunction{T}, MOI.SingleVariable}[]
    MOI_sets = Any[]
    for constraint in constraints
        process_constr!(MOI_constr_fn, MOI_sets, T, constraint, var_to_ranges, vars, id_to_variables)
    end

    # find integral and boolean variables
    for var_id in keys(var_to_ranges)
        variable = id_to_variables[var_id]
        if :Int in variable.sets
            startidx, endidx = var_to_ranges[var_id]
            for idx in startidx:endidx
                push!(MOI_constr_fn, MOI.SingleVariable(vars[idx]))
                push!(MOI_sets, MOI.Integer())
            end
        end
        if :Bin in variable.sets
            startidx, endidx = var_to_ranges[var_id]
            for idx in startidx:endidx
                push!(MOI_constr_fn, MOI.SingleVariable(vars[idx]))
                push!(MOI_sets, MOI.ZeroOne())
            end
        end
    end

    constraint_indices = MOI.add_constraints(model, MOI_constr_fn, MOI_sets)

    return var_to_ranges, id_to_variables, conic_constr_to_constr, constraints, constraint_indices
end


function solve!(problem::Problem{T}, optimizer::MOI.ModelLike;
    check_vexity = true,
    verbose = true) where {T}

    if check_vexity
        vex = vexity(problem)
    end

    model = MOIU.CachingOptimizer(MOIU.Model{T}(), MOIU.MANUAL)
    var_to_ranges, id_to_variables, conic_constr_to_constr, constraints, constraint_indices = load_MOI_model!(model, problem)

    universal_fallback = MOIU.UniversalFallback(MOIU.Model{T}())
    optimizer = MOIU.CachingOptimizer(universal_fallback, optimizer)
    optimizer = MOI.Bridges.full_bridge_optimizer(optimizer, T)

    MOIU.reset_optimizer(model, optimizer);
    MOIU.attach_optimizer(model);
    MOI.optimize!(model)
    problem.model = model

    # # populate the status, the primal (and possibly dual) solution
    # # and the primal (and possibly dual) variables with values
    moi_populate_solution!(model, problem, var_to_ranges, id_to_variables, conic_constr_to_constr, constraints, constraint_indices)
    if problem.status != MOI.OPTIMAL && verbose
        @warn "Problem status $(problem.status); solution may be inaccurate."
    end

end

function moi_populate_solution!(model::MOI.ModelLike, problem, var_to_ranges, id_to_variables, conic_constr_to_constr, constraints, constraint_indices)
    status = MOI.get(model, MOI.TerminationStatus())
    dual_status = MOI.get(model, MOI.DualStatus())
    primal_status = MOI.get(model, MOI.PrimalStatus())

    # should check when this is allowed
    objective = MOI.get(model, MOI.ObjectiveValue())

    if primal_status != MOI.NO_SOLUTION
        vars = MOI.get(model, MOI.ListOfVariableIndices())
        primal = MOI.get(model, MOI.VariablePrimal(), vars)
    else
        primal = fill(NaN, MOI.get(model, MOI.NumberOfVariables()))
    end

    if dual_status != MOI.NO_SOLUTION
        for (idx, constr) in enumerate(constraints)
            haskey(conic_constr_to_constr, constr) || continue
            MOI_constr_idx = constraint_indices[idx]
            dual_value = MOI.get(model, MOI.ConstraintDual(), MOI_constr_idx)
            if length(dual_value) == 1
                dual_value = dual_value[]
            end
            conic_constr_to_constr[constr].dual = dual_value
        end
    end


    problem.optval = problem.head == :maximize ? -objective : objective
    problem.status = status

    populate_variables_moi!(primal, var_to_ranges, id_to_variables)

end


# this is somehow working! But it's the same as MBP version even though
# MathOptInterface packs them differently than MathProgBase
# todo: figure out why this works still
function populate_variables_moi!(primal, var_to_ranges::Dict{UInt64,Tuple{Int,Int}}, id_to_variables)
    for (id, (start_index, end_index)) in var_to_ranges
        var = id_to_variables[id]
        sz = var.size
        if var.sign != ComplexSign()
            var.value = reshape(primal[start_index:end_index], sz[1], sz[2])
            if sz == (1, 1)
                var.value = var.value[1]
            end
        else
            real_value = reshape(primal[start_index:start_index + div(end_index - start_index + 1, 2) - 1], sz[1], sz[2])
            imag_value = reshape(primal[start_index + div(end_index - start_index + 1, 2):end_index], sz[1], sz[2])
            var.value = real_value + im * imag_value
            if sz == (1, 1)
                var.value = var.value[1]
            end
        end
    end
end
