import MathProgBase
export solve!

# semantics:
# * we *load* problem data from Convex into the MathProgBase model
# * we *populate* objects in Convex with their optimal values

function solve!(problem::Problem,
                s::MathProgBase.AbstractMathProgSolver;
                kwargs...)
    # TODO: warn about wiping out old model if warmstart=true?
    problem.model = MathProgBase.ConicModel(s)
    return solve!(problem; kwargs...)
end

function solve!(problem::Problem;
                warmstart=false,
                check_vexity=true,
                verbose=true)

    if problem.model === nothing
        throw(ArgumentError(
            """
            The provided problem hasn't been initialized with a conic model.
            You can resolve this by passing in `AbstractMathProgSolver` such as:
            solve!(problem, ECOSSolver())
            """
        ))
    end

    if check_vexity
        vex = vexity(problem)
    end

    c, A, b, cones, var_to_ranges, vartypes, conic_constraints = conic_problem(problem)

    # load MPB conic problem
    m = problem.model
    load_problem!(m, c, A, b, cones, vartypes)
    if warmstart
        set_warmstart!(m, problem, length(c), var_to_ranges)
    end
    # optimize MPB conic problem
    MathProgBase.optimize!(m)

    # populate the status, the primal (and possibly dual) solution
    # and the primal (and possibly dual) variables with values
    populate_solution!(m, problem, var_to_ranges, conic_constraints)
    if problem.status != :Optimal && verbose
        @warn "Problem status $(problem.status); solution may be inaccurate."
    end
end

function set_warmstart!(m::MathProgBase.AbstractConicModel,
                        problem::Problem,
                        n::Int, # length of primal (conic) solution
                        var_to_ranges)
     # use previously cached solution, if any,
     try
         primal = problem.solution.primal
     catch
         @warn """
               Unable to use cached solution to warmstart problem.
               (Perhaps this is the first time you're solving this problem?)
               Warmstart may be ineffective.
               """
         primal = zeros(n)
     end
     if length(primal) != n
         @warn """
               Unable to use cached solution to warmstart problem.
               (Perhaps the number of variables or constraints in the problem have changed since you last solved it?)
               Warmstart may be ineffective.
               """
         primal = zeros(n)
     end

     # grab any variables whose values the user may be trying to set
     load_primal_solution!(primal, var_to_ranges)

     # notify the model that we're trying to warmstart
     try
         MathProgBase.setwarmstart!(m, primal)
     catch
         @warn """
               Unable to warmstart solution.
               (Perhaps the solver doesn't support warm starts?)
               Using a cold start instead.
               """
     end
     m
end

function load_problem!(m::MathProgBase.AbstractConicModel, c, A, b, cones, vartypes)
    # no conic constraints on variables
    var_cones = fill((:Free, 1:size(A, 2)),1)
    MathProgBase.loadproblem!(m, vec(Array(c)), A, vec(Array(b)), cones, var_cones)

    # add integer and binary constraints on variables
    if !all(==(:Cont), vartypes)
        try
            MathProgBase.setvartype!(m, vartypes)
        catch
            error("model $(typeof(m)) does not support variables of some of the following types: $(unique(vartypes))")
        end
    end
    m
end

function populate_solution!(m::MathProgBase.AbstractConicModel,
                            problem::Problem,
                            var_to_ranges,
                            conic_constraints)
    dual = try
        MathProgBase.getdual(m)
    catch
        fill(NaN, MathProgBase.numconstr(m))
    end

    solution = try
        MathProgBase.getsolution(m)
    catch
        fill(NaN, MathProgBase.numvar(m))
    end

    objective = try
        MathProgBase.getobjval(m)
    catch
        NaN
    end

    if any(isnan, dual)
        problem.solution = Solution(solution, MathProgBase.status(m), objective)
    else
        problem.solution = Solution(solution, dual, MathProgBase.status(m), objective)
    end

    populate_variables!(problem, var_to_ranges)

    if problem.solution.has_dual
        populate_duals!(conic_constraints, problem.solution.dual)
    end

    # minimize -> maximize
    if (problem.head == :maximize)
        problem.solution.optval = -problem.solution.optval
    end

    # Populate the problem with the solution
    problem.optval = problem.solution.optval
    problem.status = problem.solution.status

    problem
end

function populate_variables!(problem::Problem, var_to_ranges::Dict{UInt64, Tuple{Int, Int}})
    x = problem.solution.primal
    for (id, (start_index, end_index)) in var_to_ranges
        var = id_to_variables[id]
        sz = var.size
        if var.sign != ComplexSign()
            var.value = reshape(x[start_index:end_index], sz[1], sz[2])
            if sz == (1, 1)
                var.value = var.value[1]
            end
        else
            real_value = reshape(x[start_index:start_index + div(end_index-start_index+1,2)-1], sz[1], sz[2])
            imag_value = reshape(x[start_index + div(end_index-start_index+1,2):end_index], sz[1], sz[2])
            var.value = real_value + im*imag_value
            if sz == (1, 1)
                var.value = var.value[1]
            end
        end
    end
end

# populates the solution vector from the .value fields of variables
# for use in warmstarting
# TODO: it would be super cool to grab the other expressions that appear in the primal solution vector,
# get their `expression_to_range`,
# and populate them too using `evaluate`
function load_primal_solution!(primal::Array{Float64,1}, var_to_ranges::Dict{UInt64, Tuple{Int, Int}})
    for (id, (start_index, end_index)) in var_to_ranges
        var = id_to_variables[id]
        if var.value !== nothing
            sz = size(var.value)
            if length(sz) <= 1
                primal[start_index:end_index] = var.value
            else
                primal[start_index:end_index] = reshape(var.value, sz[1]*sz[2], 1)
            end
        end
    end
end

function populate_duals!(constraints::Array{ConicConstr}, dual::Vector)
    constr_index = 1
    for constraint in constraints
        # conic_constr_to_constr only has keys for conic constraints with a single objective
        # so this will work
        if haskey(conic_constr_to_constr, constraint)
            sz = constraint.sizes[1]
            c = conic_constr_to_constr[constraint]
            c.dual = reshape(dual[constr_index:constr_index+sz-1], c.size)
            if c.size == (1, 1)
                c.dual = c.dual[1]
            end
            constr_index += sz
        else
            for i = 1:length(constraint.objs)
                sz = constraint.sizes[i]
                constr_index += sz
            end
        end
    end
end
