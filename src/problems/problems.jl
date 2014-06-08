export Problem, minimize, maximize, get_var_dict, solve!, ecos_debug, getCanonicalConstraints

Float64OrNothing = Union(Float64, Nothing)
SolutionOrNothing = Union(Solution, Nothing)

# The Problem type consists of an objective and a set of a constraints.
# The objective specifies what should be maximized/minimized whereas the constraints
# specify the different constraints on the problem.
# status, optval and solution are populated after solve! is called.
type Problem
	head::Symbol
	objective::AbstractCvxExpr
	constraints::Array{CvxConstr}
	status::ASCIIString
	optval::Float64OrNothing
	solution::SolutionOrNothing

	function Problem(head::Symbol, objective::AbstractCvxExpr, constraints::Array{CvxConstr}=CvxConstr[])
		if !all([x <= 1 for x in objective.size])
			error("Only scalar optimization problems allowed, but size(objective) = $(objective.size).")
		end

		if head == :minimize && objective.vexity == :concave
			error("Cannot minimize a concave function.")
		elseif head == :maximize && objective.vexity == :convex
			error("Cannot maximize a convex function.")
		elseif head != :maximize && head != :minimize
			error("Problem.head must be one of :minimize or :maximize.")
		end

		new(head, objective, constraints, "not yet solved", nothing, nothing)
	end
end

Problem(head::Symbol, objective::AbstractCvxExpr, constraints::CvxConstr...) =
	Problem(head, objective, [constraints...])

# Allow users to simply type minimize or maximize
minimize(objective::AbstractCvxExpr, constraints::CvxConstr...) =
	Problem(:minimize, objective, [constraints...])
minimize(objective::AbstractCvxExpr, constraints::Array{CvxConstr}=CvxConstr[]) =
	Problem(:minimize, objective, constraints)
maximize(objective::AbstractCvxExpr, constraints::CvxConstr...) =
	Problem(:maximize, objective, [constraints...])
maximize(objective::AbstractCvxExpr, constraints::Array{CvxConstr}=CvxConstr[]) =
	Problem(:maximize, objective, constraints)
satisfy(constraints::Array{CvxConstr}=CvxConstr[]) =
	Problem(:minimize, Constant(0), constraints)

# +(constraints, constraints) is overwritten in constraints.jl
add_constraints(p::Problem, constraints::Array{CvxConstr}) = +(p.constraints, constraints)
add_constraints(p::Problem, constraint::CvxConstr) = add_constraints(p, [constraint])

function solve!(p::Problem, method=:ecos)
	if method == :ecos
		solution = ecos_solve(p)
	else
		println("method $method not implemented")
	end

	# Calculate the solution
	# TODO: After switching to Julia 0.3, use dot().
	optval = c' * solution.x
	problem.optval = optval[1] # Transpose returns an array, so fetch the element
	problem.status = solution.status
	problem.solution = solution

	if problem.status == "solved"
		populate_variables!(problem, variable_index)
		populate_constraints!(problem, eq_constr_index, ineq_constr_index)
	end
end

function canonical_constraints(problem::Problem)
    objective = problem.objective

    canonical_constraints_array = CanonicalConstr[]
    for constraint in problem.constraints
        append!(canonical_constraints_array, constraint.canon_form())
    end

    append!(canonical_constraints_array, objective.canon_form())
    return canonical_constraints_array
end

# Now that the problem has been solved, populate the optimal values of the
# variables back into them
function populate_variables!(problem::Problem, variable_index::Dict{Int64, Int64})
	x = problem.solution.x
	var_dict = get_var_dict(problem.objective, problem.constraints)
	for (id, var) in var_dict
		index = variable_index[id]
		var.value = Base.reshape(x[index : index + get_vectorized_size(var) - 1], var.size)
		if var.size == (1, 1)
			# Make it a scalar
			var.value = var.value[1]
		end
	end
end

function populate_constraints!(problem::Problem, eq_constr_index::Dict{Int64, Int64},
		ineq_constr_index::Dict{Int64, Int64})

	y = problem.solution.y
	z = problem.solution.z

	for constraint in problem.constraints
		uid = constraint.canon_uid
		if constraint.head == :(==)
			index = eq_constr_index[uid]
			# constraint.dual_value = Base.reshape(y[index : index + get_vectorized_size(constraint.size) - 1], constraint.size)
		else
			index = ineq_constr_index[uid]
			constraint.dual_value = Base.reshape(z[index : index + get_vectorized_size(constraint.size) - 1], constraint.size)
		end
	end
end

# Recursively traverses the AST for the AbstractCvxExpr and finds the variables
# that were defined
# Updates var_dict with the ids of the variables as keys and variables as values
function get_var_dict!(e::AbstractCvxExpr, var_dict::Dict{Int64, Variable})
	if e.head == :variable
		var_dict[e.uid] = e
	elseif e.head == :parameter || e.head == :constant
		return
	else
		for v in e.args
			get_var_dict!(v, var_dict)
		end
	end
end

# hacky way to not crash on recursion when some arguments for some atoms are
# symbols or numbers
get_var_dict!(e, var_dict) = nothing

function get_var_dict(p::Problem)
	return get_var_dict(p.objective, p.constraints)
end

function get_var_dict(objective::AbstractCvxExpr, constraints::Array{CvxConstr})
	var_dict = Dict{Int64, Variable}()

	get_var_dict!(objective, var_dict)
	for c in constraints
		get_var_dict!(c.lhs, var_dict);
		get_var_dict!(c.rhs, var_dict);
	end

	return var_dict
end
