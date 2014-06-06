export Problem, minimize, maximize, get_var_dict, solve!, ecos_debug

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
			error("Only scalar optimization problems allowed, but size(objective) = $(objective.size)")
		end

		if head == :minimize && objective.vexity == :concave
			error("Cannot minimize a concave function")
		elseif head == :maximize && objective.vexity == :convex
			error("Cannot maximize a convex function")
		elseif head != :maximize && head != :minimize
			error("Problem.head must be one of :minimize or :maximize")
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
		ecos_solve!(p)
	else
		println("Method $method not implemented")
	end
end

function ecos_debug(problem::Problem)
  objective = problem.objective

	canonical_constraints_array = CanonicalConstr[]
	for constraint in problem.constraints
		append!(canonical_constraints_array, constraint.canon_form())
	end

	append!(canonical_constraints_array, objective.canon_form())
	return create_ecos_matrices(canonical_constraints_array)
end

# CAUTION: For now, we assume we are solving a linear program.
# Loops over the objective and constraints to get the canonical constraints array.
# It then calls create_ecos_matrices which will create the inequality/equality matrices
# and their corresponding coefficients, which are then passed to ecos_solve
function ecos_solve!(problem::Problem)
	objective = problem.objective

	canonical_constraints_array = CanonicalConstr[]
	for constraint in problem.constraints
		append!(canonical_constraints_array, constraint.canon_form())
	end

	append!(canonical_constraints_array, objective.canon_form())
	m, n, p, l, ncones, q, G, h, A, b, variable_index, eq_constr_index, ineq_constr_index =
		create_ecos_matrices(canonical_constraints_array)
	#return create_ecos_matrices(canonical_constraints_array)


	# Now, all we need to is create c
	c = zeros(n, 1)

	if objective.vexity != :constant
		uid = objective.uid
		c[variable_index[uid] : variable_index[uid] + objective.size[1] - 1] = 1
	end

	if problem.head == :maximize
		c = -c;
	end

	solution = ecos_solve(n=n, m=m, p=p, l=l, ncones=ncones, q=q, G=G, c=c, h=h, A=A, b=b)

	# Change c back to what it originally was
	if problem.head == :maximize
		c = -c;
	end

	# Calculate the optimum solution
	# TODO: After switching to Julia 0.3, use dot().
	optval = c' * solution.x

	# Transpose returns an array, so fetch the element
	problem.optval = optval[1]
	problem.status = solution.status

	problem.solution = solution
	if problem.status == "solved"
		populate_variables!(problem, variable_index)
		populate_constraints!(problem, eq_constr_index, ineq_constr_index)
	end
end


# Given the canonical_constraints_array, creates conic inequality matrix G and h
# as well as the equality matrix A and b
function create_ecos_matrices(canonical_constraints_array)
	n = 0::Int64
	variable_index = Dict{Int64, Int64}()

	eq_constr_index = Dict{Int64, Int64}()
	ineq_constr_index = Dict{Int64, Int64}()

	m = 0::Int64
	p = 0::Int64
	l = 0::Int64
	ncones = 0::Int64
	q = Int64[]

	# Loop over all the constraints to figure out the size of G and A
	for constraint in canonical_constraints_array
		# Loop over each variable in the constraint
		length_constraint_vars = length(constraint.vars)
		for i = 1:length_constraint_vars
			var = constraint.vars[i]

			# If we haven't already taken into account the size of this variable,
			# add it to the size of the variable
			if !haskey(variable_index, var)
				variable_index[var] = n + 1

				n += size(constraint.coeffs[i], 2)
			end
		end

		if constraint.is_eq
			p += size(constraint.coeffs[1], 1)
		elseif constraint.is_conic
			ncones += 1
			push!(q, size(constraint.coeffs[1], 1))
			m += size(constraint.coeffs[1], 1)
		else
			l += size(constraint.coeffs[1], 1)
			m += size(constraint.coeffs[1], 1)
		end
	end

	h = m == 0 ? nothing: zeros(m, 1)
	G = m == 0 ? nothing: spzeros(m, n)
	b = p == 0 ? nothing: zeros(p, 1)
	A = p == 0 ? nothing: spzeros(p, n)

	l_index = 1::Int64
	c_index = 1::Int64
	p_index = 1::Int64

	# Now, we actually stuff the matrices A and G
	for constraint in canonical_constraints_array
		m_var = 0::Int64

		length_constraint_vars = length(constraint.vars)
		for i = 1:length_constraint_vars
			var = constraint.vars[i]
			# Technically, the m_var size of all the variables should be the same,
			# otherwise nothing makes sense
			m_var = size(constraint.coeffs[i], 1)
			n_var = size(constraint.coeffs[i], 2)

			if constraint.is_eq
				# TODO: Julia has problems not converting ints to floats
				# An issue has been filed and should be fixed in newer versions of julia
				A[p_index : p_index + m_var - 1, variable_index[var] : variable_index[var] + n_var - 1] =
					constraint.coeffs[i] * 1.0
			elseif constraint.is_conic
        G[l + c_index : l + c_index + m_var - 1, variable_index[var] : variable_index[var] + n_var - 1] =
					constraint.coeffs[i] * 1.0
			else
				G[l_index : l_index + m_var - 1, variable_index[var] : variable_index[var] + n_var - 1] =
					constraint.coeffs[i] * 1.0
			end
		end

		if constraint.is_eq
			eq_constr_index[constraint.uid] = p_index
			b[p_index : p_index + m_var - 1] = constraint.constant
			p_index += m_var
		elseif constraint.is_conic
			ineq_constr_index[constraint.uid] = l + c_index
      h[l + c_index : l + c_index + m_var - 1] = constraint.constant
      c_index += m_var
		else
			ineq_constr_index[constraint.uid] = l_index
			h[l_index : l_index + m_var - 1] = constraint.constant
			l_index += m_var
		end
	end

	return m, n, p, l, ncones, q, G, h, A, b, variable_index, eq_constr_index, ineq_constr_index
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
			constraint.dual_value = Base.reshape(y[index : index + get_vectorized_size(constraint.size) - 1], constraint.size)
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
