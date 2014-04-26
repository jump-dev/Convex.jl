export Problem, minimize, maximize, get_variables, solve!

type Problem
	head::Symbol
	obj::AbstractCvxExpr
	constr::Array{CvxConstr}
	status
	optval
	variables

	function Problem(head::Symbol,obj::AbstractCvxExpr,constr=CvxConstr[]::Array{CvxConstr})
		if !all([x<=1 for x in obj.size])
			error("only scalar optimization problems allowed. got size(obj) = $(obj.size)")
		end

		if head == :minimize && obj.vexity == :concave
				error("cannot minimize a concave function")
		elseif head == :maximize && obj.vexity == :convex
				error("cannot maximize a convex function")
		elseif head != :maximize && head != :minimize
			error("Problem.head must be one of :minimize, :maximize")
		end

		new(head, obj, constr, nothing, nothing, get_var_dict(obj, constr))
	end
end

# infer min or max from vexity of obj
# Problem(obj::AbstractCvxExpr,constr::Array{CvxConstr}) = obj.vexity == :concave ? Problem(:maximize,obj,constr) : Problem(:minimize,obj,constr)
Problem(head::Symbol, obj::AbstractCvxExpr, constr::CvxConstr...) = Problem(head, obj, [constr...])
minimize(obj::AbstractCvxExpr, constr::CvxConstr...) = Problem(:minimize, obj, [constr...])
minimize(obj::AbstractCvxExpr, constr=CvxConstr[]::Array{CvxConstr}) = Problem(:minimize, obj, constr)
maximize(obj::AbstractCvxExpr, constr::CvxConstr...) = Problem(:maximize, obj, [constr...])
maximize(obj::AbstractCvxExpr, constr=CvxConstr[]::Array{CvxConstr}) = Problem(:maximize, obj, constr)
# display(p::Problem) = idea: loop through variables, name them sequentially if they don't already have names, and display operations prettily with variables names
# println("Problem($(p.head),$(p.obj),$(p.constr)")


function solve!(p::Problem,method=:ecos)
	if method == :ecos
		ecos_solve!(p)
	else
		println("method $method not implemented")
	end
end


# For now, assuming it is an LP, so objective is of the form c' * x
function ecos_solve!(problem::Problem)
	objective = problem.obj

	canonical_constraints_array = Any[]
	for constraint in problem.constr
		append!(canonical_constraints_array, constraint.canon_form())
	end

	append!(canonical_constraints_array, objective.canon_form())
	m, n, p, G, h, A, b, variable_index = create_ecos_matrices(canonical_constraints_array)
	# Now, all we need to is create c

	c = zeros(n, 1)

	# TODO: Handle cases such as minimize 1 or the case where objective isn't in variable index
	uid = objective.uid()
	c[variable_index[uid] : variable_index[uid] + objective.size[1] - 1] = 1

	if problem.head == :maximize
		c = -c;
	end

	sol = ecos_solve(n=n, m=m, p=p, G=G, c=c, h=h, A=A, b=b)
	if problem.head == :maximize
		c = -c;
	end
	problem.optval = c' * sol[:x]
	# TODO: After switching to Julia 0.3, use dot().
	# Transpose returns an array, so fetch the element
	problem.optval = problem.optval[1]
	problem.status = sol[:status]
end


# Given the canonical_constraints_array, creates conic inequality matrix G and h
# Also creates equality matrix A and b
function create_ecos_matrices(canonical_constraints_array)
	n = 0
	variable_index = Dict()
	m = 0
	p = 0

	# Loop over all the constraints to figure out the size of G and A
	for constraint in canonical_constraints_array
		# Loop over each variable in the constraint
		for i = 1:length(constraint[:vars])
			var = constraint[:vars][i]

			# If we haven't already taken into account the size of this variable,
			# add it to the size of the variable
			if !haskey(variable_index, var)
				variable_index[var] = n + 1

				n += size(constraint[:coeffs][i], 2)
			end

			if constraint[:is_eq]
				p += size(constraint[:coeffs][i], 1)
			else
				m += size(constraint[:coeffs][i], 1)
			end
		end
	end

	h = m == 0 ? nothing: zeros(m, 1)
	G = m == 0 ? nothing: spzeros(m, n)
	b = p == 0 ? nothing: zeros(p, 1)
	A = p == 0 ? nothing: spzeros(p, n)

	m_index = 1::Int64
	p_index = 1::Int64

	# Now, we actually stuff the matrices A and G
	for constraint in canonical_constraints_array
		m_var = 0::Int64
		for i = 1:length(constraint[:vars])
			var = constraint[:vars][i]
			# Technically, the m_var size of all the variables should be the same, otherwise nothing makes
			# sense
			m_var = size(constraint[:coeffs][i], 1)
			n_var = size(constraint[:coeffs][i], 2)

			if constraint[:is_eq]
				# TODO: Julia has problems not converting ints to floats
				# An issue has been filed and fixed in newer versions of julia
				A[p_index : p_index + m_var - 1, variable_index[var] : variable_index[var] + n_var - 1] =
					constraint[:coeffs][i] * 1.0
			else

				G[m_index : m_index + m_var - 1, variable_index[var] : variable_index[var] + n_var - 1] =
					constraint[:coeffs][i] * 1.0
			end
		end

		if constraint[:is_eq]
			b[p_index : p_index + m_var - 1] = constraint[:constant]
			p_index += m_var
		else
			h[m_index : m_index + m_var - 1] = constraint[:constant]
			m_index += m_var
		end
	end

	return m, n, p, G, h, A, b, variable_index
end


function get_variables(e::AbstractCvxExpr)
	return nothing
	# TODO: this is also broken
	if e.head == :variable
		return [(unique_id(e),e)]
	elseif e.head == :parameter
		return [(unique_id(e),e)]
	elseif e.head == :constant
		return Set()
	else
		vars = Set()
		union!(vars, [(unique_id(e),e)])

		for v in e.args
			subvars = get_variables(v)
			if subvars != nothing
				union!(vars,subvars)
			end
		end
		return vars
	end
end
# hacky way to not crash on recursion when some arguments for some atoms are symbols or numbers, as in the case of norm
get_variables(e) = Set()

function get_var_dict(p::Problem)
	# TODO: this is also broken
	return nothing
	get_var_dict(p.obj,p.constr)
end

function get_var_dict(obj::AbstractCvxExpr, constr::Array{CvxConstr})
	# TODO: this is also broken
	return nothing
	vars = get_variables(obj)
	for c in constr
		subvars = get_variables(c.lhs);
		subvars!=nothing ? union!(vars,subvars) : nothing

		subvars = get_variables(c.rhs);
		subvars!=nothing ? union!(vars,subvars) : nothing
	end

	var_dict = Dict()
	for (p,v) in vars
		var_dict[p] = v
	end

	var_dict
end
