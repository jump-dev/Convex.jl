export Problem, minimize, maximize, get_variables, solve!

type Problem
	head::Symbol
	obj::CvxExpr
	constr::Array{CvxConstr}
	status
	optval
	variables
	function Problem(head::Symbol,obj::CvxExpr,constr=CvxConstr[]::Array{CvxConstr})
		if !all([x<=1 for x in obj.size])
			error("only scalar optimization problems allowed. got size(obj) = $(obj.size)")
		end
		if head == :minimize
			if obj.vexity == :concave
				error("cannot minimize a concave function")
			end
		elseif head == :maximize
			if obj.vexity == :convex
				error("cannot maximize a convex function")
			end
		else
			error("Problem.head must be one of :minimize,:maximize")
		end
		new(head,obj,constr,nothing,nothing,get_var_dict(obj,constr))
	end
end
# infer min or max from vexity of obj
# Problem(obj::CvxExpr,constr::Array{CvxConstr}) = obj.vexity == :concave ? Problem(:maximize,obj,constr) : Problem(:minimize,obj,constr)
Problem(head::Symbol,obj::CvxExpr,constr::CvxConstr...) = Problem(head,obj,[constr...])
minimize(obj::CvxExpr,constr::CvxConstr...) = Problem(:minimize,obj,[constr...])
minimize(obj::CvxExpr,constr=CvxConstr[]::Array{CvxConstr}) = Problem(:minimize,obj,constr)
maximize(obj::CvxExpr,constr::CvxConstr...) = Problem(:maximize,obj,[constr...])
maximize(obj::CvxExpr,constr=CvxConstr[]::Array{CvxConstr}) = Problem(:maximize,obj,constr)
# display(p::Problem) = idea: loop through variables, name them sequentially if they don't already have names, and display operations prettily with variables names
# println("Problem($(p.head),$(p.obj),$(p.constr)")

function solve!(p::Problem,method=:cvxpy)
	if method == :cvxpy
		return cvxpy_solve!(p)
	else
		println("method $method not implemented")
	end
end

function get_variables(e::AbstractCvxExpr)
	if e.head == :variable
		return [(unique_id(e),e)]
	elseif e.head == :parameter
		return [(unique_id(e),e)]
	elseif e.head == :constant
		return Set()
	else
		vars = Set()
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
	get_var_dict(p.obj,p.constr)
end

function get_var_dict(obj::CvxExpr,constr::Array{CvxConstr})
	vars = get_variables(obj)
	for c in constr
		subvars = get_variables(c.lhs); subvars!=nothing ? union!(vars,subvars) : nothing
		subvars = get_variables(c.rhs); subvars!=nothing ? union!(vars,subvars) : nothing
	end
	var_dict = Dict()
	for (p,v) in vars
		var_dict[p] = v
	end
	var_dict
end
