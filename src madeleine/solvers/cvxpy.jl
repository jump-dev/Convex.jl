using PyCall
@pyimport cvxpy
@pyimport numpy
@pyimport cvxopt

# export cvxpy_solve!

pyconstrmap = {:(==) => :__eq__, :(>=) => :__ge__, :(<=) => :__le__}

function cvxpy_solve!(p::Problem)
	pyvars = [id=>make_pyvar(p.variables[id]) for id in keys(p.variables)]
	pyobj = p.head == :minimize ? cvxpy.Minimize(pyify(p.obj, pyvars)) : cvxpy.Maximize(pyify(p.obj, pyvars))
	pyconstr = [pyify(c.lhs, pyvars)[pyconstrmap[c.head]](pyify(c.rhs, pyvars)) for c in p.constr]
	pyprob = cvxpy.Problem(pyobj,pyconstr)
	p.optval = pyprob[:solve]()
	# this should work but doesn't
	# p.status = pyprob[:status]
	# assign values to variables
	for id in keys(p.variables)
		p.variables[id].value = numpy.matrix(pyvars[id][:value])
	end
	# assign dual_values to constraints
	return p.optval
end

function pyify(e::AbstractCvxExpr, pyvars)
	# println("pyifying $e")
	if e.head == :variable
		return pyvars[unique_id(e)]
	elseif e.head == :parameter
		return pyvars[unique_id(e)]
	elseif e.head == :constant
		return float64(e.value) # PyCall doesn't transform integers correctly
	elseif e.head == :-
		return (pyify(e.args[1], pyvars))[:__neg__]()
	elseif e.head == :+
		return pyify(e.args[1], pyvars)[:__add__](pyify(e.args[2], pyvars))
	elseif e.head == :*
		a = pyify(e.args[1], pyvars)
		b = pyify(e.args[2], pyvars)
		try
			out = a[:__mul__](b)
			if out != pyeval("NotImplemented")
				return out
			else
				error("Multiplication between $a and $b is not implemented")
			end
		catch
			try
				out = b[:__rmul__](a)
				if out != pyeval("NotImplemented")
					return out
				else
					error("Multiplication between $a and $b is not implemented")
				end
			catch
				error("Multiplication between $a and $b is not implemented")
			end
		end
	else # an atom!
		return cvxpy.atoms[e.head]([pyify(arg, pyvars) for arg in e.args]...)
	end
end
# hacky way to not crash on recursion when some arguments for some atoms are symbols or numbers, as in the case of norm
pyify(e, pyvars) = e

function make_2d(t::Tuple)
	if length(t) == 0
		return (1,1)
	elseif length(t) == 1
		return (t[1],1)
	else
		return t[1:2]
	end
end

function make_pyvar(e::AbstractCvxExpr)
	if e.head == :variable
		return cvxpy.Variable(make_2d(e.size)...)
	elseif e.head == :parameter
		return cvxpy.Parameter(make_2d(e.size)...)
	end
end
