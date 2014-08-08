import Base.show
export show

# A Constant is simply a wrapper around a native Julia constant
# Hence, we simply display its value
function show(io::IO, x::Constant)
  print(io, "$(x.value)")
end

# A variable, for example, Variable(3, 4), will be displayed as:
# Variable((3, 4), any)
function show(io::IO, x::Variable)
  print(io, "Variable(($(x.size[1]), $(x.size[2])), $(x.sign))")
end

# A constraint, for example, square(x) <= 4 will be displayed as:
# CvxConstr: CvxExpr((1, 1), convex) <= 4
function show(io::IO, c::CvxConstr)
  print(io, "CvxConstr: $(c.lhs) $(c.head) $(c.rhs)")
end

# An expression, for example, 2 * x, will be displayed as:
# CvxExpr((1, 1), affine, any)
function show(io::IO, e::CvxExpr)
  print(io, "CvxExpr(($(e.size[1]), $(e.size[2])), $(e.vexity), $(e.sign))")
end

# A problem, for example, p = minimize(sum(x) + 3, [x >= 0, x >= 1, x <= 10])
# will be displayed as follows:
#
# Problem: minimize CvxExpr((1, 1), affine, any)
#    subject to
#      CvxConstr: CvxExpr((2, 3), affine, any) <= 0
#      CvxConstr: CvxExpr((2, 3), affine, any) <= -1
#      CvxConstr: Variable((2, 3), any) <= 10
#    current status: not yet solved
#
# Once it is solved, the current status would look like:
#    current status: solved with optimal value of 9.0
function show(io::IO, p::Problem)
  print(io, "Problem: $(p.head) $(p.objective) \n\t subject to\n\t\t")
  print_joined(io, p.constraints, "\n\t\t")
  print(io, "\n\t current status: $(p.status)")
  if p.status == "solved"
    print(io, " with optimal value of $(round(p.optval, 4))")
  end
end

# The remaining functions are for internal use and debugging purposes since
# these types are not exposed to the end-user
# Hence, these are more detailed than/not consistent with the display functions
# above

function show(io::IO, solution::Solution)
  print(io, "Solution:\n\tvalue returned by solver: $(solution.optval)\n")
  print(io, "\tstatus: $(solution.status)\n")
  print(io, "\toptimal value of primal variables: $(solution.x)\n")
  print(io, "\toptimal value of dual equality variables: $(solution.y)\n")
  print(io, "\toptimal value of dual inequality variables: $(solution.z)")
end

function show(io::IO, canon_constr::CanonicalConstr)
  print(io, "CanonicalConstr: $(canon_constr.uid)\n")
  print(io, "\tuid of expressions: $(canon_constr.vars)\n")
  print(io, "\tcoefficients of above expressions: $(canon_constr.coeffs)\n")
  print(io, "\tconstant: $(canon_constr.constant)\n")
  print(io, "\tequality constraint: $(canon_constr.is_eq)\n")
  print(io, "\tconic constraint: $(canon_constr.is_conic)")
end
