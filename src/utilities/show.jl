import Base.show
export show

# A Constant is simply a wrapper around a native Julia constant
# Hence, we simply display its value
function show(io::IO, x::Constant)
  print(io, "$(x.value)")
end

# A variable, for example, Variable(3, 4), will be displayed as:
# Variable of
# size: (3, 4)
# sign: NoSign()
# vexity: AffineVexity()
function show(io::IO, x::Variable)
  print(io, """Variable of
    size: ($(x.size[1]), $(x.size[2]))
    sign: $(x.sign)
    vexity: $(x.vexity)""")
  if x.value != nothing
    print(io, "\nvalue: $(x.value)")
  end
end

# A constraint, for example, square(x) <= 4 will be displayed as:
# Constraint:
# <= constraint
# lhs: ...
# rhs: ...
function show(io::IO, c::Constraint)
  print(io, """Constraint:
    $(c.head) constraint
    lhs: $(c.lhs)
    rhs: $(c.rhs)
    """)
end

# SDP constraints are displayed as:
# Constraint:
# sdp constraint
# expression: ...
function show(io::IO, c::SDPConstraint)
  print(io, """Constraint:
    $(c.head) constraint
    expression: $(c.lhs)
    """)
end

# An expression, for example, 2 * x, will be displayed as:
# AbstractExpr with
# head: *
# size: (1, 1)
# sign: NoSign()
# vexity: AffineVexity()
function show(io::IO, e::AbstractExpr)
  print(io, """AbstractExpr with
    head: $(e.head)
    size: ($(e.size[1]), $(e.size[2]))
    sign: $(sign(e))
    vexity: $(vexity(e))
    """)
end

# A problem, for example, p = minimize(sum(x) + 3, [x >= 0, x >= 1, x <= 10])
# will be displayed as follows:
#
# Problem: minimize `Expression`
#    subject to
#      Constraint: ...
#      Constraint: ...
#      Constraint: ...
#    current status: not yet solved
#
# Once it is solved, the current status would look like:
#    current status: solved with optimal value of 9.0
function show(io::IO, p::Problem)
  print(io, """Problem:
    $(p.head) $(p.objective)
    subject to
    """)
  print_joined(io, p.constraints, "\n\t\t")
  print(io, "\ncurrent status: $(p.status)")
  if p.status == "solved"
    print(io, " with optimal value of $(round(p.optval, 4))")
  end
end
