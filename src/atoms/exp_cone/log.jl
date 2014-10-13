#############################################################################
# log.jl
# natural logarithm of an logression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read logressions.jl first.
#############################################################################

import Base.log
export log
export sign, curvature, monotonicity, evaluate

### Logonential

type LogAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function LogAtom(x::AbstractExpr)
    if (x.size != (1, 1))
      error("TODO: Only scalar variables supported for log as of now")
    end
    children = (x,)
    return new(:log, hash(children), children, x.size)
  end
end

function sign(x::LogAtom)
  return NoSign()
end

function monotonicity(x::LogAtom)
  return (Nondecreasing(),)
end

function curvature(x::LogAtom)
  return ConcaveVexity()
end

function evaluate(x::LogAtom)
  return log(evaluate(x.children[1]))
end

log(x::AbstractExpr) = LogAtom(x)

function conic_form(e::LogAtom, unique_constr)
  if !((e.head, e.children_hash) in keys(unique_constr))
    # log(z) \geq x  <=>  (x,ones(),z) \in ExpCone
    z = e.children[1]
    y = Constant(ones(size(z)))
    x = Variable(size(z))

    constraints = ConicConstr[]
    objective_x, constraints_x = conic_form(x, unique_constr)
    append!(constraints, constraints_x)
    objective_y, constraints_y = conic_form(y, unique_constr)
    append!(constraints, constraints_y)
    objective_z, constraints_z = conic_form(z, unique_constr)
    append!(constraints, constraints_z)
    exp_constraint = ConicConstr([objective_x, objective_y, objective_z], :ExpPrimal, [1, 1, 1])
    push!(constraints, exp_constraint)

    unique_constr[(e.head, e.children_hash)] = (objective_x, constraints)
  end
  return safe_copy(unique_constr[(e.head, e.children_hash)])
end
