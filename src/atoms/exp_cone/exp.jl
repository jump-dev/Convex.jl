#############################################################################
# exp.jl
# e raised to the power of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.exp
export exp
export sign, curvature, monotonicity, evaluate

### Exponential

type ExpAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function ExpAtom(x::AbstractExpr)
    if (x.size != (1, 1))
      error("TODO: Only scalar variables supported for exp as of now")
    end
    children = (x,)
    return new(:exp, hash(children), children, x.size)
  end
end

function sign(x::ExpAtom)
  return Positive()
end

function monotonicity(x::ExpAtom)
  return (Nondecreasing(),)
end

function curvature(x::ExpAtom)
  return ConvexVexity()
end

function evaluate(x::ExpAtom)
  return exp(evaluate(x.children[1]))
end

exp(x::AbstractExpr) = ExpAtom(x)

function conic_form(e::ExpAtom, unique_constr)
  if !((e.head, e.children_hash) in keys(unique_constr))
    # exp(x) \leq z  <=>  (x,ones(),z) \in ExpCone
    x = e.children[1]
    y = Constant(ones(size(x)))
    z = Variable(size(x))

    constraints = ConicConstr[]
    objective_x, constraints_x = conic_form(x, unique_constr)
    append!(constraints, constraints_x)
    objective_y, constraints_y = conic_form(y, unique_constr)
    append!(constraints, constraints_y)
    objective_z, constraints_z = conic_form(z, unique_constr)
    append!(constraints, constraints_z)
    exp_constraint = ConicConstr([objective_x, objective_y, objective_z], :ExpPrimal, [1, 1, 1])
    push!(constraints, exp_constraint)

    unique_constr[(e.head, e.children_hash)] = (objective_z, constraints)
  end
  return safe_copy(unique_constr[(e.head, e.children_hash)])
end
