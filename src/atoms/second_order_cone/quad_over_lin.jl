export QuadOverLinAtom, quad_over_lin
export sign, monotonicity, curvature, conic_form

type QuadOverLinAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr, AbstractExpr)
  size::(Int64, Int64)

  function QuadOverLinAtom(x::AbstractExpr, y::AbstractExpr)
    if x.size[2] != 1 && y.size != (1, 1)
      error("quad over lin arguments must be a vector and a scalar")
    end
    children = (x,y)
    return new(:qol, hash(children), children, (1, 1))
  end
end

function sign(q::QuadOverLinAtom)
  return Positive()
end

function monotonicity(q::QuadOverLinAtom)
  return (sign(x.children[1]) * Nondecreasing(), Nonincreasing())
end

function curvature(q::QuadOverLinAtom)
  return ConvexVexity()
end

function conic_form(q::QuadOverLinAtom)
  if !((q.head, q.children_hash) in keys(unique_constr))
    t = Variable()
    qol_objective, qol_constraints = conic_form(t, unique_constr)
    y_plus_t, y_plus_t_constr  = conic_form(q.children[2] + t, unique_constr)
    y_minus_t, y_minus_t_constr = conic_form(q.children[2] - t, unique_constr)
    x_obj, x_constr = conic_form(2*q.children[1], unique_constr)
    y_pos, y_pos_cosntr = conic_form(q.children[2] >= 0, unique_constr)
    soc_constraint = ConicConstr([y_plus_t, y_minus_t, x_obj], :SOC, [1, 1, get_vectorized_size(q.children[1])])
    append!(qol_constraints, y_plus_t_constr)
    append!(qol_constraints, y_minus_t_constr)
    append!(qol_constraints, x_constr)
    append!(qol_constraints, y_pos_constr)
    push!(qol_constraints, soc_constraint)
    unique_constr[(q.head, q.children_hash)] = (qol_objective, qol_constraints)
  end
  return safe_copy(unique_constr[(q.head, q.children_hash)])
end

quad_over_lin(x::AbstractExpr, y::AbstractExpr) = QuadOverLinAtom(x, y)
