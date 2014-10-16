export QuadOverLinAtom, quad_over_lin
export sign, monotonicity, curvature, conic_form

type QuadOverLinAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
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
  return (sign(q.children[1]) * Nondecreasing(), Nonincreasing())
end

function curvature(q::QuadOverLinAtom)
  return ConvexVexity()
end

function conic_form(q::QuadOverLinAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, q)
    t = Variable()
    qol_objective = conic_form(t, unique_conic_forms)
    x, y = q.children
    conic_form(SOCConstraint(y + t, y - t, 2 * x), unique_conic_forms)
    conic_form(y >= 0, unique_conic_forms)
    cache_conic_form!(unique_conic_forms, q, qol_objective)
  end
  return get_conic_form(unique_conic_forms, q)
end

quad_over_lin(x::AbstractExpr, y::AbstractExpr) = QuadOverLinAtom(x, y)
