export GeoMeanAtom, geo_mean, sqrt
export sign, monotonicity, curvature, conic_form

type GeoMeanAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr, AbstractExpr)
  size::(Int64, Int64)

  function GeoMeanAtom(x::AbstractExpr, y::AbstractExpr)
    if x.size != y.size
      error("geo mean must take two arguments of the same size")
    end
    children = (x, y)
    return new(:geo_mean, hash(children), children, x.size)
  end
end

function sign(q::GeoMeanAtom)
  return Positive()
end

function monotonicity(q::GeoMeanAtom)
  return (Nondecreasing(), Nondecreasing())
end

function curvature(q::GeoMeanAtom)
  return ConcaveVexity()
end

function conic_form(q::GeoMeanAtom, unique_constr)
  if !((q.head, q.children_hash) in keys(unique_constr))
    sz = q.children[1].size
    t = Variable(sz[1], sz[2])
    qol_objective, qol_constraints = conic_form(t, unique_constr)
    x, y = q.children
    for row in 1:sz[1]
      for col in 1:sz[2]
        y_plus_x, y_plus_x_constr  = conic_form(y[row, col] + x[row,col], unique_constr)
        y_minus_x, y_minus_x_constr = conic_form(y[row, col] - x[row,col], unique_constr)
        t_obj, t_constr = conic_form(2 * t[row, col], unique_constr)
        soc_constraint = ConicConstr([y_plus_x, y_minus_x, t_obj], :SOC, [1, 1, 1])
        append!(qol_constraints, y_plus_x_constr)
        append!(qol_constraints, y_minus_x_constr)
        append!(qol_constraints, t_constr)
        push!(qol_constraints, soc_constraint)
      end
    end
    x_pos, x_pos_constr = conic_form(x >= 0, unique_constr)
    y_pos, y_pos_constr = conic_form(y >= 0, unique_constr)
    append!(qol_constraints, x_pos_constr)
    append!(qol_constraints, y_pos_constr)
    unique_constr[(q.head, q.children_hash)] = (qol_objective, qol_constraints)
  end
  return safe_copy(unique_constr[(q.head, q.children_hash)])
end

geo_mean(x::AbstractExpr, y::AbstractExpr) = GeoMeanAtom(x, y)
sqrt(x::AbstractExpr) = GeoMeanAtom(x, ones(x.size[1], x.size[2]))
