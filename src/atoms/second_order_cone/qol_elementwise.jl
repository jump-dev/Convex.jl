export QolElemAtom, qol_elementwise, square, sum_squares, inv_pos
export sign, monotonicity, curvature, conic_form

type QolElemAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr, AbstractExpr)
  size::(Int64, Int64)

  function QolElemAtom(x::AbstractExpr, y::AbstractExpr)
    if x.size != y.size
      error("elementwise quad over lin must take two arguments of the same size")
    end
    children = (x, y)
    return new(:qol_elem, hash(children), children, x.size)
  end
end

function sign(q::QolElemAtom)
  return Positive()
end

function monotonicity(q::QolElemAtom)
  return (sign(q.children[1]) * Nondecreasing(), Nonincreasing())
end

function curvature(q::QolElemAtom)
  return ConvexVexity()
end

function conic_form(q::QolElemAtom, unique_constr)
  if !((q.head, q.children_hash) in keys(unique_constr))
    sz = q.children[1].size
    t = Variable(sz[1], sz[2])
    qol_objective, qol_constraints = conic_form(t, unique_constr)
    x, y = q.children
    for row in 1:sz[1]
      for col in 1:sz[2]
        y_plus_t, y_plus_t_constr  = conic_form(y[row, col] + t[row,col], unique_constr)
        y_minus_t, y_minus_t_constr = conic_form(y[row, col] - t[row,col], unique_constr)
        x_obj, x_constr = conic_form(2 * x[row, col], unique_constr)
        soc_constraint = ConicConstr([y_plus_t, y_minus_t, x_obj], :SOC, [1, 1, 1])
        append!(qol_constraints, y_plus_t_constr)
        append!(qol_constraints, y_minus_t_constr)
        append!(qol_constraints, x_constr)
        push!(qol_constraints, soc_constraint)
      end
    end
    y_pos, y_pos_constr = conic_form(y >= 0, unique_constr)
    append!(qol_constraints, y_pos_constr)
    unique_constr[(q.head, q.children_hash)] = (qol_objective, qol_constraints)
  end
  return safe_copy(unique_constr[(q.head, q.children_hash)])
end

qol_elementwise(x::AbstractExpr, y::AbstractExpr) = QolElemAtom(x, y)
square(x::AbstractExpr) = QolElemAtom(x, Constant(ones(x.size[1], x.size[2])))
inv_pos(x::AbstractExpr) = QolElemAtom(Constant(ones(x.size[1], x.size[2])), x)
sum_squares(x::AbstractExpr) = sum(square(x))
