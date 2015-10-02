export QolElemAtom, qol_elementwise, square, sumsquares, invpos
export sign, monotonicity, curvature, conic_form!

type QolElemAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::@compat Tuple{AbstractExpr, AbstractExpr}
  size::@compat Tuple{Int, Int}

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

function evaluate(q::QolElemAtom)
  return (evaluate(q.children[1]).^2) ./ evaluate(q.children[2])
end

function conic_form!(q::QolElemAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, q)
    sz = q.children[1].size
    t = Variable(sz[1], sz[2])
    qol_objective = conic_form!(t, unique_conic_forms)
    x, y = q.children
    conic_form!(SOCElemConstraint(y + t, y - t, 2 * x), unique_conic_forms)
    conic_form!(y >= 0, unique_conic_forms)
    cache_conic_form!(unique_conic_forms, q, qol_objective)
  end
  return get_conic_form(unique_conic_forms, q)
end

qol_elementwise(x::AbstractExpr, y::AbstractExpr) = QolElemAtom(x, y)
square(x::AbstractExpr) = QolElemAtom(x, Constant(ones(x.size[1], x.size[2])))
invpos(x::AbstractExpr) = QolElemAtom(Constant(ones(x.size[1], x.size[2])), x)
sumsquares(x::AbstractExpr) = square(norm2(x))
