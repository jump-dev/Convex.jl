export min, max, pos, neg

function max(x::AbstractCvxExpr)
  # TODO: handle signs
  if x.vexity == :concave
    error("max of concave function is not DCP compliant")
  end
  this = CvxExpr(:max, [x], :convex, :pos, (1, 1))

  # 'x <= this' will try to find the canon_form for 'this', so we need to initialize it
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (x <= this).canon_form()
  this.canon_form = ()->canon_constr_array
  return this
end

function max(x::AbstractCvxExpr, y::AbstractCvxExpr)
  # TODO: handle signs
  if x.vexity == :concave || y.vexity == :concave
    error("max of concave function is not DCP compliant")
  end

  if x.size == y.size
    sz = x.size
  elseif x.size == (1, 1)
    sz = y.size
  elseif y.size == (1, 1)
    sz = x.size
  else
    error("Got different sizes for x as $(x.size) and y as $(y.size)")
  end

  this = CvxExpr(:max, [x, y], :convex, :pos, sz)
  println(x.vexity)
  println(y.vexity)
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (x <= this).canon_form()
  append!(canon_constr_array, (y <= this).canon_form())

  this.canon_form = ()->canon_constr_array
  return this
end

max(x::AbstractCvxExpr, y::Value) = max(x, convert(CvxExpr, y))
max(y::Value, x::AbstractCvxExpr) = max(x, convert(CvxExpr, y))


function min(x::AbstractCvxExpr)
  # TODO: handle signs
  if x.vexity == :convex
    error("min of convex function is not DCP compliant")
  end
  this = CvxExpr(:min, [x], :concave, :neg, (1, 1))

  # 'x <= this' will try to find the canon_form for 'this', so we need to initialize it
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (this <= x).canon_form()
  this.canon_form = ()->canon_constr_array
  return this
end

function min(x::AbstractCvxExpr, y::AbstractCvxExpr)
  # TODO: handle signs
  if x.vexity == :convex || y.vexity == :convex
    error("min of convex function is not DCP compliant")
  end

  if x.size == y.size
    sz = x.size
  elseif x.size == (1, 1)
    sz = y.size
  elseif y.size == (1, 1)
    sz = x.size
  else
    error("Got different sizes for x as $(x.size) and y as $(y.size)")
  end

  this = CvxExpr(:min, [x], x.vexity, :pos, sz)
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (this <= x).canon_form()
  append!(canon_constr_array, (this <= y).canon_form())

  this.canon_form = ()->canon_constr_array
  return this
end

min(x::AbstractCvxExpr, y::Value) = min(x, convert(CvxExpr, y))
min(y::Value, x::AbstractCvxExpr) = min(x, convert(CvxExpr, y))

function pos(x::AbstractCvxExpr)
  this = max(x, spzeros(x.size...))
  this.sign = :pos
  return this
end

function neg(x::AbstractCvxExpr)
  this = max(-x, 0)
  this.sign = :neg
  return this
end
