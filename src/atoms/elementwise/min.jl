import Base.min
export min

function min(x::AbstractCvxExpr)
  # TODO: handle signs
  if x.vexity == :convex
    error("min of convex function is not DCP compliant")
  end
  # Fake vexity given so >= doesn't throw DCP compliance error
  this = CvxExpr(:min, [x], :linear, x.sign, (1, 1))

  # 'this <= x' will try to find the canon_form for 'this', so we need to initialize it
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (this <= x).canon_form()

  # Add back the correct vexity
  this.vexity = :concave

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.minimum(x.evaluate())
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

  # Fake vexity given so >= doesn't throw DCP compliance error
  sign = x.sign == y.sign ? x.sign : :any
  this = CvxExpr(:min, [x], :linear, sign, sz)
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (this <= x).canon_form()
  append!(canon_constr_array, (this <= y).canon_form())

  # Add back the correct vexity
  this.vexity = :concave

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.min(x.evaluate(), y.evaluate())
  return this
end

min(x::AbstractCvxExpr, y::Value) = min(x, convert(CvxExpr, y))
min(y::Value, x::AbstractCvxExpr) = min(x, convert(CvxExpr, y))
