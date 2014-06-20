import Base.max
export max

# TODO: Handle signs for max
# TODO: Find max across a specific dimension

# Maximum element of `x`
# Canonical constraint is x <= this if max(x) = this
function max(x::AbstractCvxExpr)
  if !is_convex(x.vexity)
    error("max of a non convex function is not DCP compliant")
  end
  # Fake vexity given so <= doesn't throw DCP compliance error
  this = CvxExpr(:max, [x], :linear, x.sign, (1, 1))

  # 'x <= this' will try to find the canon_form for 'this', so we need to initialize it
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (x <= this).canon_form()

  # Add back the correct vexity
  this.vexity = :convex

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.maximum(x.evaluate())
  return this
end

# Elementwise maximum between two expressions
# Canonical constraint is x <= this and y <= this if max(x, y) = this
function max(x::AbstractCvxExpr, y::AbstractCvxExpr)
  if !is_convex(x.vexity) || !is_convex(y.vexity)
    error("max of a non convex function is not DCP compliant")
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

  if x.sign == :pos || y.sign == :pos
    sign = :pos
  elseif x.sign == :neg && y.sign == :neg
    sign = :neg
  else
    sign = :any
  end

  # Fake vexity given so <= doesn't throw DCP compliance error
  this = CvxExpr(:max, [x, y], :linear, sign, sz)
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (x <= this).canon_form()
  append!(canon_constr_array, (y <= this).canon_form())
  # Add back the correct vexity
  this.vexity = :convex

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.max(x.evaluate(), y.evaluate())
  return this
end

max(x::AbstractCvxExpr, y::Value) = max(x, convert(CvxExpr, y))
max(y::Value, x::AbstractCvxExpr) = max(x, convert(CvxExpr, y))
