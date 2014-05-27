export max

function max(x::AbstractCvxExpr)
  # TODO: handle signs
  if x.vexity == :concave
    error("max of concave function is not DCP compliant")
  end
  # Fake vexity given so <= doesn't throw DCP compliance error
  this = CvxExpr(:max, [x], :linear, :pos, (1, 1))

  # 'x <= this' will try to find the canon_form for 'this', so we need to initialize it
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (x <= this).canon_form()

  # Add back the correct vexity
  this.vexity = :convex

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.maximum(x.evaluate())
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

  # Fake vexity given so <= doesn't throw DCP compliance error
  this = CvxExpr(:max, [x, y], :linear, :pos, sz)
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (x <= this).canon_form()
  append!(canon_constr_array, (y <= this).canon_form())
  # Add back the correct vexity
  this.vexity = :convex

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.max(x.evaluate(), y.evaluate())
  return this
end

# TODO: Max across dimensions

max(x::AbstractCvxExpr, y::Value) = max(x, convert(CvxExpr, y))
max(y::Value, x::AbstractCvxExpr) = max(x, convert(CvxExpr, y))
