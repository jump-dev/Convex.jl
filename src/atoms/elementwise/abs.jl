import Base.abs
export abs

# Absolute value of the expression
# Canonical constraint is -w <= x <= w if w = abs(x)
function abs(x::AbstractCvxExpr)
  if x.vexity == :constant
    this = CvxExpr(:abs, [x], :constant, :pos, x.size)
  elseif x.vexity == :linear
    if x.sign == :pos
      this = CvxExpr(:abs, [x], :linear, :pos, x.size)
    elseif x.sign == :neg
      this = CvxExpr(:abs, [x], :linear, :pos, x.size)
    else
      this = CvxExpr(:abs, [x], :convex, :pos, x.size)
    end
  elseif x.vexity == :convex && x.sign == :pos
    this = CvxExpr(:abs, [x], :convex, :pos, x.size)
  elseif x.vexity == :concave && x.sign == :neg
    this = CvxExpr(:abs, [x], :convex, :pos, x.size)
  else
    error("abs(x) is not DCP compliant when x has curvature $(x.vexity) and sign $(x.sign)")
  end

  # Change vexity temporarily to allow <=
  vexity = this.vexity
  this.vexity = :linear

  # 'x <= this' will try to find the canon_form for 'this', so we need to initialize it
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (x <= this).canon_form()
  append!(canon_constr_array, (-this <= x).canon_form())

  # Fix vexity
  this.vexity = vexity

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.abs(x.evaluate())
  return this
end

