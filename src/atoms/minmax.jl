export min, max, pos, neg

function max(x::AbstractCvxExpr)
  # TODO: handle vexities and signs
  this = CvxExpr(:max, [x], x.vexity, :pos, (1, 1))

  # -this tries to find this' canon_form
  this.canon_form = ()->[]
  canon_constr_array = (-this <= x).canon_form()
  append!(canon_constr_array, (x <= this).canon_form())
  this.canon_form = ()->canon_constr_array
  return this
end

max(x::AbstractCvxExpr, y::Value) = max(x, convert(CvxExpr, y))
max(y::Value, x::AbstractCvxExpr) = max(x, convert(CvxExpr, y))

function min(x::AbstractCvxExpr)
  return -max(-x)
end
