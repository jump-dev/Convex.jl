import Base.abs

export norm, abs, norm_inf

function norm_inf(x::AbstractCvxExpr)
  # TODO: handle vexities and signs
  this = CvxExpr(:norm_inf, [x], x.vexity, :pos, (1, 1))

  # 'x <= this' will try to find the canon_form for 'this', so we need to initialize it
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (x <= this).canon_form()
  append!(canon_constr_array, (-this <= x).canon_form())
  this.canon_form = ()->canon_constr_array
  return this
end


# TODO: Everything
function norm(x::AbstractCvxExpr, p = 2)
  norm_map = {1=>:norm1, 2=>:norm2, :inf=>:norm_inf, :nuc=>:norm_nuc}
  norm_type = norm_map[p]
  if x.vexity == :constant
    return CvxExpr(norm_type,[x],:constant,:pos,(1,1))
  elseif x.vexity == :linear
    return CvxExpr(norm_type,[x],:convex,:pos,(1,1))
  elseif x.vexity == :convex && x.sign == :pos
    return CvxExpr(norm_type,[x],:convex,:pos,(1,1))
  elseif x.vexity == :concave && x.sign == :neg
    return CvxExpr(norm_type,[x],:convex,:pos,(1,1))
  else
    error("norm(x) is not DCP compliant when x has curvature $(x.vexity) and sign $(x.sign)")
  end
end

### elementwise

function abs(x::AbstractCvxExpr)
  if x.vexity == :constant
    this = CvxExpr(:abs,[x],:constant,:pos,x.size)
  elseif x.vexity == :linear
    if x.sign == :pos
      this = CvxExpr(:abs,[x],:linear,:pos,x.size)
    elseif x.sign == :neg
      this = CvxExpr(:abs,[x],:linear,:pos,x.size)
    else
      this = CvxExpr(:abs,[x],:convex,:pos,x.size)
    end
  elseif x.vexity == :convex && x.sign == :pos
    this = CvxExpr(:abs,[x],:convex,:pos,x.size)
  elseif x.vexity == :concave && x.sign == :neg
    this = CvxExpr(:abs,[x],:convex,:pos,x.size)
  else
    error("abs(x) is not DCP compliant when x has curvature $(x.vexity) and sign $(x.sign)")
  end

  println(this.vexity)
  # 'x <= this' will try to find the canon_form for 'this', so we need to initialize it
  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = (x <= this).canon_form()
  append!(canon_constr_array, (-this <= x).canon_form())
  this.canon_form = ()->canon_constr_array
  return this
end
