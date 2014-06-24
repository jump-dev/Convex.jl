export promote_vexity_add, promote_sign_add, promote_vexity_multiply, promote_sign_multiply
export promote_vexity_norm

### Utilities for handling vexity and sign for add/subtract/hcat/vcat

function promote_vexity_add(vexity_one::Symbol, vexity_two::Symbol)
  vexities = {vexity_one, vexity_two}
  if :convex in vexities && :concave in vexities
    error("Expression not DCP compliant")
  elseif :convex in vexities
    return :convex
  elseif :concave in vexities
    return :concave
  elseif :affine in vexities
    return :affine
  else
    return :constant
  end
end

function promote_vexity_add(x::AbstractCvxExpr, y::AbstractCvxExpr)
  return promote_vexity_add(x.vexity, y.vexity)
end

# Returns the sign of x + y. Returns :any if we can't determine the sign
function promote_sign_add(sign_one::Symbol, sign_two::Symbol)
  signs = {sign_one, sign_two}
  if :any in signs || (:pos in signs && :neg in signs)
    return :any
  elseif sign_one == :zero
    return sign_two
  elseif sign_two == :zero
    return sign_one
  else
    return sign_one
  end
end

function promote_sign_add(x::AbstractCvxExpr, y::AbstractCvxExpr)
  return promote_sign_add(x.sign, y.sign)
end

### Utilities for handling vexity and sign for multiplication/division

function promote_sign_multiply(x::Constant, y::AbstractCvxExpr)
  if x.sign == :zero || y.sign == :zero
    return :zero
  elseif x.sign == :pos
    return y.sign
  elseif x.sign == :neg
    return reverse_sign(y)
  else
    return :any
  end
end

function promote_vexity_multiply(x::Constant, y::AbstractCvxExpr)
  if y.vexity == :affine
    return :affine
  elseif y.vexity == :constant || x.sign == :zero
    return :constant
  elseif x.sign == :pos
    return y.vexity
  elseif x.sign == :neg
    return reverse_vexity(y)
  else
    error("Expression not DCP compliant")
  end
end

### Utility for norm

function promote_vexity_norm(x::AbstractCvxExpr)
  if x.vexity == :constant
    return :constant
  elseif x.vexity == :affine
    return :convex
  elseif x.vexity == :convex && x.sign == :pos
    return :convex
  elseif x.vexity == :concave && x.sign == :neg
    return :convex
  else
    error("norm(x) is not DCP compliant when x has curvature $(x.vexity) and sign $(x.sign)")
  end
end

