export convert, promote_value, promote_rule, promote_size, promote_vexity, promote_sign, print_debug
export reverse_vexity, reverse_sign

### Conversion and promotion
# TODO: The difference between conversion and promotion is messy.
function convert(::Type{CvxExpr},x)
  if typeof(x) == CvxExpr
    return x
  else
    return Constant(x)
  end
end

# TODO: WTF is going on
promote_rule(::Type{CvxExpr}, ::Type{AbstractArray}) = CvxExpr
promote_rule(::Type{CvxExpr}, ::Type{Number}) = CvxExpr

### Utility functions for arithmetic
# TODO: This function is broken, but in use
function promote_size(x::AbstractCvxExpr,y::AbstractCvxExpr)
  if x.size == y.size
    size = x.size
  elseif length(x.size) == 0
    size = y.size
  elseif length(y.size) == 0
    size = x.size
  elseif length(x.size) == 1 && Set(x.size[1],1) == Set(y.size...)
    size = y.size
  elseif length(y.size) == 1 && Set(y.size[1],1) == Set(x.size...)
    size = x.size
  elseif maximum(x.size) == 1
    size = y.size
  elseif maximum(y.size) == 1
    size = x.size
  else
    error("size of arguments must be the same; got $(x.size),$(y.size)")
  end
  return size
end

function promote_vexity(x::AbstractCvxExpr,y::AbstractCvxExpr)
  v1 = x.vexity; v2 = y.vexity; vexities = Set(v1,v2)
  if vexities == Set(:convex,:concave)
    error("expression not DCP compliant")
  elseif :convex in vexities
    return :convex
  elseif :concave in vexities
    return :concave
  elseif :linear in vexities
    return :linear
  else
    return :constant
  end
end

function promote_sign(x::AbstractCvxExpr,y::AbstractCvxExpr)
  s1 = x.sign; s2 = y.sign; signs = Set(s1,s2)
  if :any in signs || signs == Set(:pos,:neg)
    return :any
  else # then s1==s2
    return s1
  end
end

function promote_value(x::Value, sz::Int64)
  if size(x, 1) < sz
    return ones(sz, 1) * x
  end
  return x
end

function reverse_vexity(x::AbstractCvxExpr)
  vexity = x.vexity
  if vexity == :convex
    return :concave
  elseif vexity == :concave
    return :convex
  else
    return vexity
  end
end

function reverse_sign(x::AbstractCvxExpr)
  sign = x.sign
  if sign == :neg
    return :pos
  elseif sign == :pos
    return :neg
  else
    return sign
  end
end

function print_debug(debug, args...)
  if (debug)
    println(args)
  end
end
