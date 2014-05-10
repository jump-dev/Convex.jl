export +, -

### Utilities for handling vexity and sign for addition/subtraction

function promote_vexity(x::AbstractCvxExpr, y::AbstractCvxExpr)
  vexities = Set(x.vexity, y.vexity)
  if vexities == Set(:convex, :concave)
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

function promote_sign(x::AbstractCvxExpr, y::AbstractCvxExpr)
  signs = Set(x.sign, y.sign)
  if :any in signs || signs == Set(:pos,:neg)
    return :any
  elseif x.sign == :zero
    return y.sign
  elseif y.sign == :zero
    return x.sign
  else
    return x.sign
  end
end

### Unary Negation

function -(x::Constant)
  # TODO this won't work once we extend constants to parameters
  return Constant(-x.value)
end

function -(x::AbstractCvxExpr)
  this = CvxExpr(:-, [x], reverse_vexity(x), reverse_sign(x), x.size)

  if x.vexity == :constant
    this.canon_form = ()->CanonicalConstr[]
  else
    coeffs = VecOrMatOrSparse[speye(get_vectorized_size(x)), speye(get_vectorized_size(x))]

    vars = [this.uid, x.uid]
    constant = zeros(get_vectorized_size(x))
    canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true)]
    append!(canon_constr_array, x.canon_form())

    this.canon_form = ()->canon_constr_array
  end

  return this
end

### Binary Addition/Subtraction

function +(x::AbstractCvxExpr, y::AbstractCvxExpr)
  x, y = promote_for_add(x, y)
  this = CvxExpr(:+, [x, y], promote_vexity(x, y), promote_sign(x, y), x.size)

  coeffs = VecOrMatOrSparse[-speye(get_vectorized_size(x)),
      -speye(get_vectorized_size(x)),
      speye(get_vectorized_size(x))]
  vars = [x.uid, y.uid, this.uid]
  constant = zeros(get_vectorized_size(x))

  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true)]

  append!(canon_constr_array, x.canon_form())
  append!(canon_constr_array, y.canon_form())

  this.canon_form = ()->begin
    return canon_constr_array
  end

  return this
end


function +(x::Constant, y::Constant)
  # TODO this won't work once we extend constants to parameters
  x, y = promote_for_add(x, y)
  this = Constant(x.value + y.value)
  return this
end


function +(x::AbstractCvxExpr, y::Constant)
  x, y = promote_for_add(x, y)
  this = CvxExpr(:+, [x, y], promote_vexity(x, y), promote_sign(x, y), x.size)

  coeffs = VecOrMatOrSparse[-speye(get_vectorized_size(x)), speye(get_vectorized_size(x))]
  vars = [x.uid, this.uid]
  constant = vec(y.value)
  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true)]

  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->begin
    return canon_constr_array
  end
  return this
end


+(y::Constant, x::AbstractCvxExpr) = +(x::AbstractCvxExpr, y::Constant)
+(x::AbstractCvxExpr, y::Value) = +(x, convert(CvxExpr, y))
+(x::Value, y::AbstractCvxExpr) = +(y, convert(CvxExpr, x))
-(x::AbstractCvxExpr, y::AbstractCvxExpr) = +(x, -y)
-(x::AbstractCvxExpr, y::Value) = +(x, -y)
-(x::Value, y::AbstractCvxExpr) = +(-y, x)
