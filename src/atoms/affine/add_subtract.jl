export +, -

### Utilities for handling vexity and sign for addition/subtraction

function promote_vexity(x::AbstractCvxExpr, y::AbstractCvxExpr)
  vexities = Set({x.vexity, y.vexity})
  if vexities == Set({:convex, :concave})
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

# Returns the sign of x + y. Returns :any if we can't determine the sign
function promote_sign(x::AbstractCvxExpr, y::AbstractCvxExpr)
  signs = Set({x.sign, y.sign})
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
  return Constant(-x.value)
end

# If the variable created for -x is w, then the canonical form is
# - (I * w) + I * -x = 0
function -(x::AbstractCvxExpr)
  this = CvxExpr(:-, [x], reverse_vexity(x), reverse_sign(x), x.size)

  if x.vexity == :constant
    this.canon_form = ()->CanonicalConstr[]
  else
    coeffs = VecOrMatOrSparse[speye(get_vectorized_size(x)), speye(get_vectorized_size(x))]

    vars = [this.uid, x.uid]
    constant = zeros(get_vectorized_size(x))
    canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]
    append!(canon_constr_array, x.canon_form())

    this.canon_form = ()->canon_constr_array
  end

  this.evaluate = ()->-x.evaluate()

  return this
end

### Binary Addition/Subtraction

# Let w = x + y be the variable created
# Addition if:
#
# 1. x and y are the same size
# # Canonical form is - x - y + w = 0
# 2. x is (1, 1) and y is a vector/ matrix of size m x n (Vectorized as of size mn x 1)
# # Canonical form is -ones(mn, 1) * x - y + w = 0
# 3. y is (1, 1). We just return y + x which then falls into case 2
function +(x::AbstractCvxExpr, y::AbstractCvxExpr)
  if x.size != y.size

    if x.size == (1, 1)
      sz_y = get_vectorized_size(y)
      this = CvxExpr(:+, [x, y], promote_vexity(x, y), promote_sign(x, y), y.size)
      coeffs = VecOrMatOrSparse[-ones(sz_y, 1), -speye(sz_y), speye(sz_y)]
      vars = [x.uid, y.uid, this.uid]
      constant = zeros(sz_y)
    elseif y.size == (1, 1)
      return y + x
    else
      error("Can't add expressions of size $(x.size) and $(y.size)")
    end
  else
    sz = get_vectorized_size(y)
    this = CvxExpr(:+, [x, y], promote_vexity(x, y), promote_sign(x, y), y.size)
    coeffs = VecOrMatOrSparse[-speye(sz), -speye(sz), speye(sz)]
    vars = [x.uid, y.uid, this.uid]
    constant = zeros(sz)
  end

  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]

  append!(canon_constr_array, x.canon_form())
  append!(canon_constr_array, y.canon_form())

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->x.evaluate() + y.evaluate()

  return this
end

# Same rules as above. Except that if y is a scalar, we can simply multiply it by
# ones(...) to promote its size to what it should be.
function +(x::AbstractCvxExpr, y::Constant)
  if x.size != y.size && x.size == (1, 1)
    sz_y = get_vectorized_size(y)
    this = CvxExpr(:+, [x, y], promote_vexity(x, y), promote_sign(x, y), y.size)
    coeffs = VecOrMatOrSparse[-ones(sz_y, 1), speye(sz_y)]
  elseif x.size != y.size && y.size != (1, 1)
    error("Can't add expressions of size $(x.size) and $(y.size)")
  else
    if y.size == (1, 1)
      y = Constant(y.value * ones(x.size...), y.sign)
    end
    sz = get_vectorized_size(x)
    this = CvxExpr(:+, [x, y], promote_vexity(x, y), promote_sign(x, y), x.size)
    coeffs = VecOrMatOrSparse[-speye(sz), speye(sz)]
  end

  vars = [x.uid, this.uid]
  constant = vec(y.value)
  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]

  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->x.evaluate() + y.evaluate()
  return this
end

# Adding two constants doesn't require canonicalization
function +(x::Constant, y::Constant)
  this = Constant(x.value + y.value)
  return this
end

# Override addition since julia doesn't allow things like [1] + eye(4)
function +(x::Array{Number,}, y::Array{Number,})
  if x.size == (1, 1)
    return x[1] + y
  elseif y.size == (1, 1)
    return x + y[1]
  else
    return x + y
  end
end

+(y::Constant, x::AbstractCvxExpr) = +(x::AbstractCvxExpr, y::Constant)
+(x::AbstractCvxExpr, y::Value) = +(x, convert(CvxExpr, y))
+(x::Value, y::AbstractCvxExpr) = +(y, convert(CvxExpr, x))
-(x::AbstractCvxExpr, y::AbstractCvxExpr) = +(x, -y)
-(x::AbstractCvxExpr, y::Value) = +(x, -y)
-(x::Value, y::AbstractCvxExpr) = +(-y, x)
