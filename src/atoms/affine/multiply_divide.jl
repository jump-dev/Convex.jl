export *, /

# TODO: Handle .* and ./

# For constants, just multiply their values
function *(x::Constant, y::Constant)
  return Constant(x.value * y.value)
end

# Multiplication of an AbstractCvxExpr `y` with a Constant `x` to the left
# 1. If the constant is 1 x 1, we simply multiply it by ones(size of `y`)
# 2. If the sizes are suitable for multiplication, we then need to handle matrix
# variables. This is done using Kroneckor products. See details here:
# http://en.wikipedia.org/wiki/Vectorization_(mathematics)
# 3. If size of `y` is 1 x 1, we remember that x * y will ultimately be vectorized
# This knowledge allows us to simply have the canonical form in this case as
# If w = x * y, I * w - vectorized(x) * y = 0
function *(x::Constant, y::AbstractCvxExpr)

  if x.size[2] != y.size[1] && x.size == (1, 1)
    x = Constant(speye(y.size[1]) * x.value[1], x.sign)
  end

  if x.size[2] == y.size[1]
    sz = (x.size[1], y.size[2])
    this = CvxExpr(:*, [x, y], promote_vexity_multiply(x, y), promote_sign_multiply(x, y), sz)

    # Kronecker product for vectorized multiplication
    vectorized_mul = kron(speye(sz[2]), x.value)

    coeffs = VecOrMatOrSparse[speye(get_vectorized_size(sz)), -vectorized_mul]
    vars = [this.uid, y.uid]
    constant = spzeros(get_vectorized_size(sz), 1)
    canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]

    this.canon_form = ()->append!(canon_constr_array, y.canon_form())

    this.evaluate = ()->x.evaluate() * y.evaluate()
    return this

  elseif y.size == (1, 1)
    this = CvxExpr(:*, [x, y], promote_vexity_multiply(x, y), promote_sign_multiply(x, y), x.size)

    coeffs = VecOrMatOrSparse[speye(get_vectorized_size(x.size)), -sparse(vec(x.value))]
    vars = [this.uid, y.uid]
    constant = spzeros(get_vectorized_size(x.size), 1)
    canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]

    this.canon_form = ()->append!(canon_constr_array, y.canon_form())

    this.evaluate = ()->x.evaluate() * y.evaluate()
    return this
  else
    error("Size of arguments cannot be multiplied; got $(x.size), $(y.size)")
  end
end

# If either of `x` or `y` is 1 x 1, we simply return `y * x`, which has been
# discussed above. Otherwise, we perform canonicalization similar to 2. above
function *(x::AbstractCvxExpr, y::Constant)
  if y.size == (1, 1) || x.size == (1, 1)
    return y * x
  end

  sz = (x.size[1], y.size[2])
  vectorized_mul = kron(y.value', -speye(sz[1]))

  this = CvxExpr(:*, [x, y], promote_vexity_multiply(y, x), promote_sign_multiply(y, x), sz)

  coeffs = VecOrMatOrSparse[speye(get_vectorized_size(sz)), vectorized_mul]
  vars = [this.uid, x.uid]
  constant = spzeros(get_vectorized_size(sz), 1)
  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]

  this.canon_form = ()->append!(canon_constr_array, x.canon_form())
  this.evaluate = ()->x.evaluate() * y.evaluate()
  return this
end

function *(x::AbstractCvxExpr, y::AbstractCvxExpr)
  error("Multiplication between two variable expressions is not DCP compliant. Perhaps
        you want to look into norm, or quad_form?")
end

*(x::AbstractCvxExpr, y::Value) = *(x, convert(CvxExpr, y))
*(x::Value, y::AbstractCvxExpr) = *(convert(CvxExpr, x), y)

# Only division by constant scalars is allowed
function inv(y::Constant)
  if y.size != (1, 1)
    error("Only division by constant scalars is allowed")
  end
  return Constant(1 / y.value)
end

/(x::AbstractCvxExpr, y::Constant) = *(x, inv(y))
/(x::AbstractCvxExpr, y::Number) = *(x, 1 / y)
