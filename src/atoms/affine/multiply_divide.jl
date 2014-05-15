export *, /

# Utilities for handling vexity and sign for multiplication/division
function promote_sign(x::Constant, y::AbstractCvxExpr)
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

function promote_vexity(x::Constant, y::AbstractCvxExpr)
  if y.vexity == :linear
    return :linear
  elseif x.sign == :pos || x.size == :zero
    return y.vexity
  elseif x.sign == :neg
    return reverse_vexity(y)
  else
    error("expression not DCP compliant")
  end
end


function *(x::Constant, y::Constant)
  # TODO this won't work once we extend constants to parameters
  return Constant(x.value * y.value)
end

function *(x::Constant, y::AbstractCvxExpr)

  if x.size[2] != y.size[1] && x.size == (1,1)
    x = Constant(speye(y.size[1])*x.value[1], x.sign)
  end

  if x.size[2] == y.size[1]
    sz = (x.size[1], y.size[2])
    this = CvxExpr(:*, [x, y], promote_vexity(x, y), promote_sign(x, y), sz)
    vectorized_mul = kron(speye(sz[2]), x.value)

    coeffs = VecOrMatOrSparse[speye(get_vectorized_size(sz)), -vectorized_mul]
    vars = [this.uid, y.uid]
    constant = spzeros(get_vectorized_size(sz), 1)
    canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]

    this.canon_form = ()->append!(canon_constr_array, y.canon_form())
    return this

  elseif y.size == (1,1)
    this = CvxExpr(:*, [x, y], promote_vexity(x, y), promote_sign(x, y), x.size)

    coeffs = VecOrMatOrSparse[speye(get_vectorized_size(x.size)), -sparse(vec(x.value))]
    vars = [this.uid, y.uid]
    constant = spzeros(get_vectorized_size(x.size), 1)
    canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]

    this.canon_form = ()->append!(canon_constr_array, y.canon_form())
    return this

  else
    error("size of arguments cannot be multiplied; got $(x.size),$(y.size)")
  end
end

function *(x::AbstractCvxExpr, y::Constant)
  if y.size == (1, 1) || x.size == (1, 1)
    return y * x
  end

  sz = (x.size[1], y.size[2])
  vectorized_mul = kron(y.value', speye(sz[1]))

  this = CvxExpr(:*, [x, y], promote_vexity(y, x), promote_sign(y, x), sz)

  coeffs = VecOrMatOrSparse[speye(get_vectorized_size(sz)), -vectorized_mul]
  vars = [this.uid, x.uid]
  constant = spzeros(get_vectorized_size(sz), 1)
  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]

  this.canon_form = ()->append!(canon_constr_array, x.canon_form())
  return this
end

*(x::AbstractCvxExpr, y::Value) = *(x, convert(CvxExpr, y))
*(x::Value, y::AbstractCvxExpr) = *(convert(CvxExpr, x), y)

# Only division by constant scalars is allowed
function inv(y::Constant)
  # determine size
  if y.size != (1, 1)
    error("only division by constant scalars is allowed.")
  end
  return Constant(1 / y.value)
end

/(x::AbstractCvxExpr, y::Constant) = *(x, inv(y))
/(x::AbstractCvxExpr, y::Number) = *(x, 1 / y)
