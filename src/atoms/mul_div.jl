export *, /

function *(x::Constant, y::Constant)
  # TODO this won't work once we extend constants to parameters
  this = Constant(x.value * y.value)
  return this
end

function *(x::Constant, y::AbstractCvxExpr)
  if y.vexity != :linear
    error("Only LPs allowed for now")
  end

  promote_for_mul!(x, y)
  sz = (x.size[1], y.size[2])

  # determine sign
  signs = Set(x.sign, y.sign)
  if :any in signs
    sign = :any
  elseif :zero in signs
    sign = :zero
  elseif length(signs) == 1
    sign = :pos
  else
    sign = :neg
  end

  this = CvxExpr(:*, [x, y], y.vexity, sign, sz)

  # TODO: Change to speye once Julia fixes kron bug
  vectorized_mul = kron(eye(sz[2]), x.value)

  canon_constr_array = Any[{
    # TODO we'll need to cache references to parameters in the future
    :coeffs => Any[speye(get_vectorized_size(sz)), -vectorized_mul],
    :vars => [this.uid(), y.uid()],
    :constant => spzeros(get_vectorized_size(sz), 1),
    :is_eq => true
  }]

  this.canon_form = ()->append!(canon_constr_array, x.canon_form())
  return this
end

function *(x::AbstractCvxExpr, y::Constant)

  if y.size == (1, 1)
    return y * x
  end

  if x.vexity != :linear
    error("Only LPs allowed for now")
  end

  lhs = x
  if x.size[2] != y.size[1]
    if x.size != (1, 1)
      error("Can't promote size of variable with size $(x.size) to $(y.size).")
    end
    return y * x
  end

  sz = (x.size[1], y.size[2])

  # determine sign
  signs = Set(x.sign, y.sign)
  if :any in signs
    sign = :any
  elseif :zero in signs
    sign = :zero
  elseif length(signs) == 1
    sign = :pos
  else
    sign = :neg
  end

  # TODO: Change to speye once Julia fixes kron bug
  vectorized_mul = kron(y.value', eye(sz[1]))

  this = CvxExpr(:*, [x, y], y.vexity, sign, sz)

  canon_constr_array = Any[{
    # TODO we'll need to cache references to parameters in the future
    :coeffs => Any[speye(get_vectorized_size(sz)), -vectorized_mul],
    :vars => [this.uid(), x.uid()],
    :constant => spzeros(get_vectorized_size(sz), 1),
    :is_eq => true
  }]

  this.canon_form = ()->append!(canon_constr_array, y.canon_form())
  return this
end

*(x::AbstractCvxExpr, y::Value) = *(x, convert(CvxExpr, y))
*(x::Value, y::AbstractCvxExpr) = *(convert(CvxExpr, x), y)

# TODO: implement inv
# only division by constant scalars is allowed, but we need to provide it for use with parameters, maybe
function inv(y::Constant)
  # determine size
  if y.size != (1,1)
    error("only division by constant scalars is allowed.")
  end

  return Constant(1/y.value)
end
/(x::AbstractCvxExpr,y::Constant) = *(x,inv(y))
/(x::AbstractCvxExpr,y::Number) = *(x,1/y)
# hcat and vcat (only vcat exists in cvxpy)
