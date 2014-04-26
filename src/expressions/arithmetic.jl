export +,-,*,/


# multiple arguments
for op = (:promote_vexity, :promote_sign, :promote_shape)
  @eval ($op)(arg1,arg2,arg3,args...) = length(args)==0 ? ($op)(($op)(arg1,arg2),arg3) : ($op)(($op)(arg1,arg2),arg3,args...)
end

# TODO: We shouldn't have a size() and an x.size
function size(x::AbstractCvxExpr)
  return x.size
end

### Unary functions on expressions

function -(x::Constant)
  # TODO this won't work once we extend constants to parameters
  return Constant(-x.value)
end

function -(x::AbstractCvxExpr)
  # evalfn = ()->
  #   begin
  #     ret_val = x.evalfn()
  #     if ret_val == nothing
  #       return nothing
  #     else
  #       return -ret_val
  #     end
  #   end

  # TODO: Implement matrix shit, only works for vectors

  this = CvxExpr(:-, [x], reverse_vexity(x), reverse_sign(x), x.size)

  # TODO: Not Any
  if x.vexity == :constant
    this.canon_form = ()->Any[]
  else
    canon_constr_array = Any[{
      :coeffs => Any[speye(x.size[1]), speye(x.size[1])],
      :vars => [this.uid(), x.uid()],
      :constant => zeros(x.size),
      :is_eq => true
    }]

    this.canon_form = ()->append!(canon_constr_array, x.canon_form())
  end

  return this
end

### Arithmetic for expressions
function +(x::AbstractCvxExpr, y::AbstractCvxExpr)
  promote_for_add!(x, y)
  this = CvxExpr(:+, [x, y], promote_vexity(x, y), promote_sign(x, y), x.size)

  # TODO: Not Any. Also deal with matrix variables
  canon_constr_array = Any[{
    :coeffs => Any[-speye(x.size[1]), -speye(x.size[1]), speye(x.size[1])],
    :vars => [x.uid(), y.uid(), this.uid()],
    :constant => zeros(x.size),
    :is_eq => true
  }]

  this.canon_form = ()->begin
    append!(canon_constr_array, x.canon_form())
    append!(canon_constr_array, y.canon_form())
    return canon_constr_array
  end

  return this
end

function +(x::Constant, y::Constant)
  # TODO this won't work once we extend constants to parameters
  this = Constant(x.value + y.value)
  return this
end

function +(x::AbstractCvxExpr, y::Constant)
  promote_for_add!(x, y)
  this = CvxExpr(:+, [x, y], promote_vexity(x, y), promote_sign(x, y), x.size)

  # TODO: Not Any. Also deal with matrix variables
  canon_constr_array = Any[{
    :coeffs => Any[-speye(x.size[1]), speye(x.size[1])],
    :vars => [x.uid(), this.uid()],
    # TODO we'll need to cache references to constants/parameters in the future
    :constant => y.value,
    :is_eq => true
  }]

  this.canon_form = ()->begin
    append!(canon_constr_array, x.canon_form())
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
  size = (x.size[1], y.size[2])

  # determine sign
  signs = Set(x.sign,y.sign)
  if :any in signs
    sign = :any
  elseif :zero in signs
    sign = :zero
  elseif length(signs) == 1
    sign = :pos
  else
    sign = :neg
  end

  this = CvxExpr(:*, [x, y], y.vexity, sign, size)

  canon_constr_array = Any[{
    # TODO we'll need to cache references to parameters in the future
    :coeffs => Any[speye(size[1]), -x.value],
    :vars => [this.uid(), y.uid()],
    :constant => spzeros(size...),
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

  size = (x.size[1], y.size[2])

  # determine sign
  signs = Set(x.sign,y.sign)
  if :any in signs
    sign = :any
  elseif :zero in signs
    sign = :zero
  elseif length(signs) == 1
    sign = :pos
  else
    sign = :neg
  end

  this = CvxExpr(:*, [x, y], y.vexity, sign, size)

  canon_constr_array = Any[{
    # TODO we'll need to cache references to parameters in the future
    :coeffs => Any[speye(size[1]), -y.value],
    :vars => [this.uid(), x.uid()],
    :constant => spzeros(size...),
    :is_eq => true
  }]

  this.canon_form = ()->append!(canon_constr_array, y.canon_form())
  return this
end

*(x::AbstractCvxExpr,y::Value) = *(x,convert(CvxExpr,y))
*(x::Value,y::AbstractCvxExpr) = *(convert(CvxExpr,x),y)

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
