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
  Constant(-x.value)
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

  this = CvxExpr(:-, [x], reverse_vexity(x), reverse_sign(x), size(x))

  # TODO: Not Any
  if x.vexity == :constant
    this.canon_form = ()->Any[]
  else
    canon_constr_array = Any[{
      :coeffs => Any[speye(x.size[1]), speye(x.size[1])],
      :vars => [this.uid(), x.uid()],
      :constant => 0,
      :is_eq => true
    }]

    this.canon_form = ()->append!(canon_constr_array, x.canon_form())
  end

  return this
end

### Arithmetic for expressions
function +(x::AbstractCvxExpr, y::AbstractCvxExpr)
  this = CvxExpr(:+, [x, y], promote_vexity(x, y), promote_sign(x, y), promote_size(x, y))

  if x.size[1] != y.size[1]
    error("TODO: Implement overloading of +")
  end

  sz = maximum([x.size[1], y.size[1]])
  # TODO: Not Any. Also deal with matrix variables
  canon_constr_array = Any[{
    :coeffs => Any[-speye(sz), -speye(sz), speye(sz)],
    :vars => [x.uid(), y.uid(), this.uid()],
    :constant => 0,
    :is_eq => true
  }]

  this.canon_form = ()->begin
    append!(canon_constr_array, x.canon_form())
    append!(canon_constr_array, y.canon_form())
    return canon_constr_array
  end

  return this
end

+(x::AbstractCvxExpr, y::Value) = +(x, convert(CvxExpr, y))
+(x, y::AbstractCvxExpr) = +(y, convert(CvxExpr, x))
-(x::AbstractCvxExpr, y::AbstractCvxExpr) = +(x, -y)
-(x::AbstractCvxExpr, y) = +(x, -y)
-(x, y::AbstractCvxExpr) = +(-y, x)

function *(x::AbstractCvxExpr, y::AbstractCvxExpr)
  # determine vexity
  # TODO: As of now, we only handle linear vexity
  if Set(x.vexity, y.vexity) == Set(:constant,:linear)
    vexity = :linear
  elseif x.vexity == :constant && y.vexity == :constant
    vexity = :constant
  elseif x.vexity == :constant && x.sign == :pos
    vexity = y.vexity
  elseif x.vexity == :constant && x.sign == :neg
    vexity = reverse_vexity(y)
  elseif y.vexity == :constant && y.sign == :pos
    vexity = x.vexity
  elseif y.vexity == :constant && y.sign == :neg
    vexity = reverse_vexity(x)
  else
    error("expression not DCP compliant; got $x * $y")
  end

  # TODO: wtf is going on
  # determine size
  # multiplication by a scalar
  if length(x.size) == 0
    size = y.size
  elseif length(y.size) == 0
    size = x.size
  # matrix multiplication
  elseif length(y.size) == 2 && length(x.size) == 2 && x.size[2] == y.size[1]
    size = (x.size[1], y.size[2])
  # cast x to a row vector
  elseif length(y.size) == 2 && length(x.size) == 1 && x.size[1] == y.size[1]
    size = (1, y.size[2])
  # cast x to a column vector
  elseif length(y.size) == 2 && length(x.size) == 1 && y.size[1] == 1
    size = (x.size[1], y.size[2])
  # cast y to a row vector
  elseif length(x.size) == 2 && length(y.size) == 1 && x.size[2] == 1
    size = (1, y.size[1])
  # cast y to a column vector
  elseif length(x.size) == 2 && length(y.size) == 1 && x.size[2] == y.size[1]
    size = (x.size[1], 1)
  else
    error("inner dimensions must agree; got $(x.size) * $(y.size)")
  end

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

  this = CvxExpr(:*, [x, y], vexity, sign, size)

  lin, cons = x.vexity == :linear ? (x, y) : (y, x)

  canon_constr_array = Any[{
    :coeffs => Any[speye(size[1]), -cons.value],
    :vars => [this.uid(), lin.uid()],
    :constant => 0,
    :is_eq => true
  }]

  this.canon_form = ()->append!(canon_constr_array, lin.canon_form())

  return this
end
*(x::AbstractCvxExpr,y::Value) = *(x,convert(CvxExpr,y))
*(x::Value,y::AbstractCvxExpr) = *(convert(CvxExpr,x),y)

# TODO: implement inv
# only division by constant scalars is allowed, but we need to provide it for use with parameters, maybe
function inv(y::AbstractCvxExpr)
  # determine size
  if y.size in Set((),(1),(1,1)) && y.vexity == :constant
    size = y.size
  else
    error("only division by constant scalars is allowed.")
  end

  CvxExpr(:/,[1,y],:constant,y.sign,size)
end
/(x::AbstractCvxExpr,y::AbstractCvxExpr) = *(x,inv(y))
/(x::AbstractCvxExpr,y::Number) = *(x,1/y)
# hcat and vcat (only vcat exists in cvxpy)
