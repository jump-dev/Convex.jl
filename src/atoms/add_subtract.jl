export +, -


function -(x::Constant)
  # TODO this won't work once we extend constants to parameters
  return Constant(-x.value)
end


function -(x::AbstractCvxExpr)
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
