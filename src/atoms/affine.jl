import Base.sum, Base.abs, Base.sqrt, Base.log
export dot
# TODO
# * max, min

# this breaks the syntax [x,y] to concatenate lists of expressions as in arguments to atoms
# extend to *args
# vcat(args::Array{AbstractCvxExpr}) = CvxExpr(:vstack,args,promote_vexity([a.vexity for a in args]...),promote_sign([a.sign for a in args]...),sizes...)

function transpose(x::AbstractCvxExpr)
  sz = get_vectorized_size(x.size)
  coeffs = spzeros(sz, sz)
  num_rows = x.size[1]
  num_cols = x.size[2]
  for r = 1:num_rows
    for c = 1:num_cols
      i = (c - 1) * num_rows + r
      j = (r - 1) * num_cols + c
      coeffs[i, j] = 1.0
    end
  end

  this = CvxExpr(:transpose, [x], x.vexity, x.sign, (x.size[2], x.size[1]))
  canon_constr_array = Any[{
    :coeffs => Any[speye(sz), -coeffs],
    :vars => [x.uid(), this.uid()],
    :constant => spzeros(sz, 1),
    :is_eq => true
  }]
  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array
  return this
end

ctranspose(x::AbstractCvxExpr) = transpose(x)

function transpose(x::Constant)
  return Constant(x.value')
end

ctranspose(x::Constant) = transpose(x)

function dot(x::AbstractCvxExpr, y::AbstractCvxExpr)
	if Set(x.vexity, y.vexity) != Set(:constant, :linear)
		error("TODO: Only LP for now")
		# TODO: Also check x and y are vectors
	end

	lin, cons = x.vexity == :linear ? (x, y) : (y, x)
	@assert cons.size == lin.size

	this = CvxExpr(:dot, [x y], :linear, :any, (1, 1))

	# TODO: Don't do any, make more efficient
	canon_constr_array = Any[{
    :coeffs => Any[1.0, -cons.value'],
    :vars => [this.uid(), lin.uid()],
    :constant => 0,
    :is_eq => true
  }]

  this.canon_form = ()->append!(canon_constr_array, lin.canon_form())
  return this
end

dot(x::Value, y::AbstractCvxExpr) = dot(convert(CvxExpr, x), y)
dot(x::AbstractCvxExpr, y::Value) = dot(x, convert(CvxExpr, y))

function sum(x::AbstractCvxExpr)
  if x.size == (1, 1)
    return x
  elseif x.size[1] == 1
    return x * ones(x.size[2], 1)
  else
    return ones(1, x.size[1]) * x
  end
end

function sum(x::Constant)
  return Constant(sum(x.value))
end

