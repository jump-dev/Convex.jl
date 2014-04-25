import Base.sum, Base.abs, Base.sqrt, Base.log
# export transpose, ctranspose, kl_div, lambda_min, lambda_max, log_det, norm, quad_form, quad_over_lin, abs, pos, square, sum, log, sqrt, inv_pos
export dot
# TODO
# * max, min

# cvxpy atoms
#['abs', 'affine', 'atom', 'elementwise', 'geo_mean', 'inv_pos', 'lambda_max', 'lambda_min', 'log', 'max', 'min', 'neg', 'nonlinear', 'norm', 'norm1', 'norm2', 'normInf', 'normNuc', 'norm_inf', 'norm_nuc', 'pos', 'quad_form', 'quad_over_lin', 'sigma_max', 'sqrt', 'square', 'sum', 'vstack']
# cvxpy atoms still to go
#['abs', 'affine', 'atom', 'elementwise', 'geo_mean', 'inv_pos', 'lambda_max', 'lambda_min', 'log', 'max', 'min', 'neg', 'nonlinear', 'norm', 'norm1', 'norm2', 'normInf', 'normNuc', 'norm_inf', 'norm_nuc', 'pos', 'quad_form', 'quad_over_lin', 'sigma_max', 'sqrt', 'square', 'sum', 'vstack']

# slicing atoms
# get the right names for these
# transpose(x::AbstractCvxExpr) = CvxExpr(:transpose,[x],x.vexity,x.sign,reverse(x.size))
# ctranspose(x::AbstractCvxExpr) = transpose(x)
# getindex(x::AbstractCvxExpr,index...) = CvxExpr(:index,[x,index...],x.vexity,x.sign,size(index))
# this breaks the syntax [x,y] to concatenate lists of expressions as in arguments to atoms
# extend to *args
# vcat(args::Array{AbstractCvxExpr}) = CvxExpr(:vstack,args,promote_vexity([a.vexity for a in args]...),promote_sign([a.sign for a in args]...),sizes...)

function dot(x::AbstractCvxExpr, y::AbstractCvxExpr)
	if Set(x.vexity, y.vexity) != Set(:constant, :linear)
		error("TODO: Only LP for now")
		# TODO: Also check x and y are vectors
	end

	lin, cons = x.vexity == :linear ? (x, y) : (y, x)
  # println(cons.size)
  # println(lin.size)
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
  return dot(x, ones(x.size))
end
