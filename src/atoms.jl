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
transpose(x::AbstractCvxExpr) = CvxExpr(:transpose,[x],x.vexity,x.sign,reverse(x.size))
ctranspose(x::AbstractCvxExpr) = transpose(x)
getindex(x::AbstractCvxExpr,index...) = CvxExpr(:index,[x,index...],x.vexity,x.sign,size)
# this breaks the syntax [x,y] to concatenate lists of expressions as in arguments to atoms
# extend to *args
# vcat(args::Array{AbstractCvxExpr}) = CvxExpr(:vstack,args,promote_vexity([a.vexity for a in args]...),promote_sign([a.sign for a in args]...),sizes...)

function dot(x::AbstractCvxExpr, y::AbstractCvxExpr)
	# TODO: error checks
	return CvxExpr(:dot, [x y], :linear, :any, (1, 1))
end

function norm(x::AbstractCvxExpr, p = 2)
	norm_map = {1=>:norm1, 2=>:norm2, :inf=>:norm_inf, :nuc=>:norm_nuc}
	norm_type = norm_map[p]
	if x.vexity == :constant
		return CvxExpr(norm_type,[x],:constant,:pos,(1,1))
	elseif x.vexity == :linear
		return CvxExpr(norm_type,[x],:convex,:pos,(1,1))
	elseif x.vexity == :convex && x.sign == :pos
		return CvxExpr(norm_type,[x],:convex,:pos,(1,1))
	elseif x.vexity == :concave && x.sign == :neg
		return CvxExpr(norm_type,[x],:convex,:pos,(1,1))
	else
		error("norm(x) is not DCP compliant when x has curvature $(x.vexity) and sign $(x.sign)")
	end
end

### elementwise
import Base.abs
function abs(x::AbstractCvxExpr)
	if x.vexity == :constant
		return CvxExpr(:abs,[x],:constant,:pos,x.size)
	elseif x.vexity == :linear
		if x.sign == :pos
			return CvxExpr(:abs,[x],:linear,:pos,x.size)
		elseif x.sign == :neg
			return CvxExpr(:abs,[x],:linear,:pos,x.size)
		else
			return CvxExpr(:abs,[x],:convex,:pos,x.size)
		end
	elseif x.vexity == :convex && x.sign == :pos
		return CvxExpr(:abs,[x],:convex,:pos,x.size)
	elseif x.vexity == :concave && x.sign == :neg
		return CvxExpr(:abs,[x],:convex,:pos,x.size)
	else
		error("abs(x) is not DCP compliant when x has curvature $(x.vexity) and sign $(x.sign)")
	end
end
