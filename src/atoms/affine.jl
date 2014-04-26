import Base.sum, Base.abs, Base.sqrt, Base.log
export dot
# TODO
# * max, min

# this breaks the syntax [x,y] to concatenate lists of expressions as in arguments to atoms
# extend to *args
# vcat(args::Array{AbstractCvxExpr}) = CvxExpr(:vstack,args,promote_vexity([a.vexity for a in args]...),promote_sign([a.sign for a in args]...),sizes...)

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
  return dot(x, ones(x.size))
end
