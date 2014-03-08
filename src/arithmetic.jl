export +,-,*,/

### Conversion and promotion
function convert(::Type{CvxExpr},x)
	if typeof(x) == CvxExpr
		return x
	else 
		return Constant(x)
	end
end
promote_rule(::Type{CvxExpr}, ::Type{AbstractArray}) = CvxExpr
promote_rule(::Type{CvxExpr}, ::Type{Number}) = CvxExpr

### Utility functions for arithmetic
function promote_size(x::AbstractCvxExpr,y::AbstractCvxExpr)
	if x.size == y.size
		size = x.size
	elseif length(x.size) == 0 
		size = y.size
	elseif length(y.size) == 0 
		size = x.size
	elseif length(x.size) == 1 && Set(x.size[1],1) == Set(y.size...)
		size = y.size
	elseif length(y.size) == 1 && Set(y.size[1],1) == Set(x.size...)
		size = x.size	
	else
		error("size of arguments must be the same; got $(x.size),$(y.size)")
	end
	return size
end

function promote_vexity(x::AbstractCvxExpr,y::AbstractCvxExpr)
	v1 = x.vexity; v2 = y.vexity; vexities = Set(v1,v2)
	if vexities == Set(:convex,:concave)
		error("expression not DCP compliant")
	elseif :convex in vexities
		return :convex
	elseif :concave in vexities
		return :concave
	elseif :linear in vexities
		return :linear
	else 
		return :constant
	end
end

function promote_sign(x::AbstractCvxExpr,y::AbstractCvxExpr)
	s1 = x.sign; s2 = y.sign; signs = Set(s1,s2)
	if :any in signs || signs == Set(:pos,:neg)
		return :any
	else # then s1==s2
		return s1
	end
end

# multiple arguments
for op = (:promote_vexity, :promote_sign, :promote_shape)
  @eval ($op)(arg1,arg2,arg3,args...) = length(args)==0 ? ($op)(($op)(arg1,arg2),arg3) : ($op)(($op)(arg1,arg2),arg3,args...)
end

function reverse_vexity(x::AbstractCvxExpr)
	vexity = x.vexity
	if vexity == :convex
		return :concave
	elseif vexity == :concave
		return :convex
	else
		return vexity
	end
end

function reverse_sign(x::AbstractCvxExpr)
	sign = x.sign
	if sign == :neg
		return :pos
	elseif sign == :pos
		return :neg
	else
		return sign
	end
end

### Unary functions on expressions
function -(x::AbstractCvxExpr)
	return CvxExpr(:-,[x],reverse_vexity(x),reverse_sign(x),size(x))
end
function size(x::AbstractCvxExpr)
	return x.size
end

### Arithmetic for expressions
+(x::AbstractCvxExpr,y::AbstractCvxExpr) = CvxExpr(:+,[x,y],promote_vexity(x,y),promote_sign(x,y),promote_size(x,y))
+(x::AbstractCvxExpr,y::Value) = +(x,convert(CvxExpr,y))
+(x,y::AbstractCvxExpr) = +(y,convert(CvxExpr,x))
-(x::AbstractCvxExpr,y::AbstractCvxExpr) = +(x,-y)
-(x::AbstractCvxExpr,y) = +(x,-y)
-(x,y::AbstractCvxExpr) = +(-y,x)

function *(x::AbstractCvxExpr,y::AbstractCvxExpr)
	# determine vexity
	if Set(x.vexity,y.vexity) == Set(:constant,:linear)
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

	# determine size
	# multiplication by a scalar
	if length(x.size) == 0
		size = y.size
	elseif length(y.size) == 0
		size = x.size
	# matrix multiplication
	elseif length(y.size) == 2 && length(x.size) == 2 && x.size[2] == y.size[1]
		size = (x.size[1],y.size[2])
	# cast x to a row vector
	elseif length(y.size) == 2 && length(x.size) == 1 && x.size[1] == y.size[1]
		size = (1,y.size[2])
	# cast x to a column vector
	elseif length(y.size) == 2 && length(x.size) == 1 && y.size[1] == 1
		size = (x.size[1],y.size[2])
	# cast y to a row vector
	elseif length(x.size) == 2 && length(y.size) == 1 && x.size[2] == 1
		size = (1,y.size[2])
	# cast y to a column vector
	elseif length(x.size) == 2 && length(y.size) == 1 && x.size[2] == y.size[1]
		size = (x.size[1],1)
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
	CvxExpr(:*,[x,y],vexity,sign,size)
end
*(x::AbstractCvxExpr,y::Value) = *(x,convert(CvxExpr,y))
*(x::Value,y::AbstractCvxExpr) = *(convert(CvxExpr,x),y)

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
