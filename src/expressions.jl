import Base.convert, Base.size
export AbstractCvxExpr, CvxExpr, Variable, Parameter, Constant, size

abstract AbstractCvxExpr
# every type inheriting from the AbstractCvxExpr type should have the following properties
# 	head --- a symbol
# 	vexity --- one of :linear, :convex, :concave, :constant
# 	sign   --- one of :pos, :neg, :any
# 	size   --- a tuple giving the size of the expression
global vexities = [:constant, :linear, :convex, :concave]
global signs    = [:pos, :neg, :any, :zero] # xxx add zero, useful for parameters
size(x::AbstractCvxExpr) = x.size
sign(x::AbstractCvxExpr) = x.sign
vexity(x::AbstractCvxExpr) = x.vexity
# function display(x::CvxExpr)
# 	println("$(x.head)($(display(y) for y in x.args))")
# end

type CvxExpr <: AbstractCvxExpr
	head
	args
	vexity
	sign
	size
	function CvxExpr(head::Symbol,args::Array,vexity::Symbol,sign::Symbol,size::Tuple)
		if !(sign in signs)
			error("sign must be one of :pos, :neg, :any; got $sign")
		elseif !(vexity in vexities)
			error("vexity must be one of :constant, :linear, :convex, :concave; got $vexity")
		else
			return new(head,args,vexity,sign,size)
		end
	end
end
CvxExpr(head::Symbol,arg,vexity::Symbol,sign::Symbol,size::Tuple) = CvxExpr(head,[arg],vexity,sign,size)

type Variable <: AbstractCvxExpr
	head
	value
	vexity
	sign
	size
	Variable(size::Tuple,sign::Symbol) = sign in signs ? new(:variable,nothing,:linear,sign,size) : error("sign must be one of :pos, :neg, :any; got $sign")
end
Variable(size::Tuple; kwargs...) = Variable(size,:any; kwargs...)
Variable(size...; kwargs...) = Variable(size,:any; kwargs...)
Variable(size::Integer,sign::Symbol; kwargs...) = Variable(Tuple(size),sign; kwargs...)

type Parameter <: AbstractCvxExpr
	head
	value
	vexity
	sign
	size
	Parameter(size::Tuple,sign::Symbol) = sign in signs ? new(:parameter,nothing,:constant,sign,size) : error("sign must be one of :pos, :neg, :any; got $sign")
end
Parameter(size::Tuple; kwargs...) = Parameter(size,:any; kwargs...)
Parameter(size...; kwargs...) = Parameter(size,:any; kwargs...)
Parameter(size::Integer,sign::Symbol; kwargs...) = Parameter(Tuple(size),sign; kwargs...)

type Constant <: AbstractCvxExpr
	head
	value
	vexity
	sign
	size
	Constant(x::Array,sign) = sign in signs ? new(:constant,x,:constant,sign,size(x)) : error("sign must be one of :pos, :neg, :any; got $sign")
	Constant(x::Number,sign) = sign in signs ? new(:constant,[x],:constant,sign,()) : error("sign must be one of :pos, :neg, :any; got $sign")
end
function Constant(x)
	# find the sign for scalar constants
	try
		if x > 0
			return Constant(x,:pos)
		elseif x < 0
			return Constant(x,:neg)	
		elseif x == 0
			return Constant(x,:zero)
		end
	catch
		return Constant(x,:any)
	end
end

### Unique ids
unique_id(x::AbstractCvxExpr) = ccall(:jl_symbol_name, Ptr{Uint8}, (Any,), x)
#convert(::Type{Bool}, x::CvxConstr) = unique_id(x.lhs) == unique_id(x.rhs)