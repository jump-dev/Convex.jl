import Base.convert, Base.size
export AbstractCvxExpr, CvxExpr, Variable, Parameter, Constant, size, unique_id, Value

abstract AbstractCvxExpr
# every type inheriting from the AbstractCvxExpr type should have the following properties
#   head --- a symbol
#   vexity --- one of :linear, :convex, :concave, :constant
#   sign   --- one of :pos, :neg, :any
#   size   --- a tuple giving the size of the expression
vexities = [:constant, :linear, :convex, :concave]
signs    = [:pos, :neg, :any, :zero] # xxx add zero, useful for parameters
# Values consist of any type that can be a Constant
Value = Union(Number,AbstractArray)
size(x::AbstractCvxExpr) = x.size
sign(x::AbstractCvxExpr) = x.sign
vexity(x::AbstractCvxExpr) = x.vexity

function eval(x::AbstractCvxExpr)
  if x.head == :sum

  end
end
# function display(x::CvxExpr)
#   println("$(x.head)($(display(y) for y in x.args))")
# end

type CvxExpr <: AbstractCvxExpr
  head::Symbol
  args::Array{AbstractCvxExpr}
  vexity::Symbol
  sign::Symbol
  size::Tuple
  uid::Function
  canon_form::Function
  # TODO: args::Array works, everything else fucks up (eg args or args::Array{AbstractCvxExpr})
  # Check why
  function CvxExpr(head::Symbol,args::Array,vexity::Symbol,sign::Symbol,size::Tuple)
    if !(sign in signs)
      error("sign must be one of :pos, :neg, :any; got $sign")
    elseif !(vexity in vexities)
      error("vexity must be one of :constant, :linear, :convex, :concave; got $vexity")
    else
      this = new(head,args,vexity,sign,size)
      this.uid = ()->unique_id(this)
      return this
    end
  end
end

CvxExpr(head::Symbol,arg,vexity::Symbol,sign::Symbol,size::Tuple) = CvxExpr(head,[arg],vexity,sign,size)

type Variable <: AbstractCvxExpr
  head::Symbol
  value
  vexity::Symbol
  sign::Symbol
  size::Tuple
  uid::Function
  canon_form::Function
  function Variable(head::Symbol,size::Tuple,sign::Symbol)
    if !(sign in signs)
      error("sign must be one of :pos, :neg, :zero, :any; got $sign")
    end
    if head == :variable
      this = new(head,nothing,:linear,sign,size)
    elseif head == :parameter
      this = new(head,nothing,:constant,sign,size)
    end
    this.uid = ()->unique_id(this)
    # Variables are already in canonical form
    this.canon_form = ()->Any[]
    return this
  end
end

Variable(size::Tuple,sign::Symbol; kwargs...) = Variable(:variable,size,sign; kwargs...)
Variable(size::Tuple; kwargs...) = Variable(size,:any; kwargs...)
Variable(size...; kwargs...) = Variable(size,:any; kwargs...)
Variable(size::Integer,sign::Symbol; kwargs...) = Variable(Tuple(size),sign; kwargs...)
Parameter(size::Tuple,sign::Symbol; kwargs...) = Variable(:parameter,size,sign; kwargs...)
Parameter(size::Tuple; kwargs...) = Parameter(size,:any; kwargs...)
Parameter(size...; kwargs...) = Parameter(size,:any; kwargs...)
Parameter(size::Integer,sign::Symbol; kwargs...) = Parameter(Tuple(size),sign; kwargs...)

function parameter!(x::Variable)
  x.head = :parameter
  x.vexity = :constant
  x
end

# TODO: Why did madeleine do this
function variable!(x::Variable)
  x.head = :variable
  x.vexity = :linear
  x
end

type Constant <: AbstractCvxExpr
  head::Symbol
  value::Value
  vexity::Symbol
  sign::Symbol
  size::Tuple
  canon_form::Function
  function Constant(x::Value,sign)
    if sign in signs
      sz = size(x) == () ? (1, 1) : (size(x, 1), size(x, 2))
      this = new(:constant,x,:constant,sign,sz)
      this.canon_form = ()->Any[]
      return this
    else
      error("sign must be one of :pos, :neg, :zero, :any; got $sign")
    end
  end
end

function Constant(x::Number)
  # find the sign for scalar constants
  if x > 0
    return Constant(x,:pos)
  elseif x < 0
    return Constant(x,:neg)
  elseif x == 0
    return Constant(x,:zero)
  end
end

# Case to catch arrays, not scalar numbers
Constant(x::Value) = Constant(x,:any)

### Unique ids
unique_id(x::AbstractCvxExpr) = ccall(:jl_symbol_name, Ptr{Uint8}, (Any,), x)
#convert(::Type{Bool}, x::CvxConstr) = unique_id(x.lhs) == unique_id(x.rhs)
