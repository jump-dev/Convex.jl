import Base.convert, Base.size
export AbstractCvxExpr, CvxExpr, Variable, Parameter, Constant, unique_id, Value

abstract AbstractCvxExpr
# Every type inheriting from the AbstractCvxExpr type should have the following properties:
#   head --- a symbol
#   vexity --- one of :linear, :convex, :concave, :constant
#   sign   --- one of :pos, :neg, :any
#   size   --- a tuple giving the size of the expression
const vexities = [:constant, :linear, :convex, :concave]
const signs = [:pos, :neg, :any, :zero]
# Values consist of any type that can be a Constant
Value = Union(Number,AbstractArray)

type CvxExpr <: AbstractCvxExpr
  head::Symbol
  args::Array{AbstractCvxExpr}
  vexity::Symbol
  sign::Symbol
  size::Tuple
  uid::Ptr{Uint8}
  canon_form::Function
  # TODO: args::Array works, everything else does not (eg args or args::Array{AbstractCvxExpr})
  # Check why
  function CvxExpr(head::Symbol, args::Array, vexity::Symbol, sign::Symbol, size::Tuple)
    if !(sign in signs)
      error("sign must be one of :pos, :neg, :any; got $sign")
    elseif !(vexity in vexities)
      error("vexity must be one of :constant, :linear, :convex, :concave; got $vexity")
    else
      this = new(head, args, vexity, sign, size)
      this.uid = unique_id(this)
      return this
    end
  end
end

CvxExpr(head::Symbol, arg, vexity::Symbol, sign::Symbol, size::Tuple) =
  CvxExpr(head, [arg], vexity, sign, size)

type Variable <: AbstractCvxExpr
  head::Symbol
  value
  vexity::Symbol
  sign::Symbol
  size::Tuple
  uid::Ptr{Uint8}
  canon_form::Function

  function Variable(head::Symbol, size::Tuple, sign::Symbol)
    if length(size) == 1
      size = (size[1], 1)
    end
    if !(sign in signs)
      error("sign must be one of :pos, :neg, :zero, :any; got $sign")
    end
    if head == :variable
      this = new(head, nothing, :linear, sign, size)
    elseif head == :parameter
      this = new(head, nothing, :constant, sign, size)
    end
    this.uid = unique_id(this)
    # Variables are already in canonical form
    this.canon_form = ()->Any[]
    return this
  end
end

Variable(size::Tuple, sign::Symbol) = Variable(:variable, size, sign)
Variable(size::Tuple) = Variable(size, :any)
Variable(size...) = Variable(size, :any)
Variable(size::Integer) = Variable((size, 1),:any)
Variable(size::Integer, sign::Symbol) = Variable((size, 1), sign)

Parameter(size::Tuple, sign::Symbol) = Variable(:parameter, size, sign)
Parameter(size::Tuple) = Parameter(size, :any)
Parameter(size...) = Parameter(size, :any)
Parameter(size::Integer, sign::Symbol) = Parameter(tuple(size), sign)

type Constant <: AbstractCvxExpr
  head::Symbol
  value::Value
  vexity::Symbol
  sign::Symbol
  size::Tuple
  canon_form::Function

  function Constant(x::Value,sign)
    if sign in signs
      sz = (size(x, 1), size(x, 2))
      # TODO: We're doing a double tranpose right now because the (1, ) causes
      # a bug in the sparse matrix implementation:
      # g = spzeros(5)
      # g[1:5, 1:1]=a
      # causes an error
      # Once julia fixes it, we can probably move back to x
      this = new(:constant,x'',:constant,sign,sz)
      this.canon_form = ()->Any[]
      return this
    else
      error("sign must be one of :pos, :neg, :zero or :any but got $sign")
    end
  end
end

function Constant(x::Number)
  # find the sign for scalar constants
  if x > 0
    return Constant(x, :pos)
  elseif x < 0
    return Constant(x, :neg)
  elseif x == 0
    return Constant(x, :zero)
  end
end

# Case to catch arrays, not scalar numbers
Constant(x::Value) = Constant(x, :any)

# Unique ids
unique_id(x::AbstractCvxExpr) = ccall(:jl_symbol_name, Ptr{Uint8}, (Any, ), x)
