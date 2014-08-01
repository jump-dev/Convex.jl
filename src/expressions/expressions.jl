import Base.convert, Base.size, Base.endof, Base.ndims
export AbstractCvxExpr, CvxExpr, Variable, Parameter, Constant, Value
export endof, size, ndims

abstract AbstractCvxExpr
# Every type inheriting from the AbstractCvxExpr type should have the following properties:
## head
## vexity - one of :affine, :convex, :concave, :constant
## sign - one of :pos, :neg, :any
## size - a tuple giving the size of the expression
## canon_form - a function that returns an array of type CanonicalConstr with the canonical form
## of itself and its descendants
## evaluate - a function that evaluates the expression and returns the result
const vexities = [:constant, :affine, :convex, :concave]
const signs = [:pos, :neg, :any, :zero]
# Values consist of any type that can be a Constant
Value = Union(Number, AbstractArray)

type CvxExpr <: AbstractCvxExpr
  head::Symbol
  args::Array{AbstractCvxExpr}
  vexity::Symbol
  sign::Symbol
  size::(Int64, Int64)
  uid::Int64
  canon_form::Function
  evaluate::Function
  # TODO: Stop using args::Array. Use a tuple instead
  function CvxExpr(head::Symbol, args::Array, vexity::Symbol, sign::Symbol, size::(Int64, Int64))
    if !(sign in signs)
      error("Sign must be one of :pos, :neg, :any; got $sign")
    elseif !(vexity in vexities)
      error("Vexity must be one of :constant, :affine, :convex, :concave; got $vexity")
    else
      this = new(head, args, vexity, sign, size)
      this.uid = unique_id(this)
      return this
    end
  end
end

CvxExpr(head::Symbol, arg, vexity::Symbol, sign::Symbol, size::(Int64, Int64)) =
  CvxExpr(head, [arg], vexity, sign, size)

type Variable <: AbstractCvxExpr
  head::Symbol
  vexity::Symbol
  sign::Symbol
  size::(Int64, Int64)
  value::ValueOrNothing
  uid::Int64
  canon_form::Function
  evaluate::Function

  function Variable(head::Symbol, size::(Int64, Int64), sign::Symbol)
    if length(size) == 1
      size = (size[1], 1)
    end
    if !(sign in signs)
      error("Sign must be one of :pos, :neg, :zero, :any; got $sign")
    end
    if head == :variable
      this = new(head, :affine, sign, size, nothing)
    elseif head == :parameter
      this = new(head, :constant, sign, size, nothing)
    end
    this.uid = unique_id(this)
    # Variables are already in canonical form
    this.canon_form = ()->CanonicalConstr[]
    this.evaluate = ()->this.value == nothing ? error("Value of the variable is yet to be calculated") : this.value
    return this
  end
end


Variable(size::(Int64, Int64), sign::Symbol) = Variable(:variable, size, sign)
Variable(size::(Int64, Int64)) = Variable(size, :any)
Variable(size...) = Variable(size, :any)
Variable() = Variable((1, 1), :any)
Variable(size::Integer) = Variable((size, 1),:any)
Variable(size::Integer, sign::Symbol) = Variable((size, 1), sign)

# TODO: Parameters are currently not in use
Parameter(size::(Int64, Int64), sign::Symbol) = Variable(:parameter, size, sign)
Parameter(size::(Int64, Int64)) = Parameter(size, :any)
Parameter(size...) = Parameter(size, :any)
Parameter(size::Integer, sign::Symbol) = Parameter((size, 1), sign)

type Constant <: AbstractCvxExpr
  head::Symbol
  value::Value
  vexity::Symbol
  sign::Symbol
  size::(Int64, Int64)
  canon_form::Function
  evaluate::Function

  function Constant(x::Value, sign)
    if sign in signs
      sz = (size(x, 1), size(x, 2))
      # TODO: We're doing a double tranpose right now because the (1, ) causes
      # a bug in the sparse matrix implementation:
      # g = spzeros(5)
      # g[1:5, 1:1]=a
      # causes an error
      # Once julia fixes it, we can probably move back to x
      this = new(:constant, x, :constant, sign, sz)
      this.canon_form = ()->CanonicalConstr[]
      this.evaluate = ()->this.value
      return this
    else
      error("Sign must be one of :pos, :neg, :zero or :any but got $sign")
    end
  end
end

function Constant(x::Number)
  # Find the sign for scalar constants
  if x > 0
    return Constant(x, :pos)
  elseif x < 0
    return Constant(x, :neg)
  elseif x == 0
    return Constant(x, :zero)
  end
end

function Constant(x::AbstractArray; check_sign::Bool=false)
  if check_sign
    if all(x .>= 0)
      return Constant(x, :pos)
    elseif all(x .<= 0)
      return Constant(x, :neg)
    end
  end
  return Constant(x, :any)
end

# The following functions are needed for indexing

endof(x::AbstractCvxExpr) = x.size[1] * x.size[2]

function size(x::AbstractCvxExpr, dim::Integer)
  if dim < 1
    error("dimension out of range")
  elseif dim > 2
    return 1
  else
    return x.size[dim]
  end
end

ndims(x::AbstractCvxExpr) = 2
