import Base.size, Base.endof, Base.ndims
export AbstractExpr, AbstractFunc, ConvexFunc, ConcaveFunc, AffineFunc, NoVexityFunc
export Constant, Variable
export vexity, sign, size, evaluate
export endof, ndims
export VecOrMatOrSparse, VecOrMatOrSparseOrNothing, Value, ArrayFloat64OrNothing, ValueOrNothing


### Abstract types

abstract AbstractExpr
abstract AbstractFunc <: AbstractExpr
abstract AffineFunc <: AbstractFunc
abstract ConvexFunc <: AbstractFunc
abstract ConcaveFunc <: AbstractFunc
abstract NoVexityFunc <: AbstractFunc


function vexity(x::AbstractExpr)
  return x.vexity
end

function vexity(x::AbstractFunc)
  monotonicities = monotonicity(x)
  vexity = intrinsic_vexity(x)
  for i = 1:length(x.children)
    vexity += monotonicities[i] * vexity(x.children[i])
  end
  return vexity
end

function sign(x::AbstractExpr)
  return x.sign
end

function size(x::AbstractExpr)
  return x.size
end



### User-defined Unions

VecOrMatOrSparse = Union(Vector, Matrix, SparseMatrixCSC)
VecOrMatOrSparseOrNothing = Union(Vector, Matrix, SparseMatrixCSC, Nothing)
ArrayFloat64OrNothing = Union(Array{Float64, }, Nothing)
Value = Union(Number, AbstractArray)
ValueOrNothing = Union(Value, Nothing)


### Constant Type

type Constant <: AbstractExpr
  head::Symbol
  value::Value
  vexity::Vexity
  sign::Sign
  size::(Int64, Int64)

  function Constant(x::Value, sign::Sign=NoSign())
    sz = (size(x, 1), size(x, 2))
    return new(:constant, x, Constant(), sign, sz)
  end

  function Constant(x::Value, check_sign::Bool=true)
    if check_sign
      if all(x .>= 0)
        return Constant(x, Positive())
      elseif all(x .<= 0)
        return Constant(x, Negative())
      end
    end
    return Constant(x, NoSign())
  end
end

function evaluate(x::Constant)
  return x.value
end
### Variable Type

type Variable <: AbstractExpr
  head::Symbol
  value::ValueOrNothing
  vexity::Vexity
  sign::Sign
  size::(Int64, Int64)

  function Variable(size::(Int64, Int64), sign::Sign=NoSign())
    if length(size) == 1
      size = (size[1], 1)
    end
    if !(sign in signs)
      error("Sign must be one of :pos, :neg, or :any; got $sign")
    end
    return new(:variable, nothing, Affine(), sign, size)
  end

  Variable(m::Integer, n::Integer, sign::Sign=NoSign()) = Variable((m,n), sign)
  Variable(sign::Sign=NoSign()) = Variable((1, 1), sign)
  Variable(size::Integer, sign::Sign=NoSign()) = Variable((size, 1), sign)
end

function evaluate(x::Variable)
  return x.value == nothing ? error("Value of the variable is yet to be calculated") : x.value
end


### Indexing Utilities

endof(x::AbstractCvxExpr) = x.size[1] * x.size[2]

function size(x::AbstractCvxExpr, dim::Integer)
  if dim < 1
    error("dimension out of range")
  elseif dim > 2
    return 1
  else
    return size(x)[dim]
  end
end

ndims(x::AbstractCvxExpr) = 2
