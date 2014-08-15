import Base.sign, Base.size, Base.endof, Base.ndims
export AbstractExpr, Constant, Variable
export vexity, sign, size, evaluate
export dual_conic_form
export endof, ndims
export Value, ValueOrNothing
export get_vectorized_size


### Abstract type

abstract AbstractExpr

function vexity(x::AbstractExpr)
  monotonicities = monotonicity(x)
  vexity = intrinsic_vexity(x)
  for i = 1:length(x.children)
    vexity += monotonicities[i] * vexity(x.children[i])
  end
  return vexity
end

function sign(x::AbstractExpr) # Constant & Vars
  return x.sign
end

function size(x::AbstractExpr)
  return x.size
end

### User-defined Unions
Value = Union(Number, AbstractArray)
ValueOrNothing = Union(Value, Nothing)


### Constant Type

type Constant <: AbstractExpr
  head::Symbol
  id::Uint64
  value::Value
  size::(Int64, Int64)
  vexity::Vexity
  sign::Sign

  function Constant(x::Value, sign::Sign=NoSign())
    sz = (size(x, 1), size(x, 2))
    return new(:constant, object_id(x), x, sz, Constant(), sign)
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

function vexity(x::Constant)
  return x.vexity
end

function evaluate(x::Constant)
  return x.value
end

function dual_conic_form(x::Constant)
  var_to_coeff = Dict{Uint64, Value}()
  var_to_coeff[object_id(:Constant)] = vec(x.value)
  return (ConicObj(var_to_coeff), ConicConstr[])
end


### Variable Type

type Variable <: AbstractExpr
  head::Symbol
  id::Uint64
  value::ValueOrNothing
  size::(Int64, Int64)
  vexity::Vexity
  sign::Sign

  function Variable(size::(Int64, Int64), sign::Sign=NoSign())
    this = new(:variable, 0, nothing, size, Affine(), sign)
    this.id = object_id(this)
    return this
  end

  Variable(m::Integer, n::Integer, sign::Sign=NoSign()) = Variable((m,n), sign)
  Variable(sign::Sign=NoSign()) = Variable((1, 1), sign)
  Variable(size::Integer, sign::Sign=NoSign()) = Variable((size, 1), sign)
end

function vexity(x::Variable)
  return x.vexity
end

function evaluate(x::Variable)
  return x.value == nothing ? error("Value of the variable is yet to be calculated") : x.value
end

function dual_conic_form(x::Variable)
  var_to_coeff = Dict{Uint64, Value}()
  var_to_coeff[x.id] = speye(get_vectorized_size(x))
  # TODO add constraints for Variable sign when needed
  return (ConicObj(var_to_coeff), ConicConstr[])
end


### Indexing Utilities

endof(x::AbstractExpr) = x.size[1] * x.size[2]

function size(x::AbstractExpr, dim::Integer)
  if dim < 1
    error("dimension out of range")
  elseif dim > 2
    return 1
  else
    return size(x)[dim]
  end
end

ndims(x::AbstractExpr) = 2

get_vectorized_size(x::AbstractExpr) = reduce(*,size(x))
