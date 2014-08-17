#############################################################################
# variable.jl
# Defines Variable, which is a subtype of AbstractExpr
#############################################################################

export Variable
export vexity, evaluate, sign, dual_conic_form

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
    id_to_variables[this.id] = this
    return this
  end

  Variable(m::Integer, n::Integer, sign::Sign=NoSign()) = Variable((m,n), sign)
  Variable(sign::Sign=NoSign()) = Variable((1, 1), sign)
  Variable(size::Integer, sign::Sign=NoSign()) = Variable((size, 1), sign)
end

# GLOBAL MAP
# TODO: Comment David.
id_to_variables = Dict{Uint64, Variable}()

function vexity(x::Variable)
  return x.vexity
end

function evaluate(x::Variable)
  return x.value == nothing ? error("Value of the variable is yet to be calculated") : x.value
end

function sign(x::Variable)
  return x.sign
end

function dual_conic_form(x::Variable)
  var_to_coeff = Dict{Uint64, Value}()
  var_to_coeff[x.id] = speye(get_vectorized_size(x))
  # TODO add constraints for Variable sign when needed
  return (ConicObj(var_to_coeff), ConicConstr[])
end
