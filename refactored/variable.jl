#############################################################################
# variable.jl
# Defines Variable, which is a subtype of AbstractExpr
#############################################################################

export Variable, Semidefinite
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

# convenience semidefinite matrix constructor
Semidefinite(m::Integer) = Variable((m,m), Semidefinite())
Semidefinite(m::Integer, n::Integer) = (m==n ? return Variable((m,m), Semidefinite()) : error("Semidefinite matrices must be square"))

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
  vec_size = get_vectorized_size(x)
  var_to_coeff[x.id] = speye(vec_size)
  constraints = sign_constraint(x.sign, var_to_coeff, vec_size)
  return (ConicObj(var_to_coeff), constraints)
end

function sign_constraint(s::NoSign, var_to_coeff, vec_size)
  return ConicConstr[]
end

function sign_constraint(s::Positive, var_to_coeff, vec_size)
  return [ConicConstr(var_to_coeff, :NonNeg, vec_size)]
end

function sign_constraint(s::Negative, var_to_coeff, vec_size)
  return [ConicConstr(-var_to_coeff, :NonNeg, vec_size)]
end

function sign_constraint(s::Semidefinite, var_to_coeff, vec_size)
  return [ConicConstr(var_to_coeff, :Semidefinite, vec_size)]
end