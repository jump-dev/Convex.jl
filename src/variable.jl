#############################################################################
# variable.jl
# Defines Variable, which is a subtype of AbstractExpr
#############################################################################

export Variable, Semidefinite
export vexity, evaluate, sign, conic_form

type Variable <: AbstractExpr
  head::Symbol
  id::Uint64
  value::ValueOrNothing
  size::(Int64, Int64)
  vexity::Vexity
  sign::Sign
  implied_constraints::Array{Constraint, 1}

  function Variable(size::(Int64, Int64), sign::Sign=NoSign())
    this = new(:variable, 0, nothing, size, AffineVexity(), sign, Constraint[])
    this.id = object_id(this)
    id_to_variables[this.id] = this
    if !(sign == NoSign())
      if sign == Positive()
        push!(this.implied_constraints, this >= 0)
      elseif sign == Negative()
        push!(this.implied_constraints, this <= 0)
      elseif sign == Semidefinite()
        push!(this.implied_constraints, SDPConstraint(this))
      end
    end
    return this
  end

  Variable(m::Integer, n::Integer, sign::Sign=NoSign()) = Variable((m,n), sign)
  Variable(sign::Sign=NoSign()) = Variable((1, 1), sign)
  Variable(size::Integer, sign::Sign=NoSign()) = Variable((size, 1), sign)
end

# convenience semidefinite matrix constructor
Semidefinite(m::Integer) = Variable((m,m), Semidefinite())
Semidefinite(m::Integer, n::Integer) = (m==n ? Variable((m,m), Semidefinite()) : error("Semidefinite matrices must be square"))

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

function conic_form(x::Variable, unique_constr)
  if !((x.head, x.id) in keys(unique_constr))
    objective = ConicObj()
    vec_size = get_vectorized_size(x)
    objective[x.id] = speye(vec_size)
    objective[object_id(:constant)] = spzeros(vec_size, 1)
    # placeholder values in unique constraints prevent infinite recursion depth
    unique_constr[(x.head, x.id)] = (objective, ConicConstr[])
    _, constraints = conic_form(x.implied_constraints, unique_constr)
    return safe_copy((objective, constraints))
  end
  return safe_copy(unique_constr[(x.head, x.id)])
end
