#############################################################################
# variable.jl
# Defines Variable, which is a subtype of AbstractExpr
#############################################################################

export Variable, Semidefinite
export vexity, evaluate, sign, conic_form

type Variable <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  value::ValueOrNothing
  size::(Int64, Int64)
  vexity::Vexity
  sign::Sign
  implied_constraints::Array{Constraint, 1}

  # is_symmetric is only needed for Semidefinite atoms. Value is ignored for everything else
  # If you wish to force symmetricity for other variables, add x == x' as a constraint
  function Variable(size::(Int64, Int64), sign::Sign=NoSign(); is_symmetric=true)
    this = new(:variable, 0, nothing, size, AffineVexity(), sign, Constraint[])
    this.id_hash = object_id(this)
    id_to_variables[this.id_hash] = this
    if !(sign == NoSign())
      if sign == Positive()
        push!(this.implied_constraints, this >= 0)
      elseif sign == Negative()
        push!(this.implied_constraints, this <= 0)
      elseif sign == Semidefinite()
        push!(this.implied_constraints, SDPConstraint(this, is_symmetric=is_symmetric))
      end
    end
    return this
  end

  Variable(m::Integer, n::Integer, sign::Sign=NoSign()) = Variable((m,n), sign)
  Variable(sign::Sign=NoSign()) = Variable((1, 1), sign)
  Variable(size::Integer, sign::Sign=NoSign()) = Variable((size, 1), sign)
end

# convenience semidefinite matrix constructor
Semidefinite(m::Integer; is_symmetric=true) = Variable((m,m), Semidefinite(), is_symmetric=is_symmetric)
Semidefinite(m::Integer, n::Integer; is_symmetric=true) = begin
  m==n ? Variable((m,m), Semidefinite(); is_symmetric=is_symmetric) : error("Semidefinite matrices must be square")
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

function conic_form(x::Variable, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    objective = ConicObj()
    vec_size = get_vectorized_size(x)
    objective[x.id_hash] = speye(vec_size)
    objective[object_id(:constant)] = spzeros(vec_size, 1)
    # placeholder values in unique constraints prevent infinite recursion depth
    add_conic_form!(unique_conic_forms, x, objective)
    for constraint in x.implied_constraints
      conic_form(constraint, unique_conic_forms)
    end
  end
  return get_conic_form(unique_conic_forms, x)
end
