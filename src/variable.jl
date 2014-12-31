#############################################################################
# variable.jl
# Defines Variable, which is a subtype of AbstractExpr
#############################################################################

export Variable, Semidefinite
export vexity, evaluate, sign, conic_form!

type Variable <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  value::ValueOrNothing
  size::(Int, Int)
  vexity::Vexity
  sign::Sign
  sets # ::Array{Symbol,1}

  # is_symmetric is only needed for Semidefinite atoms. Value is ignored for everything else
  # If you wish to force symmetricity for other variables, add x == x' as a constraint
  function Variable(size::(Int, Int), sign::Sign=NoSign(), sets::Symbol...)
    this = new(:variable, 0, nothing, size, AffineVexity(), sign, sets)
    this.id_hash = object_id(this)
    id_to_variables[this.id_hash] = this
    return this
  end

  Variable(m::Int, n::Int, sign::Sign=NoSign(), sets::Symbol...) = Variable((m,n), sign, sets...)
  Variable(sign::Sign, sets::Symbol...) = Variable((1, 1), sign, sets...)
  Variable(sets::Symbol...) = Variable((1, 1), NoSign(), sets...)
  Variable(size::(Int, Int), sets::Symbol...) = Variable(size::(Int, Int), NoSign(), sets...)
  Variable(size::Int, sign::Sign=NoSign(), sets::Symbol...) = Variable((size, 1), sign, sets...)
  Variable(size::Int, sets::Symbol...) = Variable((size, 1), sets...)
end

# convenience semidefinite matrix constructor
Semidefinite(m::Int; is_symmetric=true) = Variable((m,m), is_symmetric ? :Semidefinite : :AsymSemidefinite)
function Semidefinite(m::Int, n::Int; is_symmetric=true)
  if m==n
    return Variable((m,m), is_symmetric ? :Semidefinite : :AsymSemidefinite)
  else 
    error("Semidefinite matrices must be square")
  end
end

# global map from unique variable ids to variables. 
# the expression tree will only utilize variable ids during construction
# full information of the variables will be needed during stuffing
# and after solving to populate the variables with values
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

function conic_form!(x::Variable, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    objective = ConicObj()
    vec_size = get_vectorized_size(x)
    objective[x.id_hash] = speye(vec_size)
    objective[object_id(:constant)] = spzeros(vec_size, 1)
    # placeholder values in unique constraints prevent infinite recursion depth
    cache_conic_form!(unique_conic_forms, x, objective)
    if !(x.sign == NoSign())
      conic_form!(x.sign, x, unique_conic_forms)
    end
    for set in x.sets
      conic_form!(set, x, unique_conic_forms)
    end
  end
  return get_conic_form(unique_conic_forms, x)
end
