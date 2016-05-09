#############################################################################
# constant.jl
# Defines Constant, which is a subtype of AbstractExpr
#############################################################################
export Constant
export vexity, evaluate, sign, domain, conic_form!

type Constant <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  value::Value
  size::Tuple{Int, Int}
  vexity::Vexity
  sign::Sign
  domain::Domain

  function Constant(x::Value, sign::Sign, domain::Domain)
    sz = (size(x, 1), size(x, 2))
    return new(:constant, object_id(x), x, sz, ConstVexity(), sign, domain)
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
#### Constant Definition end   ##### 

function vexity(x::Constant)
  return x.vexity
end

function evaluate(x::Constant)
  return x.value
end

function sign(x::Constant)
  return x.sign
end

# New Code
# A function to return domain
function domain(x::Constant)
  return x.domain
end

function conic_form!(x::Constant, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    objective = ConicObj()
    objective[object_id(:constant)] = vec([x.value;])
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end
