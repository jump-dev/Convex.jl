#############################################################################
# constant.jl
# Defines Constant, which is a subtype of AbstractExpr
#############################################################################
export Constant
export vexity, evaluate, sign, conic_form

type Constant <: AbstractExpr
  head::Symbol
  value::Value
  size::(Int64, Int64)
  vexity::Vexity
  sign::Sign

  function Constant(x::Value, sign::Sign)
    sz = (size(x, 1), size(x, 2))
    return new(:constant, x, sz, ConstVexity(), sign)
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

function sign(x::Constant)
  return x.sign
end

function conic_form(x::Constant, unique_constr)
  id = object_id(x.value)
  if !((x.head, id) in keys(unique_constr))
    objective = ConicObj()
    objective[object_id(:constant)] = vec([x.value])
    unique_constr[(x.head, id)] = (objective, ConicConstr[])
  end
  return safe_copy(unique_constr[(x.head, id)])
end
