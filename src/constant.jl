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

  function Constant(x::Value, sign::Sign)
    sz = (size(x, 1), size(x, 2))
    return new(:constant, object_id(x), x, sz, ConstVexity(), sign)
  end

  function Constant(x::Value, check_sign::Bool=true)
    if check_sign
      if is_complex(x)
        return Constant(x, ComplexSign())
      elseif all(x .>= 0)
        return Constant(x, Positive())
      elseif all(x .<= 0)
        return Constant(x, Negative())
      end
    end
    return Constant(x, NoSign())
  end

  function is_complex(x::Value)
    temp = false
    y = x + 0im
    for z in y
      if z.im != 0
        temp = true
        break
      end
    end
    return temp    
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


function real_conic_form(x::Constant)
  return real(x)
end

function imag_conic_form(x::Constant)
  imag(x)
end


  
function conic_form!(x::Constant, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    real_Value = real_conic_form(x)
    imag_Value = imag_conic_form(x) 
    objective = ConicObj()
    objective[object_id(:constant)] = sign(x)==ComplexSign()?(vec(real_Value),vec(imag_Value)):(vec(real_Value),)
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end
