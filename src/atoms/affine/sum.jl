import Base.sum
export sum

function sum(x::AbstractCvxExpr)
  if x.size == (1, 1)
    return x
  elseif x.size[1] == 1
    return x * ones(x.size[2], 1)
  elseif x.size[2] == 1
    return ones(1, x.size[1]) * x
  else
    return ones(1, x.size[1]) * x * ones(x.size[2], 1)
  end
end

function sum(x::Constant)
  return Constant(sum(x.value))
end
