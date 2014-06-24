import Base.sum
export sum

# Sum all the elements in x
# Implemented as 1^T x 1, ie, by multiplying x with ones(...) on both sides
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

# Sum along a given dimension
function sum(x::AbstractCvxExpr, dimension::Int64)
  if dimension == 1
    return ones(1, x.size[1]) * x
  elseif dimension == 2
    return x * ones(x.size[2], 1)
  else
    error("Sum not implemented for dimension $dimension")
  end
end

function sum(x::Constant)
  return Constant(sum(x.value))
end

function sum(x::Constant, dimension::Int64)
  return Constant(sum(x.value, dimension))
end
