#############################################################################
# conv.jl
# Handles convolution between a constant vector and an expression vector.
#############################################################################

import Base.conv
export conv

function conv(x::Value, y::AbstractExpr)
  if length(size(x)) != 1 || y.size[2] != 1
    error("convolution only supported between two vectors")
  end
  m = length(x)
  n = y.size[1]
  X = spzeros(m+n - 1, n)
  for i = 1:n
    X[i:m+i-1, i] = x
  end
  return X*y
end
