import Base.reshape
export reshape

function reshape(x::Number, m::Int64, n::Int64)
  if m*n != 1
    error("Improper size for reshaping a number")
  end
  return x
end
