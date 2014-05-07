import Base.vec
export convert, promote_for_add, print_debug
export get_vectorized_size, full, kron, reverse_vexity, reverse_sign

### Conversion and promotion
# TODO: The difference between conversion and promotion is messy.
function convert(::Type{CvxExpr}, x)
  if typeof(x) == CvxExpr
    return x
  else
    return Constant(x)
  end
end

# Julia cannot vectorize sparse matrices. This will handle it for now
function vec(x::SparseMatrixCSC)
  return Base.vec(full(x))
end

### Utility functions for arithmetic

function get_vectorized_size(sz::(Int64, Int64))
  return sz[1] * sz[2]
end

function get_vectorized_size(x::AbstractCvxExpr)
  return x.size[1] * x.size[2]
end

function reverse_vexity(x::AbstractCvxExpr)
  vexity = x.vexity
  if vexity == :convex
    return :concave
  elseif vexity == :concave
    return :convex
  else
    return vexity
  end
end

function reverse_sign(x::AbstractCvxExpr)
  sign = x.sign
  if sign == :neg
    return :pos
  elseif sign == :pos
    return :neg
  else
    return sign
  end
end

function promote_for_add(x::Constant, sz::(Int64, Int64))
  this = Constant(x.value * ones(sz...), x.sign)
  return this
end

function promote_for_add(x::AbstractCvxExpr, sz::(Int64, Int64))
  this = Constant(ones(sz...), :pos) * x
  return this
end

function promote_for_add(x::AbstractCvxExpr, y::AbstractCvxExpr)
  if x.size != y.size
    if maximum(x.size) == 1
      x = promote_for_add(x, y.size)
    elseif maximum(y.size) == 1
      y = promote_for_add(y, x.size)
    else
      error("size of arguments cannot be added; got $(x.size),$(y.size)")
    end
  end

  return (x, y)
end

function print_debug(debug, args...)
  if (debug)
    println(args)
  end
end


# TODO: This is taken from the julia code, remove after updating to new version
function kron{Tv1,Ti1,Tv2,Ti2}(A::SparseMatrixCSC{Tv1,Ti1}, B::SparseMatrixCSC{Tv2,Ti2})
  Tv_res = promote_type(Tv1, Tv2)
  Ti_res = promote_type(Ti1, Ti2)
  A = convert(SparseMatrixCSC{Tv_res,Ti_res}, A)
  B = convert(SparseMatrixCSC{Tv_res,Ti_res}, B)
  return Base.kron(A,B)
end

kron(A::SparseMatrixCSC, B::VecOrMat) = kron(A, sparse(B))
kron(A::VecOrMat, B::SparseMatrixCSC) = kron(sparse(A), B)

kron(A::SparseMatrixCSC, B::Number) = kron(A, [B])
kron(A::Number, B::SparseMatrixCSC) = kron([A], B)
