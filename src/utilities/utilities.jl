import Base.vec
export convert
export get_vectorized_size, kron, reverse_vexity, reverse_sign
export unique_id, is_affine, is_convex, is_concave

function is_affine(vexity::Symbol)
  return vexity == :constant || vexity == :affine
end

function is_convex(vexity::Symbol)
  return vexity == :convex || is_affine(vexity)
end

function is_concave(vexity::Symbol)
  return vexity == :concave || is_affine(vexity)
end

function convert(::Type{CvxExpr}, x)
  if typeof(x) == CvxExpr
    return x
  else
    return Constant(x)
  end
end

# Julia cannot vectorize sparse matrices. This will handle it for now
function vec(x::SparseMatrixCSC)
  return Base.reshape(x, size(x, 1) * size(x, 2), 1)
end

function vec(x::Number)
  return [x]
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

# TODO: This is taken from the julia code, remove after updating to new version
function kron{Tv1,Ti1,Tv2,Ti2}(A::SparseMatrixCSC{Tv1,Ti1}, B::SparseMatrixCSC{Tv2,Ti2})
  Tv_res = promote_type(Tv1, Tv2)
  Ti_res = promote_type(Ti1, Ti2)
  A = convert(SparseMatrixCSC{Tv_res,Ti_res}, A)
  B = convert(SparseMatrixCSC{Tv_res,Ti_res}, B)
  return Base.kron(A, B)
end

kron(A::SparseMatrixCSC, B::VecOrMat) = kron(A, sparse(B))
kron(A::VecOrMat, B::SparseMatrixCSC) = kron(sparse(A), B)

kron(A::SparseMatrixCSC, B::Number) = kron(A, [B])
kron(A::Number, B::SparseMatrixCSC) = kron([A], B)

# Unique ids
unique_id(x::AbstractCvxExpr) = convert(Int64, object_id(x))
unique_id(x::CanonicalConstr) = convert(Int64, object_id(x))
