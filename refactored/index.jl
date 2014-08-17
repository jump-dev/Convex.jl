import Base.getindex
export IndexAtom, getindex

ArrayOrNothing = Union(AbstractArray, Nothing)

type IndexAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)
  rows::ArrayOrNothing
  cols::ArrayOrNothing
  inds::ArrayOrNothing

  function IndexAtom(x::AbstractExpr, rows::AbstractArray, cols::AbstractArray)
    sz = (length(rows), length(cols))
    children = (x,)
    return new(:index, hash(children), children, sz, rows, cols, nothing)
  end

  function IndexAtom(x::AbstractExpr, inds::AbstractArray)
    sz = (length(inds), 1)
    children = (x,)
    return new(:index, hash(children), children, sz, nothing, nothing, inds)
  end
end

function sign(x::IndexAtom)
  return sign(x.children[1])
end

function monotonicity(x::IndexAtom)
  return (Nondecreasing(),)
end

function curvature(x::IndexAtom)
  return ConstVexity()
end

function evaluate(x::IndexAtom)
  if x.inds == nothing
    return getindex(evaluate(x.children[1]), rows, cols)
  else
    return getindex(evaluate(x.children[1]), inds)
  end
end

function dual_conic_form(x::IndexAtom)
  index_matrix = spzeros(get_vectorized_size(x), get_vectorized_size(x.children[1]))
  if x.inds == nothing
    index_matrix_row = 1
    num_rows = x.children[1].size[1]
    for c in x.cols
      for r in x.rows
        index_matrix_col = num_rows * (convert(Int64, c) - 1) + convert(Int64, r)
        index_matrix[index_matrix_row, index_matrix_col] = 1
        index_matrix_row += 1
      end
    end
  else
    index_matrix_row = 1
    for i in x.inds
      index_matrix[index_matrix_row, i] = 1
      index_matrix_row += 1
    end
  end
  objective, constraints = dual_conic_form(x.children[1])
  objective = index_matrix * objective
  return (objective, constraints)
end

getindex{T <: Real}(x::AbstractExpr, rows::AbstractArray{T, 1}, cols::AbstractArray{T, 1}) = IndexAtom(x, rows, cols)
getindex{T <: Real}(x::AbstractExpr, inds::AbstractArray{T, 1}) = IndexAtom(x, inds)
getindex(x::AbstractExpr, ind::Real) = getindex(x, ind:ind)
getindex(x::AbstractExpr, row::Real, col::Real) = getindex(x, row:row, col:col)
getindex{T <: Real}(x::AbstractExpr, row::Real, cols::AbstractArray{T, 1}) = getindex(x, row:row, cols)
getindex{T <: Real}(x::AbstractExpr, rows::AbstractArray{T, 1}, col::Real) = getindex(x, rows, col:col)
