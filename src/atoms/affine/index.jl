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
    return new(:index, hash((children, rows, cols, nothing)), children, sz, rows, cols, nothing)
  end

  function IndexAtom(x::AbstractExpr, inds::AbstractArray)
    sz = (length(inds), 1)
    children = (x,)
    return new(:index, hash((children, nothing, nothing, inds)), children, sz, nothing, nothing, inds)
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

function conic_form(x::IndexAtom, unique_constr)
  if !haskey(unique_constr, (x.head, x.children_hash))
    m = get_vectorized_size(x)
    n = get_vectorized_size(x.children[1])

    if x.inds == nothing
      sz = length(x.cols) * length(x.rows)
      J = Array(Int64, sz)
      k = 1

      num_rows = x.children[1].size[1]
      for c in x.cols
        for r in x.rows
          J[k] = num_rows * (convert(Int64, c) - 1) + convert(Int64, r)
          k += 1
        end
      end

      index_matrix = sparse(1:sz, J, 1.0, m, n)
    else
      index_matrix = sparse(1:length(x.inds), x.inds, 1.0, m, n)
    end
    objective, constraints = conic_form(x.children[1], unique_constr)
    objective = index_matrix * objective
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return unique_constr[(x.head, x.children_hash)]
end

getindex{T <: Real}(x::AbstractExpr, rows::AbstractArray{T, 1}, cols::AbstractArray{T, 1}) = IndexAtom(x, rows, cols)
getindex{T <: Real}(x::AbstractExpr, inds::AbstractArray{T, 1}) = IndexAtom(x, inds)
getindex(x::AbstractExpr, ind::Real) = getindex(x, ind:ind)
getindex(x::AbstractExpr, row::Real, col::Real) = getindex(x, row:row, col:col)
getindex{T <: Real}(x::AbstractExpr, row::Real, cols::AbstractArray{T, 1}) = getindex(x, row:row, cols)
getindex{T <: Real}(x::AbstractExpr, rows::AbstractArray{T, 1}, col::Real) = getindex(x, rows, col:col)
