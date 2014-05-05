import Base.sum, Base.abs, Base.sqrt, Base.log
export dot, transpose, ctranspose, sum, hcat, vcat
# TODO
# * max, min

# this breaks the syntax [x,y] to concatenate lists of expressions as in arguments to atoms
# extend to *args
# vcat(args::Array{AbstractCvxExpr}) = CvxExpr(:vstack,args,promote_vexity([a.vexity for a in args]...),promote_sign([a.sign for a in args]...),sizes...)

function hcat(args::Array{AbstractCvxExpr})
  num_rows = args[1].size[1]
  num_cols = 0
  # TODO: Vexity
  for arg in args
    if arg.size[1] != num_rows
      error("Must have same number of rows for horizontal concatentation")
    end
    num_cols += arg.size[2]
  end
  this = CvxExpr(:hcat, args, args[1].vexity, :any, (num_rows, num_cols))

  this.canon_form = ()->Any[]

  canon_constr_array = Any[]
  cols_so_far = 0
  for arg in args
    coeffs = spzeros(num_cols, arg.size[2])
    coeffs[1 + cols_so_far : arg.size[2] + cols_so_far, :] = speye(arg.size[2])
    append!(canon_constr_array, (this * coeffs == arg).canon_form())
    cols_so_far += arg.size[2]
  end
  this.canon_form = ()->canon_constr_array
end

function vcat(args::Array{AbstractCvxExpr})
  num_cols = args[1].size[2]
  num_rows = 0

  for arg in args
    if arg.size[1] != num_cols
      error("Must have same number of columns for vertical concatentation")
    end
    num_rows += arg.size[2]
  end
  this = CvxExpr(:vcat, args, args[1].vexity, :any, (num_rows, num_cols))

  this.canon_form = ()->Any[]
  canon_constr_array = Any[]
  rows_so_far = 0
  for arg in args
    coeffs = spzeros(arg.size[1], num_rows)
    coeffs[:, 1 + rows_so_far : arg.size[1] + rows_so_far] = speye(arg.size[1])
    append!(canon_constr_array, (coeffs * this == arg).canon_form())
    rows_so_far += arg.size[1]
  end
  this.canon_form = ()->canon_constr_array
end

function transpose(x::AbstractCvxExpr)
  sz = get_vectorized_size(x.size)
  coeffs = spzeros(sz, sz)
  num_rows = x.size[1]
  num_cols = x.size[2]
  for r = 1:num_rows
    for c = 1:num_cols
      i = (c - 1) * num_rows + r
      j = (r - 1) * num_cols + c
      coeffs[i, j] = 1.0
    end
  end

  this = CvxExpr(:transpose, [x], x.vexity, x.sign, (x.size[2], x.size[1]))
  canon_constr_array = Any[{
    :coeffs => Any[speye(sz), -coeffs],
    :vars => [x.uid(), this.uid()],
    :constant => spzeros(sz, 1),
    :is_eq => true
  }]
  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array
  return this
end

ctranspose(x::AbstractCvxExpr) = transpose(x)

function transpose(x::Constant)
  return Constant(x.value')
end

ctranspose(x::Constant) = transpose(x)

function dot(x::AbstractCvxExpr, y::AbstractCvxExpr)
  return x' * y
end

dot(x::Value, y::AbstractCvxExpr) = dot(convert(CvxExpr, x), y)
dot(x::AbstractCvxExpr, y::Value) = dot(x, convert(CvxExpr, y))

function sum(x::AbstractCvxExpr)
  if x.size == (1, 1)
    return x
  elseif x.size[1] == 1
    return x * ones(x.size[2], 1)
  else
    return ones(1, x.size[1]) * x
  end
end

function sum(x::Constant)
  return Constant(sum(x.value))
end

hcat
