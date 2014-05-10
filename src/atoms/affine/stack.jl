export hcat, vertcat

function hcat(args::AbstractCvxExpr...)
  num_rows = args[1].size[1]
  num_cols = 0
  # TODO: Vexity
  for arg in args
    if arg.size[1] != num_rows
      error("Must have same number of rows for horizontal concatentation")
    end
    num_cols += arg.size[2]
  end
  this = CvxExpr(:hcat, [args...], args[1].vexity, :any, (num_rows, num_cols))

  this.canon_form = ()->CanonicalConstr[]

  canon_constr_array = CanonicalConstr[]
  cols_so_far = 0
  for arg in args
    coeffs = spzeros(num_cols, arg.size[2])
    coeffs[1 + cols_so_far : arg.size[2] + cols_so_far, :] = speye(arg.size[2])
    append!(canon_constr_array, (this * coeffs == arg).canon_form())
    cols_so_far += arg.size[2]
  end
  this.canon_form = ()->canon_constr_array

  return this
end

# Called vertcat since vcat would interfere with how we pass the constraints
# to the problem
function vertcat(args::AbstractCvxExpr...)
  num_cols = args[1].size[2]
  num_rows = 0

  for arg in args
    if arg.size[2] != num_cols
      error("Must have same number of columns for vertical concatentation")
    end
    num_rows += arg.size[1]
  end
  # TODO: not doing vexity
  this = CvxExpr(:vcat, [args...], args[1].vexity, :any, (num_rows, num_cols))

  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = CanonicalConstr[]
  rows_so_far = 0
  for arg in args
    coeffs = spzeros(arg.size[1], num_rows)
    coeffs[:, 1 + rows_so_far : arg.size[1] + rows_so_far] = speye(arg.size[1])
    append!(canon_constr_array, (coeffs * this == arg).canon_form())
    rows_so_far += arg.size[1]
  end
  this.canon_form = ()->canon_constr_array

  return this
end
