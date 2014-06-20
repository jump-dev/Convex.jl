import Base.hcat, Base.vcat
export hcat, vertcat

AbstractCvxExprOrValue = Union(AbstractCvxExpr, Value)
# TODO: We should be able to call hcat and vcat as [..., ...] and [...; ...]
# Right now, this conflicts with the fact that Problem type has a field
# `constraints::Array{CvxConstr}`

# TODO: In general, revisit this file and see if there is a better way of doing this

# Horizontally concatenate the list of arguments
# Suppose the result of hcat will be matrix `w` with size m x n
# Then, we loop over each argument and add a canonical constraint of the form
# `w * coeffs == arg`
# where `coeffs` is appropriately constructed
function hcat(given_args::AbstractCvxExprOrValue...)

  # Convert Values to Constants
  args = AbstractCvxExpr[]
  for arg in given_args
    if typeof(arg) <: Value
      push!(args, Constant(arg))
    else
      push!(args, arg)
    end
  end

  num_rows = args[1].size[1]
  num_cols = 0
  vexity = args[1].vexity
  sign = args[1].sign
  for arg in args
    if arg.size[1] != num_rows
      error("Must have same number of rows for horizontal concatentation")
    end
    # Same promotion rules as addition
    vexity = promote_vexity_add(vexity, arg.vexity)
    sign = promote_sign_add(sign, arg.sign)
    num_cols += arg.size[2]
  end
  this = CvxExpr(:hcat, [args...], vexity, sign, (num_rows, num_cols))

  this.canon_form = ()->CanonicalConstr[]

  canon_constr_array = CanonicalConstr[]
  cols_so_far = 0
  # Loop over each argument and add a constraint of the form this * coeffs == arg
  for arg in args
    # Reinitialize coeffs to be zeros
    coeffs = spzeros(num_cols, arg.size[2])
    # Appropriately populate coeffs
    coeffs[1 + cols_so_far : arg.size[2] + cols_so_far, :] = speye(arg.size[2])
    append!(canon_constr_array, (this * coeffs == arg).canon_form())
    cols_so_far += arg.size[2]
  end
  this.canon_form = ()->canon_constr_array

  this.evaluate = ()->begin
    result = args[1].evaluate()
    for arg in args[2 : end]
      result = hcat(result, arg.evaluate())
    end
    return result
  end

  return this
end

# Vertically concatenate the arguments
# Similar logic to hcat
function vertcat(given_args::AbstractCvxExpr...)
  # Convert Values to Constants
  args = AbstractCvxExpr[]
  for arg in given_args
    if typeof(arg) <: Value
      push!(args, Constant(arg))
    else
      push!(args, arg)
    end
  end

  num_cols = args[1].size[2]
  num_rows = 0
  vexity = args[1].vexity
  sign = args[1].sign

  for arg in args
    if arg.size[2] != num_cols
      error("Must have same number of columns for vertical concatentation")
    end
    # Same promotion rules as addition
    vexity = promote_vexity_add(vexity, arg.vexity)
    sign = promote_sign_add(sign, arg.sign)
    num_rows += arg.size[1]
  end

  this = CvxExpr(:vcat, [args...], vexity, sign, (num_rows, num_cols))

  this.canon_form = ()->CanonicalConstr[]
  canon_constr_array = CanonicalConstr[]
  rows_so_far = 0
  for arg in args
    # Reinitialize coeffs to be zeros
    coeffs = spzeros(arg.size[1], num_rows)
    # Appropriately populate coeffs
    coeffs[:, 1 + rows_so_far : arg.size[1] + rows_so_far] = speye(arg.size[1])
    append!(canon_constr_array, (coeffs * this == arg).canon_form())
    rows_so_far += arg.size[1]
  end
  this.canon_form = ()->canon_constr_array

  this.evaluate = ()->begin
    result = args[1].evaluate()
    for arg in args[2 : end]
      result = vcat(result, arg.evaluate())
    end
    return result
  end

  return this
end

# TODO: Implement hvcat
# TODO: Implement repmat
