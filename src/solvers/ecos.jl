export ecos_solve

# Calls the ECOS C solver
#
# Input
# n: is the number of variables
# m: is the number of inequality constraints (dimension 1 of the matrix G and the
# length of the vector h)
# p: is the number of equality constraints (can be 0)
# l: is the dimension of the positive orthant, i.e. in Gx+s=h, s in K, the first l
# elements of s are >=0
# ncones: is the number of second-order cones present in K
# q: is an array of integers of length ncones, where each element defines the dimension
# of the cone
# c is an array of type float of size n
# h is an array of type float of size m
# b is an array of type float of size p (can be nothing if no equalities are present)
#
# Returns:
# An object of type Solution consisting of x, y, z, status and ret_va
# where x are the primal variables, y are the multipliers for the equality constraints
# z are the multipliers for the conic inequalities
#
# Some of the keyword arguments such as n, m, c, h, G are required. Hence, their default
# values are nothing, forcing them to be provided. Other keyword arguments need not
# be provided.
#
function ecos_solve(;n::Int64=nothing, m::Int64=nothing, p::Int64=0, l::Int64=0,
    ncones::Int64=0, q::Array{Int64, }=[0], G::VecOrMatOrSparse=nothing,
    A::VecOrMatOrSparseOrNothing=nothing, c::Array{Float64, }=nothing,
    h::Array{Float64, }=nothing, b::ArrayFloat64OrNothing=nothing, debug::Bool=false)

  if l == 0
    l = m
    print_debug(debug, "Value of l=0, setting it to the same as m=$m");
  end

  if q == [0]
    q = convert(Ptr{Int64}, C_NULL)
  end

  if A == nothing
    Apr = convert(Ptr{Float64}, C_NULL)
    Ajc = convert(Ptr{Int64}, C_NULL)
    Air = convert(Ptr{Int64}, C_NULL)
  else
    sparseA = sparse(A)
    # Hack to make it a float, find a better way
    Apr = sparseA.nzval * 1.0
    # -1 since the C language is 0 indexed
    Ajc = sparseA.colptr - 1
    Air = sparseA.rowval - 1
  end

  sparseG = sparse(G)
  # Hack to make it a float, find a better way
  Gpr = sparseG.nzval * 1.0
  # -1 since the C  language is 0 indexed
  Gjc = sparseG.colptr - 1
  Gir = sparseG.rowval - 1

  if b == nothing
    b = convert(Ptr{Float64}, C_NULL)
  end

  # Call ECOS to setup the problem
  pwork = ccall((:ECOS_setup, "../ecos/ecos.so"), Ptr{Void},
      (Int64, Int64, Int64, Int64, Int64, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64},
      Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
      n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b)

  # Solve the problem
  ret_val = ccall((:ECOS_solve, "../ecos/ecos.so"), Int64, (Ptr{Void},), pwork)

  solution = get_ecos_solution(pwork, n, p, m, ret_val)

  # TODO: Check how many we need to keep
  # 3 means we keep x,y,z 2 means x,y and so on
  num_vars_to_keep = 3
  ccall((:ECOS_cleanup, "../ecos/ecos.so"), Void, (Ptr{Void}, Int64), pwork, num_vars_to_keep)
  return solution
end


# Given the arguments, returns an object of type Solution
function get_ecos_solution(pwork, n, p, m, ret_val)
  double_ptr = convert(Ptr{Ptr{Float64}}, pwork)

  # x is the 5th element in the struct
  x_ptr = unsafe_load(double_ptr, 5)
  x = pointer_to_array(x_ptr, n)

  y_ptr = unsafe_load(double_ptr, 6)
  y = pointer_to_array(y_ptr, p)

  z_ptr = unsafe_load(double_ptr, 7)
  z = pointer_to_array(z_ptr, m)

  return Solution(x, y, z, ret_val)
end
