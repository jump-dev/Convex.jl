export ecos_solve

function print_debug(args...)
  println(args)
end

# TODO: Document
function ecos_solve(;n=nothing, m=nothing, p=nothing, l=nothing, ncones=nothing,
    q=nothing, G=nothing, A=nothing, c=nothing, h=nothing, b=nothing )

  @assert n != nothing
  @assert m != nothing
  @assert p != nothing

  @assert c != nothing
  @assert h != nothing

  @assert G != nothing

  if l == nothing
    l = m
    print_debug("Value of l=nothing, setting it to the same as m=$m");
  end

  if ncones == nothing
    ncones = 0
  end

  if q == nothing
    q = convert(Ptr{Int64}, C_NULL)
  end

  if A == nothing
    Apr = convert(Ptr{Float64}, C_NULL)
    Ajc = convert(Ptr{Int64}, C_NULL)
    Air = convert(Ptr{Int64}, C_NULL)
  else
    sparseA = sparse(A)
    # TODO: hack to make it float, find a better way
    Apr = sparseA.nzval * 1.0
    # -1 since C is 0 indexed
    Ajc = sparseA.colptr - 1
    Air = sparseA.rowval - 1
  end

  sparseG = sparse(G)
  # TODO: hack to make it float, find a better way
  Gpr = sparseG.nzval * 1.0
  # -1 since C is 0 indexed
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
  # 4 means we keep x,y,z,s, 3 means x,y,z and so on
  num_vars_to_keep = 4
  ccall((:ECOS_cleanup, "../ecos/ecos.so"), Void, (Ptr{Void}, Int64), pwork, num_vars_to_keep)
  return solution
end


# Given the arguments, returns a dictionary with variables x, y, z, s and status
function get_ecos_solution(pwork, n, p, m, ret_val)
  double_ptr = convert(Ptr{Ptr{Float64}}, pwork)
  # TODO: Worry about freeing memory?

  # x is the 12th
  x_ptr = unsafe_load(double_ptr, 12)
  x = pointer_to_array(x_ptr, n)

  y_ptr = unsafe_load(double_ptr, 13)
  y = pointer_to_array(y_ptr, p)

  z_ptr = unsafe_load(double_ptr, 14)
  z = pointer_to_array(z_ptr, m)

  s_ptr = unsafe_load(double_ptr, 15)
  s = pointer_to_array(s_ptr, m)

  if ret_val == 0
    status = "solved"
  elseif ret_val == 1
    status = "primal infeasible"
  else
    status = "dual infeasible"
  end

  return ["x"=> x, "y"=> y, "z"=>s, "s"=>s, "status"=> status, "ret_val"=>ret_val]
end
