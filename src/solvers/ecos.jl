import ECOS
export solve

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
# An object of type Solution consisting of x, y, z, status and ret_val
# where x are the primal variables, y are the multipliers for the equality constraints
# z are the multipliers for the inequality constraints 
#
# Some of the keyword arguments such as n, m, c, h, G are required. Hence, their default
# values are nothing, forcing them to be provided. Other keyword arguments need not
# be provided.
#
function ecos_solve(;n::Int64=nothing, m::Int64=nothing, p::Int64=0, l::Int64=0,
    ncones::Int64=0, q::Array{Int64, }=[], G::VecOrMatOrSparse=nothing,
    A::VecOrMatOrSparseOrNothing=nothing, c::Array{Float64, }=nothing,
    h::Array{Float64, }=nothing, b::ArrayFloat64OrNothing=nothing, debug::Bool=false)
  ptr_work = ECOS.setup(n=n, m=m, p=p, l=l, ncones=ncones, q=q, G=G, c=c, h=h, A=A, b=b)
  ret_val = ECOS.solve(ptr_work)

  solution = get_ecos_solution(ptr_work, n, p, m, ret_val)

  # 4 means we keep x,y,s,z.
  num_vars_to_keep = 4
  ECOS.cleanup(ptr_work, 4)
  return solution
end
function solve(p::ECOSConicProblem)
  solution = ecos_solve(n=p.n, m=p.m, p=p.p, l=p.l, ncones=p.ncones, q=p.q, G=p.G, c=p.c, h=p.h, A=p.A, b=p.b)
  solution.optval = (p.c' * solution.x)[1] # Transpose returns an array, so fetch the element
  return solution
end
# Given the arguments, returns an object of type Solution
# x: primal variables
# y: dual variables for equality constraints
# s: slacks for Gx + s <= h, s \in K
# z: dual variables for inequality constraints s \in K
# note slacks are nonzero iff dual variables are zero, 
# by complementary slackness
function get_ecos_solution(ptr_work, n, p, m, ret_val)
  work = pointer_to_array(ptr_work, 1)[1]

  x = pointer_to_array(work.x, n)
  y = pointer_to_array(work.y, p)
  z = pointer_to_array(work.z, m)
  s = pointer_to_array(work.s, m)

  return Solution(x, y, z, ret_val)
end

function ecos_debug(problem::Problem)
  objective = problem.objective

  canonical_constraints_array = CanonicalConstr[]
  for constraint in problem.constraints
    append!(canonical_constraints_array, constraint.canon_form())
  end

  append!(canonical_constraints_array, objective.canon_form())
  return create_ecos_matrices(canonical_constraints_array)
end

# Given the canonical_constraints_array, creates conic inequality matrix G and h
# as well as the equality matrix A and b
function create_ecos_matrices(canonical_constraints_array)
  n = 0::Int64
  variable_index = Dict{Int64, Int64}()

  eq_constr_index = Dict{Int64, Int64}()
  ineq_constr_index = Dict{Int64, Int64}()

  m = 0::Int64
  p = 0::Int64
  l = 0::Int64
  ncones = 0::Int64
  q = Int64[]

  # Loop over all the constraints to figure out the size of G and A
  for constraint in canonical_constraints_array
    # Loop over each variable in the constraint
    length_constraint_vars = length(constraint.vars)
    for i = 1:length_constraint_vars
      var = constraint.vars[i]

      # If we haven't already taken into account the size of this variable,
      # add it to the size of the variable
      if !haskey(variable_index, var)
        variable_index[var] = n + 1

        n += size(constraint.coeffs[i], 2)
      end
    end

    if constraint.is_eq
      p += size(constraint.coeffs[1], 1)
    elseif constraint.is_conic
      ncones += 1
      push!(q, size(constraint.coeffs[1], 1))
      m += size(constraint.coeffs[1], 1)
    else
      l += size(constraint.coeffs[1], 1)
      m += size(constraint.coeffs[1], 1)
    end
  end

  h = m == 0 ? nothing: zeros(m, 1)
  G = m == 0 ? nothing: spzeros(m, n)
  b = p == 0 ? nothing: zeros(p, 1)
  A = p == 0 ? nothing: spzeros(p, n)

  l_index = 1::Int64
  c_index = 1::Int64
  p_index = 1::Int64

  # Now, we actually stuff the matrices A and G
  for constraint in canonical_constraints_array
    m_var = 0::Int64

    length_constraint_vars = length(constraint.vars)
    for i = 1:length_constraint_vars
      var = constraint.vars[i]
      # Technically, the m_var size of all the variables should be the same,
      # otherwise nothing makes sense
      m_var = size(constraint.coeffs[i], 1)
      n_var = size(constraint.coeffs[i], 2)

      if constraint.is_eq
        # TODO: Julia has problems not converting ints to floats
        # An issue has been filed and should be fixed in newer versions of julia
        A[p_index : p_index + m_var - 1, variable_index[var] : variable_index[var] + n_var - 1] =
          constraint.coeffs[i] * 1.0
      elseif constraint.is_conic
        G[l + c_index : l + c_index + m_var - 1, variable_index[var] : variable_index[var] + n_var - 1] =
          constraint.coeffs[i] * 1.0
      else
        G[l_index : l_index + m_var - 1, variable_index[var] : variable_index[var] + n_var - 1] =
          constraint.coeffs[i] * 1.0
      end
    end

    if constraint.is_eq
      eq_constr_index[constraint.uid] = p_index
      b[p_index : p_index + m_var - 1] = constraint.constant
      p_index += m_var
    elseif constraint.is_conic
      ineq_constr_index[constraint.uid] = l + c_index
      h[l + c_index : l + c_index + m_var - 1] = constraint.constant
      c_index += m_var
    else
      ineq_constr_index[constraint.uid] = l_index
      h[l_index : l_index + m_var - 1] = constraint.constant
      l_index += m_var
    end
  end

  return m, n, p, l, ncones, q, G, h, A, b, variable_index, eq_constr_index, ineq_constr_index
end