import ECOS
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
