export Solution

Float64OrNothing = Union(Float64, Nothing)

# Declares the Solution type, which stores the primal and dual variables as well
# as the status of the solver
# TODO: Call x, y and z primal, dual_equality and dual_inequality
# x: primal variables
# y: dual variables for equality constraints
# z: dual variables for inequality constraints s \in K
type Solution
  x::Array{Float64, 1} # x: primal variables
  y::Array{Float64, 1} # y: dual variables for equality constraints
  z::Array{Float64, 1} # z: dual variables for inequality constraints s \in K
  status::ASCIIString
  ret_val::Int64
  optval::Float64OrNothing

  const status_map = {
    0 => "solved",
    1 => "primal infeasible",
    2 => "dual infeasible",
    -1 => "max iterations reached",
    -2 => "numerical problems in solver",
    -3 => "numerical problems in solver"
  }

  function Solution(x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, ret_val::Int64, optval::Float64OrNothing=nothing)
    if haskey(status_map, ret_val)
      return new(x, y, z, status_map[ret_val], ret_val, optval)
    else
      return new(x, y, z, "unknown problem in solver", ret_val, optval)
    end
  end
end
