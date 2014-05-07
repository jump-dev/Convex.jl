export Solution

# Declares the Solution type, which stores the primal and dual variables as well
# as the status of the solver
type Solution
  x::Array{Float64, 1}
  y::Array{Float64, 1}
  z::Array{Float64, 1}
  status::ASCIIString
  ret_val::Int64

  const status_map = {
    0 => "solved",
    1 => "primal infeasible",
    2 => "dual infeasible",
    -1 => "max iterations reached",
    -2 => "numerical problems in solver",
    -3 => "numerical problems in solver"
  }

  function Solution(x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, ret_val::Int64)
    if ret_val in status_map
      return new(x, y, z, status_map[ret_val], ret_val)
    else
      return new(x, y, z, "unknown problem in solver", ret_val)
    end
  end
end
