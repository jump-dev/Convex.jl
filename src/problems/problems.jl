import ECOS, MathProgBase

export Solution, Problem, minimize, maximize, satisfy

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



Float64OrNothing = Union(Float64, Nothing)
SolutionOrNothing = Union(Solution, Nothing)

# The Problem type consists of an objective and a set of a constraints.
# The objective specifies what should be maximized/minimized whereas the constraints
# specify the different constraints on the problem.
# status, optval and solution are populated after solve! is called.
type Problem
  head::Symbol
  objective::AbstractCvxExpr
  constraints::Array{CvxConstr}
  status::ASCIIString
  optval::Float64OrNothing
  solution::SolutionOrNothing

  function Problem(head::Symbol, objective::AbstractCvxExpr, constraints::Array{CvxConstr}=CvxConstr[])
    if !all([x <= 1 for x in objective.size])
      error("Only scalar optimization problems allowed, but size(objective) = $(objective.size)")
    end

    if head == :minimize && objective.vexity == :concave
      error("Cannot minimize a concave function")
    elseif head == :maximize && objective.vexity == :convex
      error("Cannot maximize a convex function")
    elseif head != :maximize && head != :minimize
      error("Problem.head must be one of :minimize or :maximize")
    end

    new(head, objective, constraints, "not yet solved", nothing, nothing)
  end
end

Problem(head::Symbol, objective::AbstractCvxExpr, constraints::CvxConstr...) =
  Problem(head, objective, [constraints...])

# Allow users to simply type minimize or maximize
minimize(objective::AbstractCvxExpr, constraints::CvxConstr...) =
  Problem(:minimize, objective, [constraints...])
minimize(objective::AbstractCvxExpr, constraints::Array{CvxConstr}=CvxConstr[]) =
  Problem(:minimize, objective, constraints)
minimize(objective::Value, constraints::CvxConstr...) =
  minimize(convert(CvxExpr, objective), constraints)
minimize(objective::Value, constraints::Array{CvxConstr}=CvxConstr[]) =
  minimize(convert(CvxExpr, objective), constraints)

maximize(objective::AbstractCvxExpr, constraints::CvxConstr...) =
  Problem(:maximize, objective, [constraints...])
maximize(objective::AbstractCvxExpr, constraints::Array{CvxConstr}=CvxConstr[]) =
  Problem(:maximize, objective, constraints)
maximize(objective::Value, constraints::CvxConstr...) =
  maximize(convert(CvxExpr, objective), constraints)
maximize(objective::Value, constraints::Array{CvxConstr}=CvxConstr[]) =
  maximize(convert(CvxExpr, objective), constraints)

satisfy(constraints::Array{CvxConstr}=CvxConstr[]) =
  Problem(:minimize, Constant(0), constraints)
satisfy(constraint::CvxConstr) = satisfy([constraint])

# +(constraints, constraints) is overwritten in constraints.jl
add_constraints(p::Problem, constraints::Array{CvxConstr}) = +(p.constraints, constraints)
add_constraints(p::Problem, constraint::CvxConstr) = add_constraints(p, [constraint])