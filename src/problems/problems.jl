import ECOS, MathProgBase

export Solution, Problem, minimize, maximize, satisfy

Float64OrNothing = Union(Float64, Nothing)

# Declares the Solution type, which stores the primal and dual variables as well
# as the status of the solver
# TODO: Call x, y and z primal, dual_equality and dual_inequality
# x: primal variables
# y: dual variables for equality constraints
# z: dual variables for inequality constraints s \in K
type Solution{T<:Number}
  x::Array{T, 1} # x: primal variables
  y::Array{T, 1} # y: dual variables for equality constraints
  z::Array{T, 1} # z: dual variables for inequality constraints s \in K
  status::Symbol       # status is a termination status symbol, one of :Optimal, :Infeasible, :Unbounded, :UserLimit (iteration limit or timeout), :Error (and maybe others).
  optval::T
  dual::Bool           # true if dual values (y and z) have been populated
end

Solution{T}(x::Array{T, 1}, status::Symbol, optval::T) = Solution(x, T[], T[], status, optval, false)
Solution{T}(x::Array{T, 1}, y::Array{T, 1}, z::Array{T, 1}, status::Symbol, optval::T) = Solution(x, y, z, status, optval, true)

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
  status::Symbol
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