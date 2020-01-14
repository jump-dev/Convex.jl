Basic Types
===========

The basic building block of Convex.jl is called an *expression*, which
can represent a variable, a constant, or a function of another
expression. We discuss each kind of expression in turn.

Variables
---------

The simplest kind of expression in Convex.jl is a variable. Variables in
Convex.jl are declared using the `Variable` keyword, along
with the dimensions of the variable.

```julia
# Scalar variable
x = Variable()

# Column vector variable
x = Variable(5)

# Matrix variable
x = Variable(4, 6)
```

Variables may also be declared as having special properties, such as
being

-   (entrywise) positive: `x = Variable(4, Positive())`
-   (entrywise) negative: `x = Variable(4, Negative())`
-   integral: `x = Variable(4, IntVar)`
-   binary: `x = Variable(4, BinVar)`
-   (for a matrix) being symmetric, with nonnegative eigenvalues (ie,
     positive semidefinite): `z = Semidefinite(4)`

The order of the arguments is the size, the sign, and then the
[`Convex.VarType`](@ref) (i.e., integer, binary, or continuous), and any may be omitted
to use the default.
The current value of a variable `x` can be accessed with `evaluate(x)`. After
`solve!`ing a problem, the value of each variable used in the problem is set to
its optimal value.

See also [Custom Variable Types](@ref) for how to implement your own variable
types.

Constants
---------

Numbers, vectors, and matrices present in the Julia environment are
wrapped automatically into a `Constant` expression when used
in a Convex.jl expression.

Expressions
-----------

Expressions in Convex.jl are formed by applying any *atom* (mathematical
function defined in Convex.jl) to variables, constants, and other
expressions. For a list of these functions, see
[Operations](@ref). Atoms are applied to expressions using
operator overloading. For example, `2+2` calls Julia's built-in
addition operator, while `2+x` calls the Convex.jl addition method and
returns a Convex.jl expression. Many of the useful language features in
Julia, such as arithmetic, array indexing, and matrix transpose are
overloaded in Convex.jl so they may be used with variables and
expressions just as they are used with native Julia types.

Expressions that are created must be DCP-compliant. More information on
DCP can be found [here](http://dcp.stanford.edu/). :

```julia
x = Variable(5)
# The following are all expressions
y = sum(x)
z = 4 * x + y
z_1 = z[1]
```

Convex.jl allows the values of the expressions to be evaluated directly.

```julia
x = Variable()
y = Variable()
z = Variable()
expr = x + y + z
problem = minimize(expr, x >= 1, y >= x, 4 * z >= y)
solve!(problem, SCS.Optimizer)

# Once the problem is solved, we can call evaluate() on expr:
evaluate(expr)
```

Constraints
-----------

*Constraints* in Convex.jl are declared using the standard comparison
operators `<=`, `>=`, and `==`. They specify relations that must hold
between two expressions. Convex.jl does not distinguish between strict
and non-strict inequality constraints.

```julia
x = Variable(5, 5)
# Equality constraint
constraint = x == 0
# Inequality constraint
constraint = x >= 1
```

Matrices can also be constrained to be positive semidefinite.

```julia
x = Variable(3, 3)
y = Variable(3, 1)
z = Variable()
# constrain [x y; y' z] to be positive semidefinite
constraint = ([x y; y' z] in :SDP)
# or equivalently,
constraint = ([x y; y' z] ⪰ 0)
```

Constraints can also be added to variables after their construction, to automatically apply constraints
to any problem which uses the variable. For example,

```julia
x = Variable(3)
add_constraint!(x, sum(x) == 1)
```

Now, in any problem in which `x` is used, the constraint `sum(x) == 1` will be added.

Objective
---------

The objective of the problem is a scalar expression to be maximized or
minimized by using `maximize` or `minimize` respectively. Feasibility
problems can be expressed by either giving a constant as the objective,
or using `problem = satisfy(constraints)`.

Problem
-------

A *problem* in Convex.jl consists of a *sense* (minimize, maximize, or
satisfy), an *objective* (an expression to which the sense verb is to be
applied), and zero or more *constraints* that must be satisfied at the
solution. Problems may be constructed as

```julia
problem = minimize(objective, constraints)
# or
problem = maximize(objective, constraints)
# or
problem = satisfy(constraints)
```

Constraints can be added at any time before the problem is solved.

```julia
# No constraints given
problem = minimize(objective)
# Add some constraint
problem.constraints += constraint
# Add many more constraints
problem.constraints += [constraint1, constraint2, ...]
```

A problem can be solved by calling `solve!`

```julia
solve!(problem, solver)
```

passing a solver such as `SCS.Optimizer()` from the package `SCS` as the
second argument. After the problem is solved, `problem.status` records
the status returned by the optimization solver, and can be `:Optimal`,
`:Infeasible`, `:Unbounded`, `:Indeterminate` or `:Error`. If the status
is `:Optimal`, `problem.optval` will record the optimum value of the
problem. The optimal value for each variable `x` participating in the
problem can be found in `evaluate(x)`. The optimal value of an expression
can be found by calling the `evaluate()` function on the expression as
follows: `evaluate(expr)`.
