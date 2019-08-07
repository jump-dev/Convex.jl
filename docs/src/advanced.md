Advanced Features
=================

DCP warnings
------------

When an expression is created which is not of [DCP
form](https://dcp.stanford.edu/), a warning is emitted. For example,

```repl
x = Variable()
y = Variable()
x*y
```

To disable this, set the module-level parameter `DCP_WARNINGS` via

```julia
Convex.DCP_WARNINGS[] = false
```


Dual Variables
--------------

Convex.jl also returns the optimal dual variables for a problem. These
are stored in the `dual` field associated with each constraint.

```julia
using Convex, SCS

x = Variable()
constraint = x >= 0
p = minimize(x, constraint)
solve!(p, SCS.Optimizer())

# Get the dual value for the constraint
p.constraints[1].dual
# or
constraint.dual
```

Warmstarting
------------

If you're solving the same problem many times with different values of
a parameter, Convex.jl can initialize many solvers with the solution to
the previous problem, which sometimes speeds up the solution time. This
is called a **warm start**.

To use this feature, pass the optional argument
`warmstart=true` to the `solve!` method.

```julia
# initialize data
n = 1000
y = rand(n)
x = Variable(n)

# first solve
lambda = 100
problem = minimize(sumsquares(y - x) + lambda * sumsquares(x - 10))
@time solve!(problem, SCS.Optimizer)

# now warmstart
# if the solver takes advantage of warmstarts, 
# this run will be faster
lambda = 105
@time solve!(problem, SCS.Optimizer, warmstart=true)
```

Fixing and freeing variables
----------------------------

Convex.jl allows you to fix a variable `x` to a value by calling the `fix!`
method. Fixing the variable essentially turns it into a constant. Fixed
variables are sometimes also called parameters.

`fix!(x, v)` fixes the variable `x` to the value `v`.

`fix!(x)` fixes `x` to its current value, which might be the value obtained by
solving another problem involving the variable `x`.

To allow the variable `x` to vary again, call `free!(x)`.

Fixing and freeing variables can be particularly useful as a tool for performing
alternating minimization on nonconvex problems. For example, we can find an
approximate solution to a nonnegative matrix factorization problem with
alternating minimization as follows. We use warmstarts to speed up the solution.

```julia
# initialize nonconvex problem
n, k = 10, 1
A = rand(n, k) * rand(k, n)
x = Variable(n, k)
y = Variable(k, n)
problem = minimize(sum_squares(A - x*y), x>=0, y>=0)

# initialize value of y
set_value!(y, rand(k, n))
# we'll do 10 iterations of alternating minimization
for i=1:10 
    # first solve for x
    # with y fixed, the problem is convex
    fix!(y)
    solve!(problem, SCS.Optimizer, warmstart = i > 1 ? true : false)
    free!(y)

    # now solve for y with x fixed at the previous solution
    fix!(x)
    solve!(problem, SCS.Optimizer, warmstart = true)
    free!(x)
end
```


Custom Variable Types
---------------------

By making subtypes of [`AbstractVariable`](@ref) that conform to the appropriate
interface (see the [`AbstractVariable`](@ref) docstring for details), one can easily
provide custom variable types for specific constructions. For example, in
quantum mechanics, a density matrix is a complex positive semi-definite matrix
with unit trace. With the following code, we can define a `DensityMatrix` type
which encodes this information:
 
```@example 1
using Convex

mutable struct DensityMatrix <: Convex.AbstractVariable{ComplexF64}
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    value::Convex.ValueOrNothing
    vexity::Convex.Vexity
    function DensityMatrix(d)
        this = new(:DensityMatrix, 0, (d,d), nothing, Convex.AffineVexity())
        this.id_hash = objectid(this)
        this
    end
end
Convex.constraints(ρ::DensityMatrix) = [ ρ ⪰ 0, tr(ρ) == 1 ]
Convex.sign(::DensityMatrix) = Convex.ComplexSign()
Convex.vartype(::DensityMatrix) = Convex.ContVar
```

Then one can call `ρ = DensityMatrix(d)` to construct a variable which can be
used in Convex, and which already encodes the appropriate constraints (i.e.
positive semi-definite and unit trace). For example,

```@example 1
using Convex, LinearAlgebra, SCS, Test

X = randn(4,4) + im*rand(4,4); X = X+X'
# `X` is Hermitian and non-degenerate (with probability 1)
# Let us calculate the projection onto the eigenspace associated
# to the maximum eigenvalue

e_vals, e_vecs = LinearAlgebra.eigen(LinearAlgebra.Hermitian(X))
e_val, idx = findmax(e_vals)
e_vec = e_vecs[:, idx]
proj = e_vec * e_vec'

# found it!
# Now let us do it again via an SDP
ρ = DensityMatrix(4)

prob = maximize( real(tr(ρ*X)) )
solve!(prob, SCS.Optimizer(verbose=false))

@test prob.optval ≈ e_val atol = 1e-3
@test evaluate(ρ) ≈ proj atol = 1e-3
```


Printing and the tree structure
-------------------------------

A Convex problem is structured as a *tree*, with the *root* being the problem
object, with branches to the objective and the set of constraints. The objective
is an `AbstractExpr` which itself is a tree, with each atom being a node and
having `children` which are other atoms, variables, or constants. Convex
provides `children` methods from
[AbstractTrees.jl](https://github.com/Keno/AbstractTrees.jl) so that the
tree-traversal functions of that package can be used with Convex.jl problems and
structures. This is what allows powers the printing of problems, expressions,
and constraints. The depth to which the tree corresponding to a problem,
expression, or constraint is printed is controlled by the global variable
[`Convex.MAXDEPTH`](@ref), which defaults to 3. This can be changed by e.g.
setting

```julia
Convex.MAXDEPTH[] = 5
```

Likewise, [`Convex.MAXWIDTH`](@ref), which defaults to 15, controls the "width"
of the printed tree. For example, when printing a problem with 20 constraints,
only the first `MAXWIDTH` of the constraints will be printed. Vertical dots,
"⋮", will be printed indicating that some constraints were omitted in the
printing.

A related setting is [`Convex.MAXDIGITS`](@ref), which controls
printing the internal IDs of atoms: if the string representation of an
ID is longer than double the value of `MAXDIGITS`, then it is
shortened by printing only the first and last `MAXDIGITS` characters.

The AbstractTrees methods can also be used to analyze the structure
of a Convex.jl problem. For example,

```@repl 1
using Convex, AbstractTrees
x = Variable()
p = maximize( log(x), x >= 1, x <= 3 )
for leaf in AbstractTrees.Leaves(p)
    println("Here's a leaf: $(summary(leaf))")
end
```

We can also iterate over the problem in various orders. The following descriptions
are taken from the AbstractTrees.jl docstrings, which have more information.

### PostOrderDFS

Iterator to visit the nodes of a tree, guaranteeing that children
will be visited before their parents.

```@repl 1
for (i, node) in enumerate(AbstractTrees.PostOrderDFS(p))
    println("Here's node $i via PostOrderDFS: $(summary(node))")
end
```

### PreOrderDFS

Iterator to visit the nodes of a tree, guaranteeing that parents
will be visited before their children.

```@repl 1
for (i, node) in enumerate(AbstractTrees.PreOrderDFS(p))
    println("Here's node $i via PreOrderDFS: $(summary(node))")
end
```

### StatelessBFS

Iterator to visit the nodes of a tree, guaranteeing that all nodes of a level
will be visited before their children.

```@repl 1
for (i, node) in enumerate(AbstractTrees.StatelessBFS(p))
    println("Here's node $i via StatelessBFS: $(summary(node))")
end
```

Reference
---------

```@docs
Convex.MAXDEPTH
Convex.MAXWIDTH
Convex.MAXDIGITS
```
