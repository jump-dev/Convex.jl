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

To disable this, run

```julia
Convex.emit_dcp_warnings() = false
```
to redefine the method. See [`Convex.emit_dcp_warnings`](@ref) for more details.


Dual Variables
--------------

Convex.jl also returns the optimal dual variables for a problem. These
are stored in the `dual` field associated with each constraint.

```julia
using Convex, SCS

x = Variable()
constraint = x >= 0
p = minimize(x, constraint)
solve!(p, SCS.Optimizer)

# Get the dual value for the constraint
p.constraints[1].dual
# or
constraint.dual
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

By making subtypes of [`Convex.AbstractVariable`](@ref) that conform to the appropriate
interface (see the [`Convex.AbstractVariable`](@ref) docstring for details), one can
easily provide custom variable types for specific constructions. These aren't
always necessary though; for example, one can define the following function
`probabilityvector`:

```@example prob
using Convex

function probabilityvector(d::Int)
    x = Variable(d, Positive())
    add_constraint!(x, sum(x) == 1)
    return x
end
```
and then use, say, `p = probabilityvector(3)` in any Convex.jl problem. The
constraints that the entries of `p` are non-negative and sum to 1 will be
automatically added to any problem `p` is used in.

Custom types are necessary when one wants to dispatch on custom variables, use
them as callable types, or provide a different implementation. Continuing with
the probability vector example, let's say we often use probability vectors
variables in taking expectation values, and we want to use function notation for
this. To do so, we define

```@example 1
using Convex
mutable struct ProbabilityVector <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    value::Convex.ValueOrNothing
    vexity::Convex.Vexity
    function ProbabilityVector(d)
        this = new(:ProbabilityVector, 0, (d,1), nothing, Convex.AffineVexity())
        this.id_hash = objectid(this)
        this
    end
end

Convex.get_constraints(p::ProbabilityVector) = [ sum(p) == 1 ]
Convex.sign(::ProbabilityVector) = Convex.Positive()
Convex.vartype(::ProbabilityVector) = Convex.ContVar

(p::ProbabilityVector)(x) = dot(p, x)
```

Then one can call `p = ProbabilityVector(3)` to construct a our custom variable
which can be used in Convex, which already encodes the appropriate constraints
(non-negative and sums to 1), and which can act on constants via `p(x)`. For
example,

```@example 1
using SCS
p = ProbabilityVector(3)
x = [1.0, 2.0, 3.0]
prob = minimize( p(x) )
solve!(prob, SCS.Optimizer)
evaluate(p) # [1.0, 0.0, 0.0]
```

Subtypes of `AbstractVariable` must have the fields `head`, `id_hash`, and
`size`, and `id_hash` must be populated as shown in the example. Then they must also

* either have a field `value`, or implement [`Convex._value`](@ref) and
  [`Convex.set_value!`](@ref)
* either have a field `vexity`, or implement [`Convex.vexity`](@ref) and
  [`Convex.vexity!`](@ref) (though the latter is only necessary if you wish to
  support [`Convex.fix!`](@ref) and [`Convex.free!`](@ref)
* have a field `constraints` or implement [`Convex.constraints`](@ref) (optionally,
  implement [`Convex.add_constraint!`](@ref) to be able to add constraints to your
  variable after its creation),
* either have a field `sign` or implement [`Convex.sign`](@ref), and
* either have a field `vartype`, or implement [`Convex.vartype`](@ref) (optionally,
  implement [`Convex.vartype!`](@ref) to be able to change a variables' `vartype`
  after construction.)


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

```@example trees
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

```@example trees
for (i, node) in enumerate(AbstractTrees.PostOrderDFS(p))
    println("Here's node $i via PostOrderDFS: $(summary(node))")
end
```

### PreOrderDFS

Iterator to visit the nodes of a tree, guaranteeing that parents
will be visited before their children.

```@example trees
for (i, node) in enumerate(AbstractTrees.PreOrderDFS(p))
    println("Here's node $i via PreOrderDFS: $(summary(node))")
end
```

### StatelessBFS

Iterator to visit the nodes of a tree, guaranteeing that all nodes of a level
will be visited before their children.

```@example trees
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
