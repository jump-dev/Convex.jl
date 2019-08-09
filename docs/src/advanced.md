Advanced Features
=================

Dual Variables
--------------

Convex.jl also returns the optimal dual variables for a problem. These
are stored in the `dual` field associated with each constraint.

```julia
using Convex, SCS

x = Variable()
constraint = x >= 0
p = minimize(x, constraint)
solve!(p, SCSSolver())

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
@time solve!(problem, SCSSolver())

# now warmstart
# if the solver takes advantage of warmstarts, 
# this run will be faster
lambda = 105
@time solve!(problem, SCSSolver(), warmstart=true)
```

Fixing and freeing variables
----------------------------

Convex.jl allows you to fix a variable `x` to a value by
calling the `fix!` method. Fixing the variable essentially
turns it into a constant. Fixed variables are sometimes also called
parameters.

`fix!(x, v)` fixes the variable `x` to the value
`v`.

`fix!(x)` fixes `x` to the value
`x.value`, which might be the value obtained by solving
another problem involving the variable `x`.

To allow the variable `x` to vary again, call
`free!(x)`.

Fixing and freeing variables can be particularly useful as a tool for
performing alternating minimization on nonconvex problems. For example,
we can find an approximate solution to a nonnegative matrix
factorization problem with alternating minimization as follows. We use
warmstarts to speed up the solution.

```julia
# initialize nonconvex problem
n, k = 10, 1
A = rand(n, k) * rand(k, n)
x = Variable(n, k)
y = Variable(k, n)
problem = minimize(sum_squares(A - x*y), x>=0, y>=0)

# initialize value of y
y.value = rand(k, n)
# we'll do 10 iterations of alternating minimization
for i=1:10 
    # first solve for x
    # with y fixed, the problem is convex
    fix!(y)
    solve!(problem, SCSSolver(), warmstart = i > 1 ? true : false)
    free!(y)

    # now solve for y with x fixed at the previous solution
    fix!(x)
    solve!(problem, SCSSolver(), warmstart = true)
    free!(x)
end
```

Custom Variable Types
---------------------

By making subtypes of `AbstractVariable` that conform to the appropriate interface (see the `AbstractVariable` docstring for details), one can easily provide custom variable types for specific constructions. For example, in quantum mechanics, a density matrix is a complex positive semi-definite matrix with unit trace. With the following code, we can define a `DensityMatrix` type which encodes this information:
 
```@example 1
using Convex

mutable struct DensityMatrix <: Convex.AbstractVariable{ComplexF64}
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    value::Convex.ValueOrNothing
    vexity::Vexity
    function DensityMatrix(d)
        this = new(:DensityMatrix, 0, (d,d), nothing, Convex.AffineVexity())
        this.id_hash = objectid(this)
        Convex.id_to_variables[this.id_hash] = this
        this
    end
end
Convex.constraints(ρ::DensityMatrix) = [ ρ ⪰ 0, tr(ρ) == 1 ]
Convex.sign(::DensityMatrix) = Convex.ComplexSign()
Convex.vartype(::DensityMatrix) = Convex.ContVar
```

Then one can call `ρ = DensityMatrix(d)` to construct a variable which can be used in Convex, and which already encodes the appropriate constraints (i.e. positive semi-definite and unit trace). For example,

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
solve!(prob, SCSSolver())

@test prob.optval ≈ e_val atol = 1e-3
@test evaluate(ρ) ≈ proj atol = 1e-3
```