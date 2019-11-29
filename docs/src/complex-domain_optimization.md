Optimization with Complex Variables
===================================

Convex.jl also supports optimization with complex variables. Below, we
present a quick start guide on how to use Convex.jl for optimization
with complex variables, and then list the operations supported on
complex variables in Convex.jl. In general, any operation available in
Convex.jl that is well defined and DCP compliant on complex variables
should be available. We list these functions below. organized by the
type of cone (linear, second-order, or semidefinite) used to represent
that operation.

Internally, Convex.jl transforms the complex-domain problem to a larger
real-domain problem using a bijective mapping. It then solves the
real-domain problem and transforms the solution back to the complex
domain.

Complex Variables
-----------------

Complex Variables in Convex.jl are declared in the same way as the
variables are declared but using the different keyword
`ComplexVariable`.

```julia
  # Scalar complex variable
  z = ComplexVariable()

  # Column vector variable
  z = ComplexVariable(5)

  # Matrix variable
  z = ComplexVariable(4, 6)

  # Complex Positive Semidefinite variable
  z = HermitianSemidefinite(4)
```

Linear Program Representable Functions (complex variables)
----------------------------------------------------------

All of the linear functions that are listed under
[Linear Program Representable Functions](@ref) operate
on complex variables as well. In addition, several specialized functions
for complex variables are available:

| operation           | description                             | vexity | slope      | notes                        |
| ------------------- | --------------------------------------- | ------ | ---------- | ---------------------------- |
| `real(z)`           | real part of complex of variable   | affine | increasing | none                         |
| `imag(z)`           | imaginary part of complex variable | affine | increasing | none                         |
| `conj(x)`           | element-wise complex conjugate          | affine | increasing | none                         |
| `innerproduct(x,y)` | `real(trace(x'*y))`                     | affine | increasing | PR: one argument is constant |

Second-Order Cone Representable Functions (complex variables)
-------------------------------------------------------------

Most of the second order cone function listed under
[Second-Order Cone Representable Functions](@ref)
operate on complex variables as well. Notable exceptions include:

> -   inverse
> -   square
> -   quadoverlin
> -   sqrt
> -   geomean
> -   huber

One new function is available:

| operation | description      | vexity | slope      | notes |
| --------- | ---------------- | ------ | ---------- | ----- |
| `abs2(z)` | `square(abs(z))` | convex | increasing | none  |

Semidefinite Program Representable Functions (complex variables)
----------------------------------------------------------------

All SDP-representable functions listed under
[Semidefinite Program Representable Functions](@ref)
work for complex variables.

Exponential + SDP representable Functions (complex variables)
-------------------------------------------------------------

Complex variables also support logdet function.

Optimizing over quantum states
------------------------------

The complex and Hermitian matrix variables, along with the
`kron` and `partialtrace` operations, enable the
definition of a wide range of problems in quantum information theory. As
a simple example, let us consider a state $\rho$ over a composite
Hilbert space $\mathcal{H}_A\otimes\mathcal{H}_B$, where both component
spaces are isomorphic to $\mathbb{C}^2$. Assume that $\rho$ is a product
state, with its component in $\mathcal{H}_A$ given as $A$, a
complex-valued matrix. We can optimize over the second component $B$ to
meet some requirement. Here we simply fix the second component too, but
via the `partialtrace` operator:

```@example
using Convex, SCS
A = [ 0.47213595 0.11469794+0.48586827im; 0.11469794-0.48586827im  0.52786405]
B = ComplexVariable(2, 2)
ρ = kron(A, B)
constraints = [partialtrace(ρ, 1, [2; 2]) == [1 0; 0 0]
               tr(ρ) == 1
               ρ in :SDP]
p = satisfy(constraints)
solve!(p, SCSSolver(verbose=false))
p.status
```

Since we fix both components as trace-1 positive semidefinite matrices,
the last two constraints are actually redundant in this case.
