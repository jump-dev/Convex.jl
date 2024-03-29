# Contributing

We'd welcome contributions to the Convex.jl package. Here are some
short instructions on how to get started. If you don't know what you'd
like to contribute, you could

 -   take a look at the current
     [issues](https://github.com/jump-dev/Convex.jl/issues) and pick
     one. (Feature requests are probably the easiest to tackle.)
 -   add a [usage
     example](https://github.com/jump-dev/Convex.jl/tree/master/examples).

Then submit a pull request (PR). (Let us know if it's a work in
progress by putting \[WIP\] in the name of the PR.)

## Adding examples

 -   Take a look at our existing [usage
     examples](https://github.com/jump-dev/Convex.jl/tree/master/examples)
     and add another in similar style.
 -   Submit a PR. (Let us know if it's a work in progress by putting
     \[WIP\] in the name of the PR.)
 -   We'll look it over, fix up anything that doesn't work, and merge
     it.

## Adding atoms

Here are the steps to add a new function or operation (atom) to
Convex.jl. Let's say you're adding the new function $f$.

 -   Take a look at the
     [nuclear norm atom](https://github.com/jump-dev/Convex.jl/blob/master/src/atoms/sdp_cone/nuclearnorm.jl)
     for an example of how to construct atoms, and see the
     [norm atom](https://github.com/jump-dev/Convex.jl/blob/master/src/atoms/second_order_cone/norm.jl)
     for an example of an atom that depends on a parameter.
 -   Copy paste (for example, the nuclear norm file), replace anything saying
     nuclear norm with the name of the atom $f$, fill in monotonicity,
     curvature, etc. Save it in the appropriate subdirectory of
     `src/atoms/`.
 -   Ensure the atom is a mutable struct, so that `objectid` can be called.
 -   Add as a comment a description of what the atom does and its
     parameters.
 -   The most mathematically interesting part is the `new_conic_form!`
     function. Following the example in the nuclear norm atom, you'll
     see that you can just construct the problem whose optimal value is
     $f(x)$, introducing any auxiliary variables you need, exactly as
     you would normally in Convex.jl, and then call `conic_form!`
     on that problem.
 -   Add a test for the atom so we can verify it works in
     `src/problem_depot/problem/<cone>`, where `<cone>` matches the subdirectory of
     `src/atoms`. See [How to write a ProblemDepot problem](@ref) for details
     on how to write the tests.
 -   Following the other examples, add a test to `test/test_atoms.jl`.
 -   Submit a PR, including a description of what the atom does and its
     parameters. (Let us know if it's a work in progress by putting
     \[WIP\] in the name of the PR.)
 -   We'll look it over, fix up anything that doesn't work, and merge
     it.

## Fixing the guts

If you want to do a more major bug fix, you may need to understand how
Convex.jl thinks about conic form. To do this, start by reading [the
Convex.jl paper](http://arxiv.org/pdf/1410.4821.pdf). You may find our
[JuliaCon 2014 talk](https://www.youtube.com/watch?v=SoI0lEaUvTs&t=128s)
helpful as well; you can find the Jupyter notebook presented in the talk
[here](https://github.com/JuliaCon/presentations/tree/master/CVX).

Convex has been updated several times over the years however, so older information
may be out of date. Here is a brief summary of how the package works (as of July 2023).

1. A `Problem{T}` struct is created by putting together an objective function and constraints.
   This forms a tree of sorts, in which variables and constants are the leaves, and atoms form the
   intermediate branches. Here `T` refers to the numeric type of the problem. Variables and constants don't have such a type, and can
   be used in multiple problems with different types. We only have the final type when we have
   constructed the problem itself at the end.
2. When we go to `solve!` a problem, we first load it into a MathOptInterface (MOI) model.
   To do so, we need to traverse the problem and apply our extended formulations. This occurs
   via `conic_form!`. We construct a `Context{T}` associated to the problem, which holds an MOI
   model, and progressively load it by applying `conic_form!` to each object's children and then itself.
   For variables outputs of `conic_form!` are of types: `SparseTape{T}` or `ComplexTape{T}`,
   depending on the sign variable. Likewise for constant, the outputs of `conic_form!` are either `Vector{T}`
   or `ComplexStructOfVec{T}`. Here a `Tape` refers to a lazy sequence of sparse affine
   operators that will be applied to a vector of variables. The central computational task of Convex
   is to compose this sequence of operators (and thus enact it's extended formulations). For atoms,
   `conic_form!` generally either creates a new object using Convex' primitives (for example, another problem)
   and calls `conic_form!` on that, or, when that isn't possible, calls `operate` to
   manipulate the tape objects themselves (for example, to add a new operation to the composition).
   We try to minimize the amount of `operate` methods and defer to existing primitives when possible.
   `conic_form!` can also create new constraints and add them directly to the model. It is easy
   to create constraints of the form "vector-affine-function-in-cone" for any of MOI's many supported cones;
   these constraints do not need to be exposed at the level of Convex itself as `Constraint` objects, although they can be.
3. Once we have filled our `Context{T}`, we go to solve it with MOI. Then we
   recover the solution status and values of primal and dual variables, and
   populate them using dictionaries stored in the `Context`.

You're now armed and dangerous. Go ahead and open an issue (or comment
on a previous one) if you can't figure something out, or submit a PR if
you can figure it out. (Let us know if it's a work in progress by
putting \[WIP\] in the name of the PR.)

PRs that comment the code more thoroughly will also be welcomed.

## Developer notes

* `conic_form!` is allowed to mutate the context, but should never mutate the atoms or problems
* We currently construct a fresh context on every solve. It may be possible to set things up to reuse contexts for efficiency.
* Data flow: we take in user data that may be of any type.
    * At the level of problem formulation (when we construct atoms), we convert everything to an `AbstractExpr` (or `Constraint`); in particular, constants become `Constant` or `ComplexConstant`. At this time we don't know the numeric type that will be used to solve the problem.
    * Once we begin to `solve!` the problem, we recursively call `conic_form!`. The output is of type `SparseTape{T}`, `ComplexTape{T}`, `Vector{T}`, or `ComplexStructOfVec{T}`. We can call `operate` to manipulate these outputs.
    * We convert these to `MOI.VectorAffineFunction` before passing them to MOI.
