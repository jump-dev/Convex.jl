Contributing
============

We'd welcome contributions to the Convex.jl package. Here are some
short instructions on how to get started. If you don't know what you'd
like to contribute, you could

 -   take a look at the current
     [issues](https://github.com/JuliaOpt/Convex.jl/issues) and pick
     one. (Feature requests are probably the easiest to tackle.)
 -   add a [usage
     example](https://github.com/JuliaOpt/Convex.jl/tree/master/examples).

Then submit a pull request (PR). (Let us know if it's a work in
progress by putting \[WIP\] in the name of the PR.)

Adding examples
---------------

 -   Take a look at our exising [usage
     examples](https://github.com/JuliaOpt/Convex.jl/tree/master/examples)
     and add another in similar style.
 -   Submit a PR. (Let us know if it's a work in progress by putting
     \[WIP\] in the name of the PR.)
 -   We'll look it over, fix up anything that doesn't work, and merge
     it!

Adding atoms
------------

Here are the steps to add a new function or operation (atom) to
Convex.jl. Let's say you're adding the new function $f$.

 -   Take a look at the [nuclear norm
     atom](https://github.com/JuliaOpt/Convex.jl/blob/master/src/atoms/sdp_cone/nuclearnorm.jl)
     for an example of how to construct atoms, and see the [norm
     atom](https://github.com/JuliaOpt/Convex.jl/blob/master/src/atoms/second_order_cone/norm.jl)
     for an example of an atom that depends on a parameter.
 -   Copy paste (eg) the nuclear norm file, replace anything saying
     nuclear norm with the name of the atom $f$, fill in monotonicity,
     curvature, etc. Save it in the appropriate subfolder of
     `src/atoms/`.
 -   Add as a comment a description of what the atom does and its
     parameters.
 -   The most mathematically interesting part is the `conic_form!`
     function. Following the example in the nuclear norm atom, you'll
     see that you can just construct the problem whose optimal value is
     $f(x)$, introducing any auxiliary variables you need, exactly as
     you would normally in Convex.jl, and then call `cache_conic_form!`
     on that problem.
 -   Add a test for the atom so we can verify it works in
     `src/problem_depot/problem/<cone>`, where `<cone>` matches the subfolder of
     `src/atoms`. See [How to write a ProblemDepot problem](@ref) for details
     on how to write the tests.
 -   Submit a PR, including a description of what the atom does and its
     parameters. (Let us know if it's a work in progress by putting
     \[WIP\] in the name of the PR.)
 -   We'll look it over, fix up anything that doesn't work, and merge
     it!

Fixing the guts
---------------

If you want to do a more major bug fix, you may need to understand how
Convex.jl thinks about conic form. To do this, start by reading [the
Convex.jl paper](http://arxiv.org/pdf/1410.4821.pdf). You may find our
[JuliaCon 2014 talk](https://www.youtube.com/watch?v=SoI0lEaUvTs&t=128s)
helpful as well; you can find the ipython notebook presented in the talk
[here](https://github.com/JuliaCon/presentations/tree/master/CVX).

Then read the conic form code:

 -   We define data structures for conic objectives and conic
     constraints, and simple ways of combining them, in
     [conic\_form.jl](https://github.com/JuliaOpt/Convex.jl/blob/master/src/conic_form.jl)
 -   We load the internal conic form representation into the
     [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl)
     model in the function
     [load\_MOI\_model!](https://github.com/JuliaOpt/Convex.jl/blob/master/src/solution.jl#L151).
 -   We solve problems (that is, pass the standard form of the problem
     to a solver, and put the solution back into the values of the
     appropriate variables) in
     [solve!](https://github.com/JuliaOpt/Convex.jl/blob/master/src/solution.jl#L205).

You're now armed and dangerous. Go ahead and open an issue (or comment
on a previous one) if you can't figure something out, or submit a PR if
you can figure it out. (Let us know if it's a work in progress by
putting \[WIP\] in the name of the PR.)

PRs that comment the code more thoroughly will also be welcomed.
