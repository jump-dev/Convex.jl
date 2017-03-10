=====================================
Contributing
=====================================

We'd welcome contributions to the Convex.jl package. Here are some short instructions on how to get started. If you don't know what you'd like to contribute, you could 

	* take a look at the current `issues <https://github.com/JuliaOpt/Convex.jl/issues>`_ and pick one. (Feature requests are probably the easiest to tackle.)
	* add a `usage example <https://github.com/JuliaOpt/Convex.jl/tree/master/examples>`_.

Then submit a pull request (PR). (Let us know if it's a work in progress by putting [WIP] in the name of the PR.)

Adding examples
***************

	* Take a look at our exising `usage examples <https://github.com/JuliaOpt/Convex.jl/tree/master/examples>`_ and add another in similar style. 
	* Submit a PR. (Let us know if it's a work in progress by putting [WIP] in the name of the PR.)
	* We'll look it over, fix up anything that doesn't work, and merge it!

Adding  atoms
*************************************

Here are the steps to add a new function or operation (atom) to Convex.jl. Let's say you're
adding the new function :math:`f`.

	* Take a look at the `nuclear norm atom <https://github.com/JuliaOpt/Convex.jl/blob/master/src/atoms/sdp_cone/nuclearnorm.jl>`_ for an example of how to construct atoms, and see the `norm atom <https://github.com/JuliaOpt/Convex.jl/blob/master/src/atoms/norm.jl>`_ for an example of an atom that depends on a parameter.
	* Copy paste (eg) the nuclear norm file, replace anything saying nuclear norm with the name of the atom :math:`f`, fill in monotonicity, curvature, etc. Save it in the appropriate subfolder of :code:`src/atoms/`. 
	* Add as a comment a description of what the atom does and its parameters.
	* The most mathematically interesting part is the :code:`conic_form!` function. Following the example in the nuclear norm atom, you'll see that you can just construct the problem whose optimal value is :math:`f(x)`, introducing any auxiliary variables you need, exactly as you would normally in Convex.jl, and then call :code:`cache_conic_form!` on that problem.
	* Add a test for the atom so we can verify it works in :code:`test/test_<cone>`, where :code:`<cone>` matches the subfolder of :code:`src/atoms`.
	* Submit a PR, including a description of what the atom does and its parameters. (Let us know if it's a work in progress by putting [WIP] in the name of the PR.)
	* We'll look it over, fix up anything that doesn't work, and merge it!

Fixing the guts
***************

If you want to do a more major bug fix, you may need to understand how Convex.jl 
thinks about conic form. To do this, start by reading 
`the Convex.jl paper <http://arxiv.org/pdf/1410.4821.pdf>`_.
You may find our `JuliaCon 2014 talk <https://www.youtube.com/watch?v=SoI0lEaUvTs&t=128s>`_ helpful as well; you can find the ipython notebook presented in the talk `here <https://github.com/JuliaCon/presentations/tree/master/CVX>`_.

Then read the conic form code: 

	* We define data structures for conic objectives and conic constraints, and simple ways of combining them, in `conic_form.jl <https://github.com/JuliaOpt/Convex.jl/blob/master/src/conic_form.jl>`_
	* We convert the internal conic form representation into the `standard form for conic solvers <http://mathprogbasejl.readthedocs.io/en/latest/conic.html>`_ in the function `conic_problem <https://github.com/JuliaOpt/Convex.jl/blob/master/src/problems.jl#L97>`_.
	* We solve problems (that is, pass the standard form of the problem to a solver, and put the solution back into the values of the appropriate variables) in `solution.jl <https://github.com/JuliaOpt/Convex.jl/blob/master/src/solution.jl#L8>`_.

You're now armed and dangerous. Go ahead and open an issue (or comment on a previous one) if you can't figure something out, or submit a PR if you can figure it out. (Let us know if it's a work in progress by putting [WIP] in the name of the PR.) 

PRs that comment the code more thoroughly will also be welcomed.
