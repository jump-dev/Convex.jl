Advanced Features
=================

Dual Variables
--------------

Convex.jl also returns the optimal dual variables for a problem. These
are stored in the `dual` field associated with each constraint. :

    using Convex, SCS

    x = Variable()
    constraint = x >= 0
    p = minimize(x, constraint)
    solve!(p, SCSSolver())

    # Get the dual value for the constraint
    p.constraints[1].dual
    # or
    constraint.dual

Warmstarting
------------

If you're solving the same problem many times with different values of
a parameter, Convex.jl can initialize many solvers with the solution to
the previous problem, which sometimes speeds up the solution time. This
is called a **warm start**.

To use this feature, pass the optional argument
[warmstart=true]{.title-ref} to the [solve!]{.title-ref} method. :

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

Fixing and freeing variables
----------------------------

Convex.jl allows you to fix a variable `x` to a value by
calling the [fix!]{.title-ref} method. Fixing the variable essentially
turns it into a constant. Fixed variables are sometimes also called
parameters.

[fix(x, v)]{.title-ref} fixes the variable `x` to the value
`v`.

[fix(x)]{.title-ref} fixes `x` to the value
[x.value]{.title-ref}, which might be the value obtained by solving
another problem involving the variable `x`.

To allow the variable `x` to vary again, call
[free!(x)]{.title-ref}.

Fixing and freeing variables can be particularly useful as a tool for
performing alternating minimization on nonconvex problems. For example,
we can find an approximate solution to a nonnegative matrix
factorization problem with alternating minimization as follows. We use
warmstarts to speed up the solution. :

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
