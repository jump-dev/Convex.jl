Optimizing in a Loop
====================

Memory Management
-----------------

Convex uses a module-level dictionary to store the conic forms of every
variable and expression created in the same Julia session. These
variables and expressions persist even after they are out of scope. If
you create large numbers of variables inside a loop, this dictionary can
eat a considerable amount of memory.

To flush the memory, you can call :: Convex.clearmemory()

This will remove every variable and expression you've formed before
from the memory cache, so that you're starting as fresh as if you'd
just reimported Convex.

Caching Expressions
-------------------

Better yet, take advantage of this cache of variables and expressions!
Create variables and expressions outside the loop, and reuse them inside
the loop as you tweak parameters. Doing this will allow Convex to reuse
the conic forms it has already calculated for previously used
expressions.

For example, the following **bad** code will create a new instance of a
variable and of the expression `square(x)` for each value of `i`. Don't
do this: :: for i=1:10 x = Variable() p = minimize(square(x), x \>= i)
solve!(p, SCSSolver()) end

Contrast this with the following **good** code, which will reuse the
cached conic form for `square(x)` for each `i`, reducing the memory
footprint and speeding up the computation. Do this instead: :: x =
Variable() obj = square(x) for i=1:10 p = minimize(obj, x \>= i)
solve!(p, SCSSolver()) end

Warmstarts, Parameters, Fixing and Freeing Variables
----------------------------------------------------

If you're solving many problems of the same form, or many similar
problems, you may also want to use warmstarts, or to dynamically fix and
free variables. The former is particularly good for a family of problems
related by a parameter; the latter allows for easy implementation of
alternating minimization for nonconvex problems. See the [Advanced
Features](@ref) section of the documentation for more
information.
