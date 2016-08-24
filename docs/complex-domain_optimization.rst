=====================================
Optimization with Complex Variables
=====================================

Convex.jl also supports optimization with complex variables. Below, we present a quick start guide on how to use Convex.jl for optimization with complex variables, and then list the operations supported on complex variables in Convex.jl. In general, any operation available in Convex.jl that is well defined and DCP compliant on complex variables should be available. We list these functions below. organized by the type of cone (linear, second-order, or semidefinite) used to represent that operation.

Internally, Convex.jl transforms the complex-domain problem to a larger real-domain problem using a bijective mapping. It then solves the real-domain problem and transforms the solution back to the complex domain.

Complex Variables
*****************
Complex Variables in Convex.jl are declared in the same way as the variables are declared but using the different keyword `ComplexVariable`.
::

    # Scalar complex variable
    z = ComplexVariable()

    # Column vector variable
    z = ComplexVariable(5)

    # Matrix variable
    z = ComplexVariable(4, 6)

    # Complex Positive Semidefinite variable
    z = HermitianSemidefinite(4)


Linear Program Representable Functions
**************************************

Apart from all the linear functions that are listed `here <operations.html#linear-program-representable-functions>`_, we have added new functions:

+--------------------------+-------------------------+------------+---------------+---------------------------------+
|operation                 | description             | vexity     | slope         | notes                           |
+==========================+=========================+============+===============+=================================+
|:code:`real(z)`           | real part of complex    | affine     |increasing     | none                            |
|                          | variable                |            |               |                                 |
+--------------------------+-------------------------+------------+---------------+---------------------------------+
|:code:`imag(z)`           | imaginary part of       | affine     |increasing     | none                            |
|                          | complex variable        |            |               |                                 |
+--------------------------+-------------------------+------------+---------------+---------------------------------+
|:code:`x'`                | ctranspose              | affine     |increasing     | none                            |
+--------------------------+-------------------------+------------+---------------+---------------------------------+
|:code:`innerproduct(x,y)` | real(trace(x'*y))       | affine     |increasing     | PR: one argument is constant    |
+--------------------------+-------------------------+------------+---------------+---------------------------------+


Second-Order Cone Representable Functions
*****************************************

Most of the second order cone function listed here operate on complex variables as well except few listed below:

  * inverse 
  * square 
  * quadoverlin
  * sqrt
  * geomean
  * huber

We added a new function  

+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|operation                   | description                         | vexity     | slope         | notes                    |
+============================+=====================================+============+===============+==========================+
|:code:`squaremodulus(z)`    | `square(abs(z))`                    | convex     |increasing     | none                     |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+


Semidefinite Program Representable Functions
********************************************

Complex variables support all the semidefinite functions same as real variables listed `here <operations.html#semidefinite-program-representable-functions>`_.


Exponential + SDP representable Functions
********************************************

Complex variables also support logdet function.

