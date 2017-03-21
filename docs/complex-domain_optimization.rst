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

All of the linear functions that are listed `here <operations.html#linear-program-representable-functions>`_ operate on
complex variables as well. In addition, several specialized functions for complex variables are available:

+--------------------------+-------------------------+------------+---------------+---------------------------------+
|operation                 | description             | vexity     | slope         | notes                           |
+==========================+=========================+============+===============+=================================+
|:code:`real(z)`           | real part of complex    | affine     |increasing     | none                            |
|                          | variable                |            |               |                                 |
+--------------------------+-------------------------+------------+---------------+---------------------------------+
|:code:`imag(z)`           | imaginary part of       | affine     |increasing     | none                            |
|                          | complex variable        |            |               |                                 |
+--------------------------+-------------------------+------------+---------------+---------------------------------+
|:code:`conj(x)'`          | complex conjugate       | affine     |increasing     | none                            |
+--------------------------+-------------------------+------------+---------------+---------------------------------+
|:code:`innerproduct(x,y)` | real(trace(x'*y))       | affine     |increasing     | PR: one argument is constant    |
+--------------------------+-------------------------+------------+---------------+---------------------------------+


Second-Order Cone Representable Functions
*****************************************

Most of the second order cone function listed `here <operations.html#second-order-cone-representable-functions>`_ operate on
complex variables as well. Notable exceptions include:

  * inverse 
  * square 
  * quadoverlin
  * sqrt
  * geomean
  * huber

One new function is available:  

+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|operation                   | description                         | vexity     | slope         | notes                    |
+============================+=====================================+============+===============+==========================+
|:code:`abs2(z)`             | `square(abs(z))`                    | convex     |increasing     | none                     |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+


Semidefinite Program Representable Functions
********************************************

All SDP-representable functions listed `here <operations.html#semidefinite-program-representable-functions>`_ work for complex variables.


Exponential + SDP representable Functions
********************************************

Complex variables also support logdet function.

