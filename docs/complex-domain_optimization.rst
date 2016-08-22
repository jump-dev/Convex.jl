=====================================
Complex-domain Optimization
=====================================

Convex.jl also supports optimization over the complex domain.
The idea is to transform the complex-domain problem to the corresponding real-domain problem using a bijective mapping and then use the existing machinery to solve the reduced real-domain problem and then combine the result to return the complex solution.
We support the all the conic constraints such as linear, second-order, or semidefinite constraints for complex variables as well.

Below, we present a quick start guide on how to use Convex.jl for complex-domain optimization and then list the operation supported on complex variables in Convex.jl organized by the (most complex) type of cone used to represent that operation.

Complex Variables
**************************************
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

Apart from all the linear functions that are listed `here <operations.html#linear-program-representable-functionsl>`_, we have added new functions:

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



+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|operation                   | description                         | vexity     | slope         | notes                    |
+============================+=====================================+============+===============+==========================+
|:code:`norm(x, p)`          | :math:`(\sum x_i^p)^{1/p}`          | convex     |increasing on  | PR: :code:`p >= 1`       |
|                            |                                     |            |:math:`x \ge 0`|                          |
|                            |                                     |            |               |                          |
|                            |                                     |            |decreasing on  |                          |
|                            |                                     |            |:math:`x \le 0`|                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`vecnorm(x, p)`       | :math:`(\sum x_{ij}^p)^{1/p}`       | convex     |increasing on  | PR: :code:`p >= 1`       |
|                            |                                     |            |:math:`x \ge 0`|                          |
|                            |                                     |            |               |                          |
|                            |                                     |            |decreasing on  |                          |
|                            |                                     |            |:math:`x \le 0`|                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`quadform(x, P)`      | :math:`x^T P x`                     | convex in  |increasing on  | PR: either :math:`x` or  |
|                            |                                     | :math:`x`  |:math:`x \ge 0`| :math:`P`                |
|                            |                                     |            |               |                          |
|                            |                                     | affine in  |decreasing on  | must be constant;        |
|                            |                                     | :math:`P`  |:math:`x \le 0`| if :math:`x` is not      |
|                            |                                     |            |               | constant, then :math:`P` |
|                            |                                     |            |increasing in  | must be symmetric and    |
|                            |                                     |            |:math:`P`      | positive semidefinite    |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`quadoverlin(x, y)`   | :math:`x^T x/y`                     | convex     |increasing on  |                          |
|                            |                                     |            |:math:`x \ge 0`| IC: :math:`y > 0`        |
|                            |                                     |            |               |                          |
|                            |                                     |            |decreasing on  |                          |
|                            |                                     |            |:math:`x \le 0`|                          |
|                            |                                     |            |               |                          |
|                            |                                     |            |decreasing in  |                          |
|                            |                                     |            |:math:`y`      |                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`sumsquares(x)`       | :math:`\sum x_i^2`                  | convex     |increasing on  | none                     |
|                            |                                     |            |:math:`x \ge 0`|                          |
|                            |                                     |            |               |                          |
|                            |                                     |            |decreasing on  |                          |
|                            |                                     |            |:math:`x \le 0`|                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`sqrt(x)`             | :math:`\sqrt{x}`                    | concave    |decreasing     | IC: :math:`x>0`          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`square(x), x^2`      | :math:`x^2`                         | convex     |increasing on  | none                     |
|                            |                                     |            |:math:`x \ge 0`|                          |
|                            |                                     |            |               |                          |
|                            |                                     |            |decreasing on  |                          |
|                            |                                     |            |:math:`x \le 0`|                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`geomean(x, y)`       | :math:`\sqrt{xy}`                   | concave    |increasing     | IC: :math:`x\ge0`,       |
|                            |                                     |            |               | :math:`y\ge0`            |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`huber(x)`            | :math:`\begin{cases}                | convex     |increasing on  | PR: :math:`M>=1`         |
|                            | x^2 &|x| \leq                       |            |:math:`x \ge 0`|                          |
|:code:`huber(x, M)`         | M  \\                               |            |               |                          |
|                            | 2M|x| - M^2                         |            |               |                          |
|                            | &|x| >  M                           |            |decreasing on  |                          |
|                            | \end{cases}`                        |            |:math:`x \le 0`|                          |
|                            |                                     |            |               |                          |
|                            |                                     |            |               |                          |
|                            |                                     |            |               |                          |
|                            |                                     |            |               |                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+




Semidefinite Program Representable Functions
********************************************

Complex variables support all the semidefinite functions same as real variables listed `here <operations.html#semidefinite-program-representable-functions>`_.


Exponential + SDP representable Functions
********************************************

Complex variables also support logdet function.

