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

    # Complex Positive Semidefinite matrix
    z = HermitianSemidefinite(4)


Linear Program Representable Functions
**************************************

Apart from all the linear functions that are listed `here <http://convexjl.readthedocs.io/en/latest/operations.html#linear-program-representable-functionsl>`_, we have added new functions:

+--------------------------+-------------------------+------------+---------------+---------------------------------+
|operation                 | description             | vexity     | slope         | notes                           |
+==========================+=========================+============+===============+=================================+
|:code:`real(z)`           | real part of complex    | affine     |increasing     | none                            |
|                          | variable                |            |               |                                 |
+--------------------------+-------------------------+------------+---------------+---------------------------------+
|:code:`imag(z)`           | imaginary part of       | affine     |increasing in  | none                            |
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


Exponential Cone  Representable Functions
******************************************

An optimization problem using these functions can be solved by any exponential cone solver (SCS).

+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|operation                   | description                         | vexity     | slope         | notes                    |
+============================+=====================================+============+===============+==========================+
|:code:`logsumexp(x)`        | :math:`\log(\sum_i \exp(x_i))`      | convex     |increasing     |none                      |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`exp(x)`              | :math:`\exp(x)`                     | convex     |increasing     | none                     |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`log(x)`              | :math:`\log(x)`                     | concave    |increasing     | IC: :math:`x>0`          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`entropy(x)`          | :math:`\sum_{ij}                    | concave    |not monotonic  | IC: :math:`x>0`          |
|                            | -x_{ij} \log (x_{ij})`              |            |               |                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`logisticloss(x)`     | :math:`\log(1 + \exp(x_i))`         | convex     |increasing     | none                     |
|                            |                                     |            |               |                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+


Semidefinite Program Representable Functions
********************************************

An optimization problem using these functions can be solved by any SDP solver (including SCS and Mosek).

+---------------------------+-------------------------------------+------------+---------------+------------------------------+
|operation                  | description                         | vexity     | slope         | notes                        |
+===========================+=====================================+============+===============+==============================+
|:code:`nuclearnorm(x)`     | sum of singular values of :math:`x` | convex     |not monotonic  | none                         |
+---------------------------+-------------------------------------+------------+---------------+------------------------------+
|:code:`operatornorm(x)`    | max of singular values of :math:`x` | convex     |not monotonic  | none                         |
+---------------------------+-------------------------------------+------------+---------------+------------------------------+
|:code:`lambdamax(x)`       | max eigenvalue of :math:`x`         | convex     |not monotonic  |IC: x is positive semidefinite|
+---------------------------+-------------------------------------+------------+---------------+------------------------------+
|:code:`lambdamin(x)`       | min eigenvalue of :math:`x`         | concave    |not monotonic  |IC: x is positive semidefinite|
+---------------------------+-------------------------------------+------------+---------------+------------------------------+
|:code:`matrixfrac(x, P)`   | :math:`x^TP^{-1}x`                  | convex     |not monotonic  |IC: P is positive semidefinite|
+---------------------------+-------------------------------------+------------+---------------+------------------------------+

Exponential + SDP representable Functions
********************************************

An optimization problem using these functions can be solved by any solver that supports exponential constraints *and* semidefinite constraints simultaneously (SCS).

+----------------------------+-------------------------------------+------------+---------------+------------------------------+
|operation                   | description                         | vexity     | slope         | notes                        |
+============================+=====================================+============+===============+==============================+
|:code:`logdet(x)`           | log of determinant of :math:`x`     | concave    |increasing     |IC: x is positive semidefinite|
+----------------------------+-------------------------------------+------------+---------------+------------------------------+

Promotions
***********

When an atom or constraint is applied to a scalar and a higher dimensional variable, the scalars are promoted. For example, we can do :code:`max(x, 0)` gives an expression with the shape of :code:`x` whose elements are the maximum of the corresponding element of :code:`x` and :code:`0`.
