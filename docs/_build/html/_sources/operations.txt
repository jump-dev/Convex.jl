=====================================
Operations
=====================================

Convex.jl currently supports the following operations. These functions ("atoms") may be composed according to the `DCP <http://dcp.stanford.edu>`_ composition rules to form new convex, concave, or affine expressions.

Affine atoms
*************

These atoms are affine in their arguments.

  * addition, subtraction, multiplication, division: :code:`+, -, /, *`
  * indexing into vectors and matrices: :code:`x[1:4, 2:3]`
  * k-th diagonal of a matrix: :code:`diag(x, k)`
  * transpose: :code:`x'`
  * dot product: :code:`x' * y` or :code:`dot(x, y)`
  * reshape, vec: :code:`reshape(x, 2, 3)` or :code:`vec(x)`
  * minimum, maximum element of a vector or matrix: :code:`maximum(x)`
  * horizontal and vertical stacking: :code:`hcat(x, y); vcat(x, y)`
  * trace: :code:`trace(X)`

Elementwise atoms
*******************

These atoms are applied elementwise to their arguments, returning an expression of the same size.

+----------------------+-------------------------+------------+---------------+-------------------------+
|atom                  | description             | vexity     | slope         | implicit constraint     |
+======================+=========================+============+===============+=========================+
|:code:`min(x,y)`      | :math:`\min(x,y)`       | concave    |increasing     | none                    |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`max(x,y)`      | :math:`\max(x,y)`       | convex     |increasing     | none                    |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`pos(x,y)`      | :math:`\max(x,0)`       | convex     |increasing     | none                    |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`neg(x,y)`      | :math:`\max(-x,0)`      | convex     |decreasing     | none                    |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`inv_pos(x)`    | :math:`1/\max(x,0)`     | convex     |decreasing     | :math:`x>0`             |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`sqrt(x)`       | :math:`\sqrt{x}`        | convex     |decreasing     | :math:`x>0`             |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`square(x), x^2`| :math:`x^2`             | convex     |increasing on  | none                    |
|                      |                         |            |:math:`x \ge 0`|                         |
|                      |                         |            |decreasing on  |                         |
|                      |                         |            |:math:`x \le 0`|                         |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`abs(x)`        | :math:`\left|x\right|`  | convex     |increasing on  | none                    |
|                      |                         |            |:math:`x \ge 0`|                         |
|                      |                         |            |decreasing on  |                         |
|                      |                         |            |:math:`x \le 0`|                         |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`geo_mean(x, y)`| :math:`\sqrt{xy}`       | concave    |increasing     | :math:`x\ge0`,          |
|                      |                         |            |               | :math:`y\ge0`           |
|                      |                         |            |               |                         |
|                      |                         |            |               |                         |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`exp(x)`        | :math:`\exp(x)`         | convex     |increasing     | none                    |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`log(x)`        | :math:`\log(x)`         | concave    |increasing     | :math:`x>0`             |
+----------------------+-------------------------+------------+---------------+-------------------------+
|:code:`huber(x)`      | :math:`\begin{cases}    | convex     |increasing on  | :math:`M>=1`            |
|:code:`huber(x, M)`   | x^2 &|x| \leq           |            |:math:`x \ge 0`|                         |
|                      | M  \\                   |            |decreasing on  |                         |
|                      | 2M|x| - M^2             |            |:math:`x \le 0`|                         |
|                      | &|x| >                  |            |               |                         |
|                      | M                       |            |               |                         |
|                      | \end{cases}`            |            |               |                         |
+----------------------+-------------------------+------------+---------------+-------------------------+

Vector and Matrix atoms
**************************

These atoms take vector or matrix arguments and return scalar expressions.

+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|atom                        | description                         | vexity     | slope         | implicit constraint      |
+============================+=====================================+============+===============+==========================+
|:code:`norm(x, p)`          | :math:`(\sum x_i^p)^{1/p}`          | convex     |increasing on  | :code:`p = 1, 2, Inf`    |
|                            |                                     |            |:math:`x \ge 0`|                          |
|                            |                                     |            |decreasing on  |                          |
|                            |                                     |            |:math:`x \le 0`|                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`vecnorm(x, p)`       | :math:`(\sum x_i^p)^{1/p}`          | convex     |increasing on  | :code:`p = 1, 2, Inf`    |
|                            |                                     |            |:math:`x \ge 0`|                          |
|                            |                                     |            |decreasing on  |                          |
|                            |                                     |            |:math:`x \le 0`|                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`quad_form(x, P)`     | :math:`x^T P x`                     | convex in  |increasing on  | either :math:`x` or      |
|                            |                                     | :math:`x`  |:math:`x \ge 0`| :math:`P` must be        |
|                            |                                     | affine in  |decreasing on  | constant                 |
|                            |                                     | :math:`P`  |:math:`x \le 0`|                          |
|                            |                                     |            |increasing in  |                          |
|                            |                                     |            |:math:`P`      |                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`quad_over_lin(x, y)` | :math:`x^T x/y`                     | convex     |increasing on  |                          |
|                            |                                     |            |:math:`x \ge 0`| :math:`y > 0`            |
|                            |                                     |            |decreasing on  |                          |
|                            |                                     |            |:math:`x \le 0`|                          |
|                            |                                     |            |decreasing in  |                          |
|                            |                                     |            |:math:`y`      |                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`sum_squares(x)`      | :math:`\sum x_i^2`                  | convex     |increasing on  | none                     |
|                            |                                     |            |:math:`x \ge 0`|                          |
|                            |                                     |            |decreasing on  |                          |
|                            |                                     |            |:math:`x \le 0`|                          |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`nuclear_norm(x)`     | sum of singular values of :math:`x` | convex     |not monotonic  | none                     |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`operator_norm(x)`    | max of singular values of :math:`x` | convex     |not monotonic  | none                     |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`lambda_max(x)`       | max eigenvalue of :math:`x`         | convex     |increasing     |x is positive semidefinite|
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`lambda_min(x)`       | min of singular values of :math:`x` | concave    |increasing     |x is positive semidefinite|
+----------------------------+-------------------------------------+------------+---------------+--------------------------+
|:code:`logsumexp(x)`        | :math:`\log(\sum_i \exp(x_i))`      | convex     |increasing     |none                      |
+----------------------------+-------------------------------------+------------+---------------+--------------------------+

Promotions
***********

When an atom or constraint is applied to a scalar and a higher dimensional variable, the scalars are promoted. For example, we can do :code:`max(x, 0)` gives an expression with the shape of :code:`x` whose elements are the maximum of the corresponding element of :code:`x` and :code:`0`.
