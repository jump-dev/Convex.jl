Operations
==========

Convex.jl currently supports the following functions. These functions
may be composed according to the [DCP](http://dcp.stanford.edu)
composition rules to form new convex, concave, or affine expressions.
Convex.jl transforms each problem into an equivalent [cone
program](http://mathprogbasejl.readthedocs.org/en/latest/conic.html) in
order to pass the problem to a specialized solver. Depending on the
types of functions used in the problem, the conic constraints may
include linear, second-order, exponential, or semidefinite constraints,
as well as any binary or integer constraints placed on the variables.
Below, we list each function available in Convex.jl organized by the
(most complex) type of cone used to represent that function, and
indicate which solvers may be used to solve problems with those cones.
Problems mixing many different conic constraints can be solved by any
solver that supports every kind of cone present in the problem.

In the notes column in the tables below, we denote implicit constraints
imposed on the arguments to the function by IC, and parameter
restrictions that the arguments must obey by PR. (Convex.jl will
automatically impose ICs; the user must make sure to satisfy PRs.)
Elementwise means that the function operates elementwise on vector
arguments, returning a vector of the same size.

Linear Program Representable Functions
--------------------------------------

An optimization problem using only these functions can be solved by any
LP solver.

+---------------+----------------+------+--------+-------------------+
| operation     | description    | vexi | slope  | notes             |
|               |                | ty   |        |                   |
+===============+================+======+========+===================+
| `x+y` or      | addition       | affi | increa | none              |
| `x.+y`        |                | ne   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `x-y` or      | subtraction    | affi | increa | none              |
| `x.-y`        |                | ne   | sing   |                   |
|               |                |      | in $x$ | none              |
|               |                |      |        |                   |
|               |                |      | decrea |                   |
|               |                |      | sing   |                   |
|               |                |      | in $y$ |                   |
+---------------+----------------+------+--------+-------------------+
| `x*y`         | multiplication | affi | increa | PR: one argument  |
|               |                | ne   | sing   | is constant       |
|               |                |      | if     |                   |
|               |                |      |        |                   |
|               |                |      | consta |                   |
|               |                |      | nt     |                   |
|               |                |      | term   |                   |
|               |                |      | $\ge 0 |                   |
|               |                |      | $      |                   |
|               |                |      |        |                   |
|               |                |      | decrea |                   |
|               |                |      | sing   |                   |
|               |                |      | if     |                   |
|               |                |      |        |                   |
|               |                |      | consta |                   |
|               |                |      | nt     |                   |
|               |                |      | term   |                   |
|               |                |      | $\le 0 |                   |
|               |                |      | $      |                   |
|               |                |      |        |                   |
|               |                |      | not    |                   |
|               |                |      | monoto |                   |
|               |                |      | nic    |                   |
|               |                |      |        |                   |
|               |                |      | otherw |                   |
|               |                |      | ise    |                   |
+---------------+----------------+------+--------+-------------------+
| `x/y`         | division       | affi | increa | PR: $y$ is scalar |
|               |                | ne   | sing   | constant          |
+---------------+----------------+------+--------+-------------------+
| `dot(*)(x, y) | elementwise    | affi | increa | PR: one argument  |
| `             | multiplication | ne   | sing   | is constant       |
+---------------+----------------+------+--------+-------------------+
| `dot(/)(x, y) | elementwise    | affi | increa | PR: one argument  |
| `             | division       | ne   | sing   | is constant       |
+---------------+----------------+------+--------+-------------------+
| `x[1:4, 2:3]` | indexing and   | affi | increa | none              |
|               | slicing        | ne   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `diag(x, k)`  | $k$-th         | affi | increa | none              |
|               | diagonal of a  | ne   | sing   |                   |
|               | matrix         |      |        |                   |
+---------------+----------------+------+--------+-------------------+
| `diagm(x)`    | construct      | affi | increa | PR: $x$ is a      |
|               | diagonal       | ne   | sing   | vector            |
|               | matrix         |      |        |                   |
+---------------+----------------+------+--------+-------------------+
| `x'`          | transpose      | affi | increa | none              |
|               |                | ne   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `vec(x)`      | vector         | affi | increa | none              |
|               | representation | ne   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `dot(x,y)`    | $\sum_i x_i y_ | affi | increa | PR: one argument  |
|               | i$             | ne   | sing   | is constant       |
+---------------+----------------+------+--------+-------------------+
| `kron(x,y)`   | Kronecker      | affi | increa | PR: one argument  |
|               | product        | ne   | sing   | is constant       |
+---------------+----------------+------+--------+-------------------+
| `vecdot(x,y)` | `dot(vec(x),ve | affi | increa | PR: one argument  |
|               | c(y))`         | ne   | sing   | is constant       |
+---------------+----------------+------+--------+-------------------+
| `sum(x)`      | $\sum_{ij} x_{ | affi | increa | none              |
|               | ij}$           | ne   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `sum(x, k)`   | sum elements   | affi | increa | none              |
|               | across         | ne   | sing   |                   |
|               |                |      |        |                   |
|               | dimension $k$  |      |        |                   |
+---------------+----------------+------+--------+-------------------+
| `sumlargest(x | sum of $k$     | conv | increa | none              |
| , k)`         | largest        | ex   | sing   |                   |
|               |                |      |        |                   |
|               | elements of    |      |        |                   |
|               | $x$            |      |        |                   |
+---------------+----------------+------+--------+-------------------+
| `sumsmallest( | sum of $k$     | conc | increa | none              |
| x, k)`        | smallest       | ave  | sing   |                   |
|               |                |      |        |                   |
|               | elements of    |      |        |                   |
|               | $x$            |      |        |                   |
+---------------+----------------+------+--------+-------------------+
| `dotsort(a, b | `dot(sort(a),s | conv | increa | PR: one argument  |
| )`            | ort(b))`       | ex   | sing   | is constant       |
+---------------+----------------+------+--------+-------------------+
| `reshape(x, m | reshape into   | affi | increa | none              |
| , n)`         | $m \times n$   | ne   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `minimum(x)`  | $\min(x)$      | conc | increa | none              |
|               |                | ave  | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `maximum(x)`  | $\max(x)$      | conv | increa | none              |
|               |                | ex   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `[x y]` or    | stacking       | affi | increa | none              |
| [\[x;         |                | ne   | sing   |                   |
| y\]]{.title-r |                |      |        |                   |
| ef}           |                |      |        |                   |
|               |                |      |        |                   |
| `hcat(x, y)`  |                |      |        |                   |
| or            |                |      |        |                   |
|               |                |      |        |                   |
| `vcat(x, y)`  |                |      |        |                   |
+---------------+----------------+------+--------+-------------------+
| `trace(x)`    | $\mathrm{tr}   | affi | increa | none              |
|               | \left(X \right | ne   | sing   |                   |
|               | )$             |      |        |                   |
+---------------+----------------+------+--------+-------------------+
| partialtrace( | Partial trace  | affi | increa | none              |
| x,sys,dims)   |                | ne   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `conv(h,x)`   | $h \in         | affi | increa | PR: $h$ is        |
|               | \mathbb{R}^m$  | ne   | sing   | constant          |
|               |                |      | if     |                   |
|               | $x \in         |      | $h\ge  |                   |
|               | \mathbb{R}^m$  |      | 0$     |                   |
|               |                |      |        |                   |
|               | $h*x           |      | decrea |                   |
|               | \in \mathbb{R} |      | sing   |                   |
|               | ^{m+n-1}$      |      | if     |                   |
|               |                |      | $h\le  |                   |
|               | entry $i$ is   |      | 0$     |                   |
|               | given by       |      |        |                   |
|               |                |      | not    |                   |
|               | $\sum_{j=1}^m  |      | monoto |                   |
|               | h_jx_{i-j}$    |      | nic    |                   |
|               |                |      |        |                   |
|               |                |      | otherw |                   |
|               |                |      | ise    |                   |
+---------------+----------------+------+--------+-------------------+
| `min(x,y)`    | $\min(x,y)$    | conc | increa | none              |
|               |                | ave  | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `max(x,y)`    | $\max(x,y)$    | conv | increa | none              |
|               |                | ex   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `pos(x)`      | $\max(x,0)$    | conv | increa | none              |
|               |                | ex   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `neg(x)`      | $\max(-x,0)$   | conv | decrea | none              |
|               |                | ex   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `invpos(x)`   | $1/x$          | conv | decrea | IC: $x>0$         |
|               |                | ex   | sing   |                   |
+---------------+----------------+------+--------+-------------------+
| `abs(x)`      | $\left|x\right | conv | increa | none              |
|               | |$             | ex   | sing   |                   |
|               |                |      | on     |                   |
|               |                |      | $x \ge |                   |
|               |                |      |  0$    |                   |
|               |                |      |        |                   |
|               |                |      | decrea |                   |
|               |                |      | sing   |                   |
|               |                |      | on     |                   |
|               |                |      | $x \le |                   |
|               |                |      |  0$    |                   |
+---------------+----------------+------+--------+-------------------+

Second-Order Cone Representable Functions
-----------------------------------------

An optimization problem using these functions can be solved by any SOCP
solver (including ECOS, SCS, Mosek, Gurobi, and CPLEX). Of course, if an
optimization problem has both LP and SOCP representable functions, then
any solver that can solve both LPs and SOCPs can solve the problem.

+---------------+---------------------+------+--------+--------------+
| operation     | description         | vexi | slope  | notes        |
|               |                     | ty   |        |              |
+===============+=====================+======+========+==============+
| `norm(x, p)`  | $(\sum x_i^p)^{1/p} | conv | increa | PR: `p >= 1` |
|               | $                   | ex   | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \ge |              |
|               |                     |      |  0$    |              |
|               |                     |      |        |              |
|               |                     |      | decrea |              |
|               |                     |      | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \le |              |
|               |                     |      |  0$    |              |
+---------------+---------------------+------+--------+--------------+
| `vecnorm(x, p | $(\sum x_{ij}^p)^{1 | conv | increa | PR: `p >= 1` |
| )`            | /p}$                | ex   | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \ge |              |
|               |                     |      |  0$    |              |
|               |                     |      |        |              |
|               |                     |      | decrea |              |
|               |                     |      | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \le |              |
|               |                     |      |  0$    |              |
+---------------+---------------------+------+--------+--------------+
| `quadform(x,  | $x^T P x$           | conv | increa | PR: either   |
| P)`           |                     | ex   | sing   | $x$ or $P$   |
|               |                     | in   | on     |              |
|               |                     | $x$  | $x \ge | must be      |
|               |                     |      |  0$    | constant; if |
|               |                     | affi |        | $x$ is not   |
|               |                     | ne   | decrea | constant,    |
|               |                     | in   | sing   | then $P$     |
|               |                     | $P$  | on     | must be      |
|               |                     |      | $x \le | symmetric    |
|               |                     |      |  0$    | and positive |
|               |                     |      |        | semidefinite |
|               |                     |      | increa |              |
|               |                     |      | sing   |              |
|               |                     |      | in $P$ |              |
+---------------+---------------------+------+--------+--------------+
| `quadoverlin( | $x^T x/y$           | conv | increa | IC: $y > 0$  |
| x, y)`        |                     | ex   | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \ge |              |
|               |                     |      |  0$    |              |
|               |                     |      |        |              |
|               |                     |      | decrea |              |
|               |                     |      | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \le |              |
|               |                     |      |  0$    |              |
|               |                     |      |        |              |
|               |                     |      | decrea |              |
|               |                     |      | sing   |              |
|               |                     |      | in $y$ |              |
+---------------+---------------------+------+--------+--------------+
| `sumsquares(x | $\sum x_i^2$        | conv | increa | none         |
| )`            |                     | ex   | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \ge |              |
|               |                     |      |  0$    |              |
|               |                     |      |        |              |
|               |                     |      | decrea |              |
|               |                     |      | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \le |              |
|               |                     |      |  0$    |              |
+---------------+---------------------+------+--------+--------------+
| `sqrt(x)`     | $\sqrt{x}$          | conc | decrea | IC: $x>0$    |
|               |                     | ave  | sing   |              |
+---------------+---------------------+------+--------+--------------+
| `square(x), x | $x^2$               | conv | increa | PR : $x$ is  |
| ^2`           |                     | ex   | sing   | scalar       |
|               |                     |      | on     |              |
|               |                     |      | $x \ge |              |
|               |                     |      |  0$    |              |
|               |                     |      |        |              |
|               |                     |      | decrea |              |
|               |                     |      | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \le |              |
|               |                     |      |  0$    |              |
+---------------+---------------------+------+--------+--------------+
| `dot(^)(x,2)` | $x.^2$              | conv | increa | elementwise  |
|               |                     | ex   | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \ge |              |
|               |                     |      |  0$    |              |
|               |                     |      | decrea |              |
|               |                     |      | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \le |              |
|               |                     |      |  0$    |              |
+---------------+---------------------+------+--------+--------------+
| `geomean(x, y | $\sqrt{xy}$         | conc | increa | IC: $x\ge0$, |
| )`            |                     | ave  | sing   | $y\ge0$      |
+---------------+---------------------+------+--------+--------------+
| `huber(x, M=1 | $\begin{cases}      | conv | increa | PR: $M>=1$   |
| )`            | x^2 &|x| \leq       | ex   | sing   |              |
|               | M  \\               |      | on     |              |
|               | 2M|x| - M^2         |      | $x \ge |              |
|               | &|x| >  M           |      |  0$    |              |
|               | \end{cases}$        |      |        |              |
|               |                     |      | decrea |              |
|               |                     |      | sing   |              |
|               |                     |      | on     |              |
|               |                     |      | $x \le |              |
|               |                     |      |  0$    |              |
+---------------+---------------------+------+--------+--------------+

Exponential Cone Representable Functions
----------------------------------------

An optimization problem using these functions can be solved by any
exponential cone solver (SCS).

  -------------------------------------------------------------------------------------
  operation           description                vexity    slope        notes
  ------------------- -------------------------- --------- ------------ ---------------
  `logsumexp(x)`      $\log(\sum_i \exp(x_i))$   convex    increasing   none

  `exp(x)`            $\exp(x)$                  convex    increasing   none

  `log(x)`            $\log(x)$                  concave   increasing   IC: $x>0$

  `entropy(x)`        $\sum_{ij}                 concave   not          IC: $x>0$
                      -x_{ij} \log (x_{ij})$               monotonic    

  `logisticloss(x)`   $\log(1 + \exp(x_i))$      convex    increasing   none
  -------------------------------------------------------------------------------------

Semidefinite Program Representable Functions
--------------------------------------------

An optimization problem using these functions can be solved by any SDP
solver (including SCS and Mosek).

  ----------------------------------------------------------------------------------
  operation            description           vexity    slope       notes
  -------------------- --------------------- --------- ----------- -----------------
  `nuclearnorm(x)`     sum of singular       convex    not         none
                       values of $x$                   monotonic   

  `operatornorm(x)`    max of singular       convex    not         none
                       values of $x$                   monotonic   

  `lambdamax(x)`       max eigenvalue of $x$ convex    not         none
                                                       monotonic   

  `lambdamin(x)`       min eigenvalue of $x$ concave   not         none
                                                       monotonic   

  `matrixfrac(x, P)`   $x^TP^{-1}x$          convex    not         IC: P is positive
                                                       monotonic   semidefinite
  ----------------------------------------------------------------------------------

Exponential + SDP representable Functions
-----------------------------------------

An optimization problem using these functions can be solved by any
solver that supports exponential constraints *and* semidefinite
constraints simultaneously (SCS).

  -------------------------------------------------------------------------------
  operation        description           vexity    slope        notes
  ---------------- --------------------- --------- ------------ -----------------
  `logdet(x)`      log of determinant of concave   increasing   IC: x is positive
                   $x$                                          semidefinite

  -------------------------------------------------------------------------------

Promotions
----------

When an atom or constraint is applied to a scalar and a higher
dimensional variable, the scalars are promoted. For example, we can do
`max(x, 0)` gives an expression with the shape of `x` whose elements are
the maximum of the corresponding element of `x` and `0`.
