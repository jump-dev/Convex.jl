Operations
==========

Convex.jl currently supports the following functions. These functions
may be composed according to the [DCP](http://dcp.stanford.edu)
composition rules to form new convex, concave, or affine expressions.
Convex.jl transforms each problem into an equivalent conic program in
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

| operation                                        | description                                                                                                                                           | vexity  | slope                                                                                                     | notes                          |
| ------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------- | ------- | --------------------------------------------------------------------------------------------------------- | ------------------------------ |
| `x+y` or `x.+y`                                  | addition                                                                                                                                              | affine  | increasing                                                                                                | none                           |
| `x-y` or `x.-y`                                  | subtraction                                                                                                                                           | affine  | increasing in $x$ decreasing in $y$                                                               | none none                      |
| `x*y`                                            | multiplication                                                                                                                                        | affine  | increasing if constant term $\ge 0$ decreasing if constant term $\le 0$ not monotonic otherwise | PR: one argument is constant   |
| `x/y`                                            | division                                                                                                                                              | affine  | increasing                                                                                                | PR: $y$ is scalar constant |
| `dot(*)(x, y)`                                   | elementwise multiplication                                                                                                                            | affine  | increasing                                                                                                | PR: one argument is constant   |
| `dot(/)(x, y)`                                   | elementwise division                                                                                                                                  | affine  | increasing                                                                                                | PR: one argument is constant   |
| `x[1:4, 2:3]`                                    | indexing and slicing                                                                                                                                  | affine  | increasing                                                                                                | none                           |
| `diag(x, k)`                                     | $k$-th diagonal of a matrix                                                                                                                       | affine  | increasing                                                                                                | none                           |
| `diagm(x)`                                       | construct diagonal matrix                                                                                                                             | affine  | increasing                                                                                                | PR: $x$ is a vector        |
| `x'`                                             | transpose                                                                                                                                             | affine  | increasing                                                                                                | none                           |
| `vec(x)`                                         | vector representation                                                                                                                                 | affine  | increasing                                                                                                | none                           |
| `dot(x,y)`                                       | $\sum_i x_i y_i$                                                                                                                              | affine  | increasing                                                                                                | PR: one argument is constant   |
| `kron(x,y)`                                      | Kronecker product                                                                                                                                     | affine  | increasing                                                                                                | PR: one argument is constant   |
| `vecdot(x,y)`                                    | `dot(vec(x),vec(y))`                                                                                                                                  | affine  | increasing                                                                                                | PR: one argument is constant   |
| `sum(x)`                                         | $\sum_{ij} x_{ij}$                                                                                                                             | affine  | increasing                                                                                                | none                           |
| `sum(x, k)`                                      | sum elements across dimension $k$                                                                                                                 | affine  | increasing                                                                                                | none                           |
| `sumlargest(x, k)`                               | sum of $k$ largest elements of $x$                                                                                                            | convex  | increasing                                                                                                | none                           |
| `sumsmallest(x, k)`                              | sum of $k$ smallest elements of $x$                                                                                                           | concave | increasing                                                                                                | none                           |
| `dotsort(a, b)`                                  | `dot(sort(a),sort(b))`                                                                                                                                | convex  | increasing                                                                                                | PR: one argument is constant   |
| `reshape(x, m, n)`                               | reshape into $m \times n$                                                                                                                        | affine  | increasing                                                                                                | none                           |
| `minimum(x)`                                     | $\min(x)$                                                                                                                                        | concave | increasing                                                                                                | none                           |
| `maximum(x)`                                     | $\max(x)$                                                                                                                                        | convex  | increasing                                                                                                | none                           |
| `[x y]` or `[x; y]` `hcat(x, y)` or `vcat(x, y)` | stacking                                                                                                                                              | affine  | increasing                                                                                                | none                           |
| `tr(x)`                                          | $\mathrm{tr} \left(X \right)$                                                                                                                  | affine  | increasing                                                                                                | none                           |
| `partialtrace(x,sys,dims)`                       | Partial trace                                                                                                                                         | affine  | increasing                                                                                                | none                           |
| `partialtranspose(x,sys,dims)`                   | Partial transpose                                                                                                                                     | affine  | increasing                                                                                                | none                           |
| `conv(h,x)`                                      | $h \in \mathbb{R}^m$, $x \in \mathbb{R}^n$, $h\star x \in \mathbb{R}^{m+n-1}$; entry $i$ is given by $\sum_{j=1}^m h_jx_{i-j+1}$ with $x_k=0$ for $k$ out of bounds | affine  | increasing if $h\ge 0$ decreasing if $h\le 0$ not monotonic otherwise                           | PR: $h$ is constant        |
| `min(x,y)`                                       | $\min(x,y)$                                                                                                                                      | concave | increasing                                                                                                | none                           |
| `max(x,y)`                                       | $\max(x,y)$                                                                                                                                      | convex  | increasing                                                                                                | none                           |
| `pos(x)`                                         | $\max(x,0)$                                                                                                                                      | convex  | increasing                                                                                                | none                           |
| `neg(x)`                                         | $\max(-x,0)$                                                                                                                                     | convex  | decreasing                                                                                                | none                           |
| `invpos(x)`                                      | $1/x$                                                                                                                                             | convex  | decreasing                                                                                                | IC: $x>0$                 |
| `abs(x)`                                         | $\left\|x\right\|$                                                                                                                              | convex  | increasing on $x \ge 0$ decreasing on $x \le 0$                                                 | none                           |

Second-Order Cone Representable Functions
-----------------------------------------

An optimization problem using these functions can be solved by any SOCP
solver (including ECOS, SCS, Mosek, Gurobi, and CPLEX). Of course, if an
optimization problem has both LP and SOCP representable functions, then
any solver that can solve both LPs and SOCPs can solve the problem.

| operation           | description                                                                         | vexity                              | slope                                                                           | notes                                                                                                                                |
| ------------------- | ----------------------------------------------------------------------------------- | ----------------------------------- | ------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `norm(x, p)`        | $(\sum x_i^p)^{1/p}$                                                          | convex                              | increasing on $x \ge 0$ decreasing on $x \le 0$                       | PR: `p >= 1`                                                                                                                         |
| `vecnorm(x, p)`     | $(\sum x_{ij}^p)^{1/p}$                                                       | convex                              | increasing on $x \ge 0$ decreasing on $x \le 0$                       | PR: `p >= 1`                                                                                                                         |
| `quadform(x, P)`    | $x^T P x$                                                                       | convex in $x$ affine in $P$ | increasing on $x \ge 0$ decreasing on $x \le 0$ increasing in $P$ | PR: either $x$ or $P$ must be constant; if $x$ is not constant, then $P$ must be symmetric and positive semidefinite |
| `quadoverlin(x, y)` | $x^T x/y$                                                                       | convex                              | increasing on $x \ge 0$ decreasing on $x \le 0$ decreasing in $y$ | IC: $y > 0$                                                                                                                     |
| `sumsquares(x)`     | $\sum x_i^2$                                                                  | convex                              | increasing on $x \ge 0$ decreasing on $x \le 0$                       | none                                                                                                                                 |
| `sqrt(x)`           | $\sqrt{x}$                                                                     | concave                             | decreasing                                                                      | IC: $x>0$                                                                                                                       |
| `square(x), x^2`    | $x^2$                                                                           | convex                              | increasing on $x \ge 0$ decreasing on $x \le 0$                       | PR : $x$ is scalar                                                                                                               |
| `dot(^)(x,2)`       | $x.^2$                                                                          | convex                              | increasing on $x \ge 0$ decreasing on $x \le 0$                       | elementwise                                                                                                                          |
| `geomean(x, y)`     | $\sqrt{xy}$                                                                    | concave                             | increasing                                                                      | IC: $x\ge0$, $y\ge0$                                                                                                       |
| `huber(x, M=1)`     | $\begin{cases} x^2 &\|x\| \leq M \\ 2M\|x\| - M^2 &\|x\| > M \end{cases}$ | convex                              | increasing on $x \ge 0$ decreasing on $x \le 0$                       | PR: $M>=1$                                                                                                                      |


Exponential Cone Representable Functions
----------------------------------------

An optimization problem using these functions can be solved by any
exponential cone solver (SCS).

| operation         | description                                | vexity  | slope         | notes          |
| ----------------- | ------------------------------------------ | ------- | ------------- | -------------- |
| `logsumexp(x)`    | $\log(\sum_i \exp(x_i))$          | convex  | increasing    | none           |
| `exp(x)`          | $\exp(x)$                             | convex  | increasing    | none           |
| `log(x)`          | $\log(x)$                             | concave | increasing    | IC: $x>0$ |
| `entropy(x)`      | $\sum_{ij} -x_{ij} \log (x_{ij})$ | concave | not monotonic | IC: $x>0$ |
| `logisticloss(x)` | $\log(1 + \exp(x_i))$               | convex  | increasing    | none           |

Semidefinite Program Representable Functions
--------------------------------------------

An optimization problem using these functions can be solved by any SDP
solver (including SCS and Mosek).

| operation          | description                       | vexity  | slope         | notes                          |
| ------------------ | --------------------------------- | ------- | ------------- | ------------------------------ |
| `nuclearnorm(x)`   | sum of singular values of $x$ | convex  | not monotonic | none                           |
| `operatornorm(x)`  | max of singular values of $x$ | convex  | not monotonic | none                           |
| `eigmax(x)`     | max eigenvalue of $x$         | convex  | not monotonic | none                           |
| `eigmin(x)`     | min eigenvalue of $x$         | concave | not monotonic | none                           |
| `matrixfrac(x, P)` | $x^TP^{-1}x$                  | convex  | not monotonic | IC: P is positive semidefinite |

Exponential + SDP representable Functions
-----------------------------------------

An optimization problem using these functions can be solved by any
solver that supports exponential constraints *and* semidefinite
constraints simultaneously (SCS).

| operation   | description                   | vexity  | slope      | notes                          |
| ----------- | ----------------------------- | ------- | ---------- | ------------------------------ |
| `logdet(x)` | log of determinant of $x$ | concave | increasing | IC: x is positive semidefinite |

Promotions
----------

When an atom or constraint is applied to a scalar and a higher
dimensional variable, the scalars are promoted. For example, we can do
`max(x, 0)` gives an expression with the shape of `x` whose elements are
the maximum of the corresponding element of `x` and `0`.
