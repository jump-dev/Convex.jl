# # Time Series Analysis
# A time series is a sequence of data points, each associated with a time. In our example, we will work with a time series of daily temperatures in the city of Melbourne, Australia over a period of a few years. Let $x$ be the vector of the time series, and $x_i$ denote the temperature in Melbourne on day $i$. Here is a picture of the time series:

using Gadfly, Convex, SCS
temps = readdlm("melbourne_temps.txt", ',')
n = size(temps)[1]
p = plot(
  x=1:1500, y=temps[1:1500], Geom.line,
  Theme(panel_fill=color("white"))
)

# We can quickly compute the mean of the time series to be $11.2$. If we were to always guess the mean as the temperature of Melbourne on a given day, the RMS error of our guesswork would be $4.1$. We'll try to lower this RMS error by coming up with better ways to model the temperature than guessing the mean.
#
# A simple way to model this time series would be to find a smooth curve that approximates the yearly ups and downs.
# We can represent this model as a vector $s$ where $s_i$ denotes the temperature on the $i$-th day.
# To force this trend to repeat yearly, we simply want
#
# $$s_i = s_{i + 365}$$
#
# for each applicable $i$.
#
# We also want our model to have two more properties:
#
# - The first is that the temperature on each day in our model should be relatively close to the actual temperature of that day.
# - The second is that our model needs to be smooth, so the change in temperature from day to day should be relatively small. The following objective would capture both properties:
#
# $$\sum_{i = 1}^n (s_i - x_i)^2 + \lambda \sum_{i = 2}^n(s_i - s_{i - 1})^2$$
#
# where $\lambda$ is the smoothing parameter. The larger $\lambda$ is, the smoother our model will be.
#
# The following code uses Convex to find and plot the model:

yearly = Variable(n)
eq_constraints = []
for i in 365 + 1 : n
  eq_constraints += yearly[i] == yearly[i - 365]
end

smoothing = 100
smooth_objective = sumsquares(yearly[1 : n - 1] - yearly[2 : n])
problem = minimize(sumsquares(temps - yearly) + smoothing * smooth_objective, eq_constraints)
solve!(problem, SCSSolver(max_iters=5000, verbose=0))
residuals = temps - evaluate(yearly)

## Plot smooth fit
p = plot(
  layer(x=1:1500, y=evaluate(yearly)[1:1500], Geom.line, Theme(default_color=color("red"), line_width=2px)),
  layer(x=1:1500, y=temps[1:1500], Geom.line),
  Theme(panel_fill=color("white"))
)

# We can also plot the residual temperatures, $r$, defined as $r = x - s$.

## Plot residuals for a few days
p = plot(
    x=1:100, y=residuals[1:100], Geom.line,
    Theme(default_color=color("green"), panel_fill=color("white"))
)

# Our smooth model has a RMS error of $2.7$, a significant improvement from just guessing the mean, but we can do better.
#
# We now make the hypothesis that the residual temperature on a given day is some linear combination of the previous $5$ days. Such a model is called autoregressive. We are essentially trying to fit the residuals as a function of other parts of the data itself. We want to find a vector of coefficients $a$ such that
#
# $$\mbox{r}(i) \approx \sum_{j = 1}^5 a_j \mbox{r}(i - j)$$
#
# This can be done by simply minimizing the following sum of squares objective
#
# $$\sum_{i = 6}^n \left(\mbox{r}(i) - \sum_{j = 1}^5 a_j \mbox{r}(i - j)\right)^2$$
#
# The following Convex code solves this problem and plots our autoregressive model against the actual residual temperatures:

## Generate the residuals matrix
ar_len = 5
residuals_mat = residuals[ar_len : n - 1]
for i = 1:ar_len - 1
  residuals_mat = [residuals_mat residuals[ar_len - i : n - i - 1]]
end

## Solve autoregressive problem
ar_coef = Variable(ar_len)
problem = minimize(sumsquares(residuals_mat * ar_coef - residuals[ar_len + 1 : end]))
solve!(problem, SCSSolver(max_iters=5000, verbose=0))

## plot autoregressive fit of daily fluctuations for a few days
ar_range = 1:145
day_range = ar_range + ar_len
p = plot(
  layer(x=day_range, y=residuals[day_range], Geom.line, Theme(default_color=color("green"))),
  layer(x=day_range, y=residuals_mat[ar_range, :] * evaluate(ar_coef), Geom.line, Theme(default_color=color("red"))),
  Theme(panel_fill=color("white"))
)

# Now, we can add our autoregressive model for the residual temperatures to our smooth model to get an better fitting model for the daily temperatures in the city of Melbourne:

total_estimate = evaluate(yearly)
total_estimate[ar_len + 1 : end] += residuals_mat * evaluate(ar_coef)

## plot final fit of data
p = plot(
  layer(x=1:1500, y=total_estimate[1:1500], Geom.line, Theme(default_color=color("red"))),
  layer(x=1:1500, y=temps[1:1500], Geom.line),
  Theme(panel_fill=color("white"))
)

# The RMS error of this final model is $2.3$.

