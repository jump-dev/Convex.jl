# # Time Series Analysis
# A time series is a sequence of data points, each associated with a time. In our example, we will work with a time series of daily temperatures in the city of Melbourne, Australia over a period of a few years. Let $x$ be the vector of the time series, and $x_i$ denote the temperature in Melbourne on day $i$. Here is a picture of the time series:

using Plots, Convex, ECOS, DelimitedFiles
aux(str) = joinpath(@__DIR__, "aux_files", str) # path to auxiliary files

temps = readdlm(aux("melbourne_temps.txt"), ',')
n = size(temps, 1)
plot(
    1:n,
    temps[1:n],
    ylabel = "Temperature (°C)",
    label = "data",
    xlabel = "Time (days)",
    xticks = 0:365:n,
)

# We can quickly compute the mean of the time series to be $11.2$. If we were to always guess the mean as the temperature of Melbourne on a given day, the RMS error of our guesswork would be $4.1$. We'll try to lower this RMS error by coming up with better ways to model the temperature than guessing the mean.
#
# A simple way to model this time series would be to find a smooth curve that approximates the yearly ups and downs.
# We can represent this model as a vector $s$ where $s_i$ denotes the temperature on the $i$-th day.
# To force this trend to repeat yearly, we simply want
#
# $$
#  s_i = s_{i + 365}
# $$
#
# for each applicable $i$.
#
# We also want our model to have two more properties:
#
# - The first is that the temperature on each day in our model should be relatively close to the actual temperature of that day.
# - The second is that our model needs to be smooth, so the change in temperature from day to day should be relatively small. The following objective would capture both properties:
#
# $$
#  \sum_{i = 1}^n (s_i - x_i)^2 + \lambda \sum_{i = 2}^n(s_i - s_{i - 1})^2
# $$
#
# where $\lambda$ is the smoothing parameter. The larger $\lambda$ is, the smoother our model will be.
#
# The following code uses Convex to find and plot the model:

yearly = Variable(n)
eq_constraints = [yearly[i] == yearly[i-365] for i in 365+1:n]

smoothing = 100
smooth_objective = sumsquares(yearly[1:n-1] - yearly[2:n])
problem = minimize(
    sumsquares(temps - yearly) + smoothing * smooth_objective,
    eq_constraints,
);
solve!(
    problem,
    MOI.OptimizerWithAttributes(ECOS.Optimizer, "maxit" => 200, "verbose" => 0),
)
residuals = temps - evaluate(yearly)

## Plot smooth fit
plot(1:n, temps[1:n], label = "data")
plot!(
    1:n,
    evaluate(yearly)[1:n],
    linewidth = 2,
    label = "smooth fit",
    ylabel = "Temperature (°C)",
    xticks = 0:365:n,
    xlabel = "Time (days)",
)

# We can also plot the residual temperatures, $r$, defined as $r = x - s$.

## Plot residuals for a few days
plot(1:100, residuals[1:100], ylabel = "Residuals", xlabel = "Time (days)")

#-

root_mean_square_error = sqrt(sum(x -> x^2, residuals) / length(residuals))
# Our smooth model has a RMS error of $2.7$, a significant improvement from just guessing the mean, but we can do better.
#
# We now make the hypothesis that the residual temperature on a given day is some linear combination of the previous $5$ days. Such a model is called autoregressive. We are essentially trying to fit the residuals as a function of other parts of the data itself. We want to find a vector of coefficients $a$ such that
#
# $$
#  \text{r}(i) \approx \sum_{j = 1}^5 a_j \text{r}(i - j)
# $$
#
# This can be done by simply minimizing the following sum of squares objective
#
# $$
#  \sum_{i = 6}^n \left(\text{r}(i) - \sum_{j = 1}^5 a_j \text{r}(i - j)\right)^2
# $$
#
# The following Convex code solves this problem and plots our autoregressive model against the actual residual temperatures:

## Generate the residuals matrix
ar_len = 5

residuals_mat = Matrix{Float64}(undef, length(residuals) - ar_len, ar_len)
for i in 1:ar_len
    residuals_mat[:, i] = residuals[ar_len-i+1:n-i]
end

## Solve autoregressive problem
ar_coef = Variable(ar_len)
problem =
    minimize(sumsquares(residuals_mat * ar_coef - residuals[ar_len+1:end]))
solve!(
    problem,
    MOI.OptimizerWithAttributes(ECOS.Optimizer, "maxit" => 200, "verbose" => 0),
)

## plot autoregressive fit of daily fluctuations for a few days
ar_range = 1:145
day_range = ar_range .+ ar_len
plot(
    day_range,
    residuals[day_range],
    label = "fluctuations from smooth fit",
    ylabel = "Temperature difference (°C)",
)
plot!(
    day_range,
    residuals_mat[ar_range, :] * evaluate(ar_coef),
    label = "autoregressive estimate",
    xlabel = "Time (days)",
)

# Now, we can add our autoregressive model for the residual temperatures to our smooth model to get an better fitting model for the daily temperatures in the city of Melbourne:

total_estimate = evaluate(yearly)
total_estimate[ar_len+1:end] += residuals_mat * evaluate(ar_coef)

# We can plot the final fit of data across the whole time range:
plot(1:n, temps, label = "data", ylabel = "Temperature (°C)")
plot!(
    1:n,
    total_estimate,
    label = "estimate",
    xticks = 0:365:n,
    xlabel = "Time (days)",
)

# The RMS error of this final model is $\sim 2.3$:
root_mean_square_error =
    sqrt(sum(x -> x^2, total_estimate - temps) / length(temps))
