using Gadfly, Convex, SCS
temps = readdlm("melbourne_temps.txt", ',')
n = size(temps)[1]
p = plot(
  x=1:1500, y=temps[1:1500], Geom.line,
  Theme(panel_fill=color("white"))
)

yearly = Variable(n)
eq_constraints = []
for i in 365 + 1 : n
  eq_constraints += yearly[i] == yearly[i - 365]
end

smoothing = 100
smooth_objective = sum_squares(yearly[1 : n - 1] - yearly[2 : n])
problem = minimize(sum_squares(temps - yearly) + smoothing * smooth_objective, eq_constraints)
solve!(problem, SCSSolver(max_iters=5000, verbose=0))
residuals = temps - evaluate(yearly)

# Plot smooth fit
p = plot(
  layer(x=1:1500, y=evaluate(yearly)[1:1500], Geom.line, Theme(default_color=color("red"), line_width=2px)),
  layer(x=1:1500, y=temps[1:1500], Geom.line),
  Theme(panel_fill=color("white"))
)

# Plot residuals for a few days
p = plot(
    x=1:100, y=residuals[1:100], Geom.line,
    Theme(default_color=color("green"), panel_fill=color("white"))
)

# Generate the residuals matrix
ar_len = 5
residuals_mat = residuals[ar_len : n - 1]
for i = 1:ar_len - 1
  residuals_mat = [residuals_mat residuals[ar_len - i : n - i - 1]]
end

# Solve autoregressive problem
ar_coef = Variable(ar_len)
problem = minimize(sum_squares(residuals_mat * ar_coef - residuals[ar_len + 1 : end]))
solve!(problem, SCSSolver(max_iters=5000, verbose=0))

# plot autoregressive fit of daily fluctuations for a few days
ar_range = 1:145
day_range = ar_range + ar_len
p = plot(
  layer(x=day_range, y=residuals[day_range], Geom.line, Theme(default_color=color("green"))),
  layer(x=day_range, y=residuals_mat[ar_range, :] * evaluate(ar_coef), Geom.line, Theme(default_color=color("red"))),
  Theme(panel_fill=color("white"))
)

total_estimate = evaluate(yearly)
total_estimate[ar_len + 1 : end] += residuals_mat * evaluate(ar_coef)

# plot final fit of data
p = plot(
  layer(x=1:1500, y=total_estimate[1:1500], Geom.line, Theme(default_color=color("red"))),
  layer(x=1:1500, y=temps[1:1500], Geom.line),
  Theme(panel_fill=color("white"))
)
