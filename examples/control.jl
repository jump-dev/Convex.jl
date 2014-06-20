# Simple control problem
# Code was initially written by Jenny Hong for EE103
# Translated to lsqpy by Keegan Go
# Written in CVX.jl by Karanveer Mohan and David Zeng
#
# In this control problem, the object starts from the origin

using CVX

# Some constraints on our motion
# The object should start from the origin, and end at rest
initial_velocity = [-20; 20]
final_position = [10; 0]

T = 100 # The number of timesteps
h = 0.1 # The time between time intervals
mass = 0.4 # Mass of object
drag = 0.02 # Drag on object


# Declare the variables we need
position = Variable(2, T)
velocity = Variable(2, T)
force = Variable(2, T - 1)

# Create the list of constraints on our variables
constraints = CvxConstr[]
for i in 1:T - 1
  constraints += position[:, i + 1] == position[:, i] + h * velocity[:, i]
end

for i in 1:T - 1
  constraints += velocity[:, i + 1] == velocity[:, i] + h / mass *
      force[:, i] - drag * velocity[:, i]
end

# Add position constraints
constraints += position[:, 1] == 0
constraints += position[:, T] == final_position

# Add velocity constraints
constraints += velocity[:, 1] == initial_velocity
constraints += velocity[:, T] == 0

# Solve the problem
mu = 1
p = minimize(mu * sum_squares(vec(velocity)) + sum_squares(vec(force)), constraints)
solve!(p)


import PyPlot.plt
plt.plot(position.value[1, 1:2:T]', position.value[2, 1:2:T]', "r-", linewidth=1.5)
plt.quiver(position.value[1, 1:2:T], position.value[2, 1:2:T], force.value[1, 1:2:T-1], force.value[2, 1:2:T-1], width=0.002)
plt.plot(0, 0, "bo", markersize=10)
plt.plot(final_position[1], final_position[2], "bo", markersize=10)
plt.xlim([-15, 16])
plt.ylim([-10, 16])
plt.show()
