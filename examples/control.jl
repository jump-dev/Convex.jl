# Simple control problem
# Code was initially written by Jenny Hong for EE103
# Translated to lsqpy by Keegan Go
# Written in CVX.jl by Karanveer Mohan and David Zeng
#
# In this control problem, the object starts from the origin

# Some constraints on our motion
# The object should start from the origin, and end at rest
initial_velocity = [-20; 20]
final_position = [10; 0]

T = 100 # The number of timesteps
h = 0.1 # The time between time intervals
mass = 1 # Mass of object
drag = 0.01 # Drag on object


# Declare the variables we need
position = Variable(2, T)
velocity = Variable(2, T)
force = Variable(2, T - 1)

# CVX.jl is yet to support x[:, idx], so we use all_rows instead here
all_rows = 1:2

# Create the list of constraints on our variables
constraints = CvxConstr[]
for i in 1:T - 1
  constraints += position[all_rows, i + 1] == position[all_rows, i] + h * velocity[all_rows, i]
end

for i in 1:T - 1
  constraints += velocity[all_rows, i + 1] == velocity[all_rows, i] + h / mass *
      force[all_rows, i] - drag * velocity[all_rows, i]
end

# Add position constraints
constraints += position[all_rows, 1] == 0
constraints += position[all_rows, T] == final_position

# Add velocity constraints
constraints += velocity[all_rows, 1] == initial_velocity
constraints += velocity[all_rows, T] == 0

# Solve the problem
mu = 1
p = minimize(mu * sum(square(velocity)) + sum(square(force)), constraints)
solve!(p)

# TODO: Test plotting
using PyPlot
figure(0, (4,4))
quiver(position.value[0,0:T-1:2],position.value[1,0:T-1:2],
  force.value[0,::2],force.value[1,::2])
plot(position.value[0,:],position.value[1,:],'r')
plot(0,0,'bo')
plot(final_position[0],final_position[1],'bo')
xlim([-2,12])
ylim([-1,4])
show()
