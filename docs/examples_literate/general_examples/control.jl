# # Control

#-

# A simple control problem on a system usually involves a variable $x(t)$
# that denotes the state of the system over time, and a variable $u(t)$ that
# denotes the input into the system over time. Linear constraints are used to
# capture the evolution of the system over time:
#
# $$
# x(t) = Ax(t - 1) + Bu(t), \ \text{for} \ t = 1,\ldots, T,
# $$
#
# where the numerical matrices $A$ and $B$ are called the dynamics and input matrices,
# respectively.
#
# The goal of the control problem is to find a sequence of inputs
# $u(t)$ that will allow the state $x(t)$ to achieve specified values
# at certain times. For example, we can specify initial and final states of the system:
#
# $$
# \begin{aligned}
#   x(0) &= x_i \\
#   x(T) &= x_f
# \end{aligned}
# $$
#
# Additional states between the initial and final states can also be specified. These
# are known as waypoint constraints. Often, the input and state of the system will
# have physical meaning, so we often want to find a sequence inputs that also
# minimizes a least squares objective like the following:
#
# $$
#   \sum_{t = 0}^T \|Fx(t)\|^2_2 + \sum_{t = 1}^T\|Gu(t)\|^2_2,
# $$
#
# where $F$ and $G$ are numerical matrices.
#
# We'll now apply the basic format of the control problem to an example of controlling
# the motion of an object in a fluid over $T$ intervals, each of $h$ seconds.
# The state of the system at time interval $t$ will be given by the position and the velocity of the
# object, denoted $p(t)$ and $v(t)$, while the input will be forces
# applied to the object, denoted by $f(t)$.
# By the basic laws of physics, the relationship between force, velocity, and position
# must satisfy:
#
# $$
#   \begin{aligned}
#     p(t+1) &= p(t) + h v(t) \\
#     v(t+1) &= v(t) + h a(t)
#   \end{aligned}.
# $$
#
# Here, $a(t)$ denotes the acceleration at time $t$, for which we we use
# $a(t) = f(t) / m + g - d v(t)$,
# where $m$, $d$, $g$ are constants for the mass of the object, the drag
# coefficient of the fluid, and the acceleration from gravity, respectively.
#
# Additionally, we have our initial/final position/velocity conditions:
#
# $$
#   \begin{aligned}
#     p(1) &= p_i\\
#     v(1) &= v_i\\
#     p(T+1) &= p_f\\
#     v(T+1) &= 0
#   \end{aligned}
# $$
#
# One reasonable objective to minimize would be
#
# $$
#   \text{objective} = \mu \sum_{t = 1}^{T+1} (v(t))^2 + \sum_{t = 1}^T (f(t))^2
# $$
#
# We would like to keep both the forces small to perhaps save fuel, and keep
# the velocities small for safety concerns.
# Here $\mu$ serves as a parameter to control which part of the objective we
# deem more important, keeping the velocity small or keeping the force small.
#
# The following code builds and solves our control example:

using Convex, SCS, Plots

## Some constraints on our motion
## The object should start from the origin, and end at rest
initial_velocity = [-20; 100]
final_position = [100; 100]

T = 100 # The number of timesteps
h = 0.1 # The time between time intervals
mass = 1 # Mass of object
drag = 0.1 # Drag on object
g = [0, -9.8] # Gravity on object

## Declare the variables we need
position = Variable(2, T)
velocity = Variable(2, T)
force = Variable(2, T - 1)

## Create a problem instance
mu = 1

## Add constraints on our variables
constraints = Constraint[
    position[:, i+1] == position[:, i] + h * velocity[:, i] for i in 1:T-1
]

for i in 1:T-1
    acceleration = force[:, i] / mass + g - drag * velocity[:, i]
    push!(constraints, velocity[:, i+1] == velocity[:, i] + h * acceleration)
end

## Add position constraints
push!(constraints, position[:, 1] == 0)
push!(constraints, position[:, T] == final_position)

## Add velocity constraints
push!(constraints, velocity[:, 1] == initial_velocity)
push!(constraints, velocity[:, T] == 0)

## Solve the problem
problem = minimize(sumsquares(force), constraints)
solve!(problem, SCS.Optimizer; silent_solver = true)

# We can plot the trajectory taken by the object.

pos = evaluate(position)
plot([pos[1, 1]], [pos[2, 1]], st = :scatter, label = "initial point")
plot!([pos[1, T]], [pos[2, T]], st = :scatter, label = "final point")
plot!(pos[1, :], pos[2, :], label = "trajectory")

# We can also see how the magnitude of the force changes over time.

plot(vec(sum(evaluate(force) .^ 2, dims = 1)), label = "force (magnitude)")
