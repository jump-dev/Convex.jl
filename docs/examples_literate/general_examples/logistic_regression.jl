# # Logistic regression
using DataFrames
using Plots
using RDatasets
using Convex
using SCS

#-

## we'll use iris data
## predict whether the iris species is versicolor using the sepal length and width and petal length and width
iris = dataset("datasets", "iris")
## outcome variable: +1 for versicolor, -1 otherwise
Y = [species == "versicolor" ? 1.0 : -1.0 for species in iris.Species]
## create data matrix with one column for each feature (first column corresponds to offset)
X = hcat(ones(size(iris, 1)), iris.SepalLength, iris.SepalWidth, iris.PetalLength, iris.PetalWidth);

#-

## solve the logistic regression problem
n, p = size(X)
beta = Variable(p)
problem = minimize(logisticloss(-Y.*(X*beta)))

solve!(problem, SCSSolver(verbose=false))

#-

## let's see how well the model fits
using Plots
logistic(x::Real) = inv(exp(-x) + one(x))
perm = sortperm(vec(X*beta.value))
plot(1:n, (Y[perm] .+ 1)/2, st=:scatter)
plot!(1:n, logistic.(X*beta.value)[perm])
