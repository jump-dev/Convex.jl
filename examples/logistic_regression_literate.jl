using DataFrames
using RDatasets
using Convex
using SCS

#-

## we'll use iris data
## predict whether the iris species is versicolor using the sepal length and width and petal length and width
iris = dataset("datasets", "iris")
## outcome variable: +1 for versicolor, -1 otherwise
iris[:Y] = [species == "versicolor" ? 1.0 : -1.0 for species in iris[:Species]]
Y = array(iris[:Y])
## create data matrix with one column for each feature (first column corresponds to offset)
X = [ones(size(iris, 1)) iris[:SepalLength] iris[:SepalWidth] iris[:PetalLength] iris[:PetalWidth]];

#-

## solve the logistic regression problem
n, p = size(X)
beta = Variable(p)
problem = minimize(logisticloss(-Y.*(X*beta)))

solve!(problem, SCSSolver(verbose=false))

#-

## let's see how well the model fits
using Gadfly
perm = Base.Sort.sortperm(vec(X*beta.value))
set_default_plot_size(25cm, 12cm)
plot(layer(x=1:n,y=(Y[perm]+1)/2,Geom.point),layer(x=1:n,y=logistic(X*beta.value)[perm],Geom.line))

