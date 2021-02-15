# # Lasso, Ridge and Elastic Net Regressions
#
# This notebook presents a simple implementation of Lasso and elastic net regressions.

#-

# ## Load Packages and Extra Functions

using DelimitedFiles, LinearAlgebra, Statistics, Plots, Convex, SCS

import MathOptInterface
const MOI = MathOptInterface

# # Loading Data
#
# We use the diabetes data from Efron et al, downloaded from https://web.stanford.edu/~hastie/StatLearnSparsity_files/DATA/diabetes.html and then converted from a tab to a comma delimited file.
#
# All data series are standardised (see below) to have zero means and unit standard deviation, which improves the numerical stability. (Efron et al do not standardise the scale of the response variable.)

(x,header) = readdlm("aux_files/diabetes.csv",',',header=true)
#display(header)
#display(x)

x = (x .- mean(x,dims=1))./std(x,dims=1)          #standardise

(Y,X) = (x[:,end],x[:,1:end-1]);                  #to get traditional names
xNames = header[1:end-1];

# # Lasso, Ridge and Elastic Net Regressions
#
# (a)  The regression is $Y = Xb + u$,
# where $Y$ and $u$ are $T \times 1$, $X$ is $T \times K$, and $b$ is the $K$-vector of regression coefficients.
#
# (b) We want to minimize $(Y-Xb)'(Y-Xb) + \gamma \sum |b_i| + \lambda \sum b_i^2$.
#
# (c) We can equally well minimise $b'Qb - 2c'b + \gamma \sum |b_i| + \lambda \sum b_i^2$,
# where $Q = X'X$ and $c=X'Y$
#
# (d) Lasso: $\gamma>0,\lambda=0$; Ridge: $\gamma=0,\lambda>0$; elastic net: $\gamma>0,\lambda>0$.

"""
    LassoEN(Y,X,γ,λ)

Do Lasso (set γ>0,λ=0), ridge (set γ=0,λ>0) or elastic net regression (set γ>0,λ>0).


## Input
- `Y::Vector`:     T-vector with the response (dependent) variable
- `Y::VecOrMat`:   TxK matrix of covariates (regressors)
- `γ::Number`:     penalty on sum(abs.(b))
- `λ::Number`:     penalty on sum(b.^2)

"""
function LassoEN(Y,X,γ,λ=0.0)

    K = size(X,2)

    b_ls = X\Y                    #LS estimate of weights, no restrictions

    Q  = X'X
    c  = X'Y                      #c'b = Y'X*b

    b  = Variable(K)              #define variables to optimize over
    L1 = quadform(b,Q)            #b'Q*b
    L2 = dot(c,b)                 #c'b
    L3 = norm(b,1)                #sum(|b|)
    L4 = sumsquares(b)            #sum(b^2)

    Sol = minimize(L1-2*L2+γ*L3+λ*L4)      #u'u + γ*sum(|b|) + λsum(b^2), where u = Y-Xb
    solve!(Sol,()->SCS.Optimizer(verbose = false))
    Sol.status == MOI.OPTIMAL ? b_i = vec(evaluate(b)) : b_i = NaN

    return b_i, b_ls

end

# The next cell makes a Lasso regression for a single value of γ.


K = size(X,2)
γ = 100

(b,b_ls) = LassoEN(Y,X,γ)

println("OLS and Lasso coeffs (with γ=$γ)")
display([["" "OLS" "Lasso"];xNames b_ls b])

# # Redo the Lasso Regression with Different Gamma Values
#
#
# We now loop over $\gamma$ values.
#
# Remark: it would be quicker to put this loop inside the `LassoEN()` function so as to not recreate `L1`-`L4`.

nγ = 101
γM = range(0; stop=600, length=nγ)             #different γ values

bLasso   = fill(NaN,size(X,2),nγ)       #results for γM[i] are in bLasso[:,i]
for i = 1:nγ
    bLasso[:,i], = LassoEN(Y,X,γM[i])
end

#-

plot(log.(γM),bLasso',
     title = "Lasso regression coefficients",
     xlabel = "log(γ)",
     label = permutedims(xNames),
     size = (600,400))

# # Ridge Regression
#
# We use the same function to do a ridge regression. Alternatively, do `b = inv(X'X + λ*I)*X'Y`.

nλ = 101
λM = range(0; stop=3000, length=nλ)

bRidge  = fill(NaN,size(X,2),nλ)
for i = 1:nλ
    bRidge[:,i], = LassoEN(Y,X,0,λM[i])
end

#-

plot(log.(λM),bRidge',
     title = "Ridge regression coefficients",
     xlabel = "log(λ)",
     label = permutedims(xNames),
     size = (600,400))

# # Elastic Net Regression

λ = 200
println("redo the Lasso regression, but with λ=$λ: an elastic net regression")

bEN   = fill(NaN,size(X,2),nγ)
for i = 1:nγ
    bEN[:,i], = LassoEN(Y,X,γM[i],λ)
end

#-

plot(log.(γM),bEN',
     title = "Elastic Net regression coefficients",
     xlabel = "log(γ)",
     label = permutedims(xNames),
     size = (600,400))
#-
