# # Phase recovery using MaxCut
# 
# In this example, we relax the phase retrieval problem similar to the classical
# [MaxCut](http://www-math.mit.edu/~goemans/PAPERS/maxcut-jacm.pdf) semidefinite
# program and recover the phase of the signal given the magnitude of the linear
# measurements.
#
# Phase recovery has wide applications such as  in X-ray and crystallography
# imaging, diffraction imaging or microscopy and audio signal processing. In all
# these applications, the detectors cannot measure the phase of the incoming wave
# and only record its amplitude i.e complex measurements of a signal
# $x \in \mathbb{C}^p$ are obtained from a linear injective operator $A$, **but we
# can only measure the magnitude vector $Ax$, not the phase of $Ax$**.
#
# Recovering the phase of $Ax$ from $|Ax|$ is a **nonconvex optimization problem**. Using results from [this paper](https://arxiv.org/abs/1206.0102), the problem can be relaxed to a (complex) semidefinite program (complex SDP).
#
# The original reprsentation of the problem is as follows:
#
# $$
# \begin{array}{ll}
#   \text{find} & x \in \mathbb{C}^p \\
#     \text{subject to} & |Ax| = b
# \end{array}
# $$
#
# where $A \in \mathbb{C}^{n \times p}$ and $b \in \mathbb{R}^n$.

# In this example, **the problem is to find the phase of $Ax$ given the value $|Ax|$**.
# Given a linear operator $A$ and a vector $b= |Ax|$ of measured amplitudes,
# in the noiseless case, we can write $Ax = \text{diag}(b)u$ where
# $u \in \mathbb{C}^n$  is a phase vector, satisfying
# $|\mathbb{u}_i| = 1$ for $i = 1,\ldots, n$. 
#
# We relax this problem as Complex Semidefinite Programming.
#
# ### Relaxed Problem similar to [MaxCut](http://www-math.mit.edu/~goemans/PAPERS/maxcut-jacm.pdf)
#
# Define the positive semidefinite hermitian matrix
# $M = \text{diag}(b) (I - A A^*) \text{diag}(b)$. The problem is:
#
# $$
# \begin{array}{ll}
#   \text{minimize} & \langle U, M \rangle \\
#     \text{subject to} & \text{diag}(U) = 1\\
#     & U \succeq 0
# \end{array}
# $$
#
# Here the variable $U$ must be hermitian ($U \in \mathbb{H}_n $),
# and we have a solution to the phase recovery problem if $U = u u^*$
# has rank one. Otherwise, the leading singular vector of $U$ can be used
# to approximate the solution.

using Convex, SCS, LinearAlgebra
if VERSION < v"1.2.0-DEV.0"
    (I::UniformScaling)(n::Integer) = Diagonal(fill(I.λ, n))
     LinearAlgebra.diagm(v::AbstractVector) = diagm(0 => v)
end

n = 20
p = 2
A = rand(n,p) + im*randn(n,p)
x = rand(p) + im*randn(p)
b = abs.(A*x) + rand(n)

M = diagm(b)*(I(n)-A*A')*diagm(b)
U = ComplexVariable(n,n)
objective = inner_product(U,M)
c1 = diag(U) == 1 
c2 = U in :SDP
p = minimize(objective,c1,c2)
solve!(p, () -> SCS.Optimizer(verbose=0))
evaluate(U)


#-

# Verify if the rank of $U$ is 1:
B, C = eigen(evaluate(U));
length([e for e in B if(abs(real(e))>1e-4)])

#- 

# Decompose $U = uu^*$ where $u$ is the phase of $Ax$
u = C[:,1];
for i in 1:n
    u[i] = u[i]/abs(u[i])
end
u
