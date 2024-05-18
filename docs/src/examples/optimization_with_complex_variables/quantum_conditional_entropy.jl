# # Continuity of the quantum conditional entropy
#
# The quantum conditional entropy is given by
#
# ```math
# S(A|B)_\rho := S(\rho^{AB}) - S(\rho^{B})
# ```
#
# where $S$ is the von Neumann entropy,
#
# ```math
# S(\rho) := - \text{tr}(\rho \log \rho)
# ```
#
# and $\rho$ is a positive semidefinite trace-1 matrix (density matrix).
#
# Here, $\rho^{AB}$ represents an operator on the tensor product of two finite-dimensional Hilbert spaces $A$ and $B$ (with dimensions $d_A$ and $d_B$ respectively), so we can regard $\rho_{AB}$ as a matrix on the vector space $\mathbb{C}^{d_Ad_B}$. Moreover, $\rho^B$ denotes the partial trace of $\rho^{AB}$ over the system $A$, so $\rho^B$ is a matrix on $\mathbb{C}^{d_B}$.
#
# One question is how much can $S(A|B)_\rho$ vary between two density matrices $\rho$ and $\sigma$ as a function of the trace-distance $\text{trdist}(\rho, \sigma) := \|\rho-\sigma\|_1 = \frac{1}{2} \text{tr}\left(\sqrt{(\rho-\sigma)^\dagger (\rho-\sigma)}\right)$ (i.e. 1/2 of the nuclear norm). Here the trace distance is meaningful as it is the quantum analog to the total variation distance, and has an inteprepration in terms of the maximal possible probability to distinguish between $\rho$ and $\sigma$ by measurement.
#
# The Alicki-Fannes-Winter (AFW) bound (https://arxiv.org/abs/1507.07775v6 Lemma 2) states that if $\rho$ and $\sigma$ are density matrices, then $\text{trdist}(\rho, \sigma) \leq \varepsilon \leq 1$, then
#
# ```math
# | S(A|B)_\rho - S(A|B)_\sigma| \leq 2 \varepsilon \log d_A + (1 + \varepsilon) h \left(\frac{\varepsilon}{1+\varepsilon}\right)
# ```
#
# where $h(x) = -x\log x  - (1-x)\log(1-x)$ is the binary entropy.
#
# We can illustrate this bound by computing
#
# ```math
#  \max_{\rho}  S(A|B)_\rho - S(A|B)_\sigma
# ```
#
# for a fixed state $\sigma$, and comparing to the AFW bound.
#
# We will choose $d_A=d_B=2$, and $\sigma$ as the maximally entangled state:
#
# ```math
# \sigma = \frac{1}{2}\begin{pmatrix}1 & 0 & 0 & 1\\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 1 & 0 & 0 & 1\\\end{pmatrix}
# ```math
#
# First, we can formulate the conditional entropy in terms of the relative entropy using the relationship
#
# ```math
# S(A|B)_\rho = - D(\rho_^{AB} \| I_A \otimes \rho^B)
# ```math
#
# where $D$ is the quantum relative entropy. Thus:

using Convex
using LinearAlgebra: I

function quantum_conditional_entropy(ρ_AB, d_A, d_B)
    ρ_B = partialtrace(ρ_AB, 1, [d_A, d_B])
    return -quantum_relative_entropy(ρ_AB, kron(I(d_A), ρ_B))
end

# Now we setup the problem data:

ϵ = 0.1
d_A = d_B = 2
σ_AB = 0.5 * [
    1 0 0 1
    0 0 0 0
    0 0 0 0
    1 0 0 1
]

# And we build and solve problem itself

using SCS

ρ_AB = HermitianSemidefinite(d_A * d_B)
add_constraint!(ρ_AB, tr(ρ_AB) == 1)

add_constraint!(ρ_AB, 0.5 * nuclearnorm(ρ_AB - σ_AB) ≤ ϵ)
problem = maximize(quantum_conditional_entropy(ρ_AB, d_A, d_B))
solve!(problem, SCS.Optimizer; silent=true)

# We can then check the observed difference in relative entropies:

difference = evaluate(quantum_conditional_entropy(ρ_AB, d_A, d_B) - quantum_conditional_entropy(σ_AB, d_A, d_B))

# We can compare to the bound:
h(x) = -x*log(x)  - (1-x)*log(1-x)
2 * ϵ *  log(d_A) + (1 + ϵ) * h(ϵ/(1+ϵ))
