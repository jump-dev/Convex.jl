# # POVM simulation
# This notebook shows how we can check how much depolarizing noise a qubit positive operator-valued measure (POVM) can take before it becomes simulable by projective measurements. The general method is described in [arXiv:1609.06139](https://arxiv.org/abs/1609.06139). The question of simulability by projective measurements boils down to an SDP problem. Eq. (8) from the paper defines the noisy POVM that we obtain subjecting a POVM $\mathbf{M}$ to a depolarizing channel $\Phi_t$:
#
# $$
# \left[\Phi_t\left(\mathbf{M}\right)\right]_i := t M_i + (1-t)\frac{\mathrm{tr}(M_i)}{d} \mathbb{1}.
# $$
#
# If this visibility $t\in[0,1]$ is one, the POVM $\mathbf{M}$ is simulable.
#
# We will use Convex.jl to solve the SDP problem.

using Convex, SCS, LinearAlgebra
if VERSION < v"1.2.0-DEV.0"
    (I::UniformScaling)(n::Integer) = Diagonal(fill(I.λ, n))
     LinearAlgebra.diagm(v::AbstractVector) = diagm(0 => v)
end

# For the qubit case, a four outcome qubit POVM $\mathbf{M} \in\mathcal{P}(2,4)$ is simulable if and only if 
#
# $M_{1}=N_{12}^{+}+N_{13}^{+}+N_{14}^{+},$
#
# $M_{2}=N_{12}^{-}+N_{23}^{+}+N_{24}^{+},$
#
# $M_{3}=N_{13}^{-}+N_{23}^{-}+N_{34}^{+},$
#
# $M_{4}=N_{14}^{-}+N_{24}^{-}+N_{34}^{-},$
#
# where Hermitian operators $N_{ij}^{\pm}$ satisfy $N_{ij}^{\pm}\geq0$ and $N_{ij}^{+}+N_{ij}^{-}=p_{ij}\mathbb{1}$, where $i<j$ , $i,j=1,2,3,4$ and $p_{ij}\geq0$ as well as $\sum_{i<j}p_{ij}=1$, that is, the $p_{ij}$ values form a probability vector. This forms an SDP feasibility problem, which we can rephrase as an optimization problem by adding depolarizing noise to the left-hand side of the above equations and maximizing the visibility $t$:
#
# $\max_{t\in[0,1]} t$
#
# such that
#
# $t\,M_{1}+(1-t)\,\mathrm{tr}(M_{1})\frac{\mathbb{1}}{2}=N_{12}^{+}+N_{13}^{+}+N_{14}^{+},$
#
# $t\,M_{2}+(1-t)\,\mathrm{tr}(M_{2})\frac{\mathbb{1}}{2}=N_{12}^{-}+N_{23}^{+}+N_{24}^{+},$
#
# $t\,M_{3}+(1-t)\,\mathrm{tr}(M_{3})\frac{\mathbb{1}}{2}=N_{13}^{-}+N_{23}^{-}+N_{34}^{+},$
#
# $t\,M_{4}+(1-t)\,\mathrm{tr}(M_{4})\frac{\mathbb{1}}{2}=N_{14}^{-}+N_{24}^{-}+N_{34}^{-}$.
#
# We organize these constraints in a function that takes a four-output qubit POVM as its argument:

function get_visibility(K)
    noise = real([tr(K[i])*I(2)/2 for i=1:size(K, 1)])
    P = [[ComplexVariable(2, 2) for i=1:2] for j=1:6]
    q = Variable(6, Positive())
    t = Variable(1, Positive())
    constraints = [P[i][j] in :SDP for i=1:6 for j=1:2]
    constraints += sum(q)==1
    constraints += t<=1
    constraints += [P[i][1]+P[i][2] == q[i]*I(2) for i=1:6]
    constraints += t*K[1] + (1-t)*noise[1] == P[1][1] + P[2][1] + P[3][1]
    constraints += t*K[2] + (1-t)*noise[2] == P[1][2] + P[4][1] + P[5][1]
    constraints += t*K[3] + (1-t)*noise[3] == P[2][2] + P[4][2] + P[6][1]
    constraints += t*K[4] + (1-t)*noise[4] == P[3][2] + P[5][2] + P[6][2]
    p = maximize(t, constraints)
    solve!(p, SCSSolver(verbose=0))
    return p.optval
end

# We check this function using the tetrahedron measurement (see Appendix B in [arXiv:quant-ph/0702021](https://arxiv.org/abs/quant-ph/0702021)). This measurement is non-simulable, so we expect a value below one.

function dp(v)
    I(2) + v[1]*[0 1; 1 0] + v[2]*[0 -im; im 0] + v[3]*[1 0; 0 -1]
end
b = [ 1  1  1; 
     -1 -1  1; 
     -1  1 -1;  
      1 -1 -1]/sqrt(3)
M = [dp(b[i, :]) for i=1:size(b,1)]/4;
get_visibility(M)

# This value matches the one [we obtained](https://github.com/peterwittek/ipython-notebooks/blob/master/Simulating_POVMs.ipynb) using [PICOS](http://picos.zib.de/).
