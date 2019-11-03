# # Fidelity in quantum information theory
# This example is inspired from a lecture of John Watrous in the [course on Theory of Quantum Information](https://cs.uwaterloo.ca/~watrous/CS766/LectureNotes/08.pdf).
#
# The Fidelity between two Hermitian semidefinite matrices P and Q is defined as:
#
# >>>> $F(P,Q) = {||{P}^{1/2}{Q}^{1/2} ||}_{tr}$ = max $|trace({P}^{1/2}U{Q}^{1/2})|$
#
# where the trace norm $||.||_{tr}$ is the sum of the singular values, and the maximization goes over the set of all unitary matrices U. This quantity can be expressed as the optimal value of the following complex-valued SDP:
#
# >>>> **maximize** 1/2 trace(Z+Z*)
#
# >>>> subject to $\left[\begin{array}{cc}P&Z\\{Z}^{*}&Q\end{array}\right] \succeq 0$
#
# >>>> where $Z \in \mathbf {C}^{n \times n}$

using Convex, SCS, LinearAlgebra
if VERSION < v"1.2.0-DEV.0"
     LinearAlgebra.diagm(v::AbstractVector) = diagm(0 => v)
end

n = 20
P = randn(n,n) + im*randn(n,n)
P = P*P'
Q = randn(n,n) + im*randn(n,n)
Q = Q*Q'
Z = ComplexVariable(n,n)
objective = 0.5*real(tr(Z+Z'))
constraint = [P Z;Z' Q] ⪰ 0
problem = maximize(objective,constraint)
solve!(problem, SCSSolver(verbose=0))
computed_fidelity = evaluate(objective)

#-

## Verify that computer fidelity is equal to actual fidelity
P1,P2 = eigen(P)
sqP = P2 * diagm([p1^0.5 for p1 in P1]) * P2'
Q1,Q2 = eigen(Q)
sqQ = Q2 * diagm([q1^0.5 for q1 in Q1]) * Q2'

#-

actual_fidelity = sum(svdvals(sqP * sqQ))

# We can see that the actual fidelity value is very close the computed fidelity value.
