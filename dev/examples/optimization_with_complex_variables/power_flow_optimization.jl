# # Power flow optimization
# The data for example is taken from [MATPOWER](http://www.pserc.cornell.edu/matpower/) website. MATPOWER is Matlab package for solving power flow and optimal power flow problems.

using Convex, SCS
using Test
using MAT   #Pkg.add("MAT")
aux(str) = joinpath(@__DIR__, "aux_files", str) # path to auxiliary files

TOL = 1e-2;
input = matopen(aux("Data.mat"))
varnames = names(input)
Data = read(input, "inj", "Y");

n = size(Data[2], 1);
Y = Data[2];
inj = Data[1];
W = ComplexVariable(n, n);
objective = real(sum(diag(W)));
c1 = Constraint[];
for i in 2:n
    push!(c1, sum(W[i, :] .* (Y[i, :]')) == inj[i])
end
c2 = W in :SDP
c3 = real(W[1, 1]) == 1.06^2;
push!(c1, c2)
push!(c1, c3)
p = maximize(objective, c1);
solve!(p, SCS.Optimizer; silent_solver = true)
p.optval
#15.125857662600703
evaluate(objective)
#15.1258578588357

output = matopen(joinpath(@__DIR__, "Res.mat"))
names(output)
outputData = read(output, "Wres");
Wres = outputData
real_diff = real(evaluate(W)) - real(Wres);
imag_diff = imag(evaluate(W)) - imag(Wres);
@test real_diff ≈ zeros(n, n) atol = TOL
@test imag_diff ≈ zeros(n, n) atol = TOL

real_diff = real(evaluate(W)) - (real(evaluate(W)))';
imag_sum = imag(evaluate(W)) + (imag(evaluate(W)))';
@test real_diff ≈ zeros(n, n) atol = TOL
@test imag_diff ≈ zeros(n, n) atol = TOL
