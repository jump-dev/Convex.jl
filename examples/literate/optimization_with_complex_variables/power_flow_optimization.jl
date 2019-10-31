# The data for example is taken from [MATPOWER](http://www.pserc.cornell.edu/matpower/) website. MATPOWER is Matlab package for solving power flow and optimal power flow problems.  

using Convex
using FactCheck
using MAT   #Pkg.add("MAT")
TOL = 1e-2;
input = matopen("Data.mat")
varnames = names(input)
Data = read(input, "inj","Y");

n=size(Data[2],1);
Y=Data[2];
inj=Data[1];
W = ComplexVariable(n,n);
objective = real(sum(diag(W)));
c1 = Constraint[];
for i=2:n
    push!(c1,sum(W[i,:].*(Y[i,:]'))==inj[i]);
end
c2 = W in :SDP
c3 = real(W[1,1])==1.06^2;
push!(c1, c2)
push!(c1, c3)
p = maximize(objective,c1);
solve!(p)
p.optval
#15.125857662600703
evaluate(objective)
#15.1258578588357


output = matopen("Res.mat")
names(output)
outputData = read(output, "Wres");
Wres = outputData
real_diff = real(W.value) - real(Wres);
imag_diff = imag(W.value) - imag(Wres);
@fact real_diff => roughly(zeros(n,n), TOL)
@fact imag_diff => roughly(zeros(n,n), TOL)

real_diff = real(W.value) - (real(W.value))';
imag_sum = imag(W.value) + (imag(W.value))';
@fact real_diff => roughly(zeros(n,n), TOL)
@fact imag_sum => roughly(zeros(n,n), TOL)

