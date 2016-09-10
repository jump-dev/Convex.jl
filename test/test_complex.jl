using Convex
using FactCheck

TOL = 1e-3

facts("Optimization with complex variables") do

    context("Real Variables with complex equality constraints") do
    n = 10 # variable dimension (parameter)
    m = 5 # number of constraints (parameter)
    xo = rand(n)
    A = randn(m,n) + im*randn(m,n)
    b = A * xo
    x = Variable(n)
    p1 = minimize(sum(x), A*x == b, x>=0)
    solve!(p1)
    x1 = x.value
    
    p2 = minimize(sum(x), real(A)*x == real(b), imag(A)*x==imag(b), x>=0)
    solve!(p2)
    x2 = x.value
    @fact x1==x2 => true
  end

  context("Complex Variable with complex equality constraints") do 
    n = 10 # variable dimension (parameter)
    m = 5 # number of constraints (parameter)
    xo = rand(n)+im*rand(n)
    A = randn(m,n) + im*randn(m,n)
    b = A * xo
    x = ComplexVariable(n)
    p1 = minimize(real(sum(x)), A*x == b, real(x)>=0, imag(x)>=0)
    solve!(p1)
    x1 = x.value
    
    xr = Variable(n)
    xi = Variable(n)
    p2 = minimize(sum(xr), real(A)*xr-imag(A)*xi == real(b), imag(A)*xr+real(A)*xi == imag(b), xr>=0, xi>=0)
    solve!(p2)
    #x2 = xr.value + im*xi.value
    real_diff = real(x1) - xr.value

    @fact real_diff => roughly(zeros(10,1), TOL)
    imag_diff = imag(x1) - xi.value 
    @fact imag_diff => roughly(zeros(10,1), TOL)
    #@fact x1==x2 => true
  end


  context("norm2 atom") do
    a = 2+4im
    x = ComplexVariable()
    objective = norm2(a-x)
    c1 = real(x)>=0
    p = minimize(objective,c1)
    solve!(p)
    @fact p.optval => roughly(0, TOL)
    @fact evaluate(objective) => roughly(0, TOL)
    real_diff = real(x.value) - real(a);
    imag_diff = imag(x.value) - imag(a);
    @fact real_diff => roughly(0, TOL)
    @fact imag_diff => roughly(0, TOL)
    end

    context("sumsquares atom") do
    a = [2+4im;4+6im]
    x = ComplexVariable(2)
    objective = sumsquares(a-x)
    c1 = real(x)>=0
    p = minimize(objective,c1)
    solve!(p)
    @fact p.optval => roughly(0, TOL)
    @fact evaluate(objective) => roughly(zeros(1,1), TOL)
    real_diff = real(x.value) - real(a);
    imag_diff = imag(x.value) - imag(a);
    @fact real_diff => roughly(zeros(2,1), TOL)
    @fact imag_diff => roughly(zeros(2,1), TOL)
    end

    context("abs atom") do
    a = [5-4im]
    x = ComplexVariable()
    objective = abs(a-x)
    c1 = real(x)>=0
    p = minimize(objective,c1)
    solve!(p)
    @fact p.optval => roughly(0, TOL)
    @fact evaluate(objective) => roughly(zeros(1), TOL)
    real_diff = real(x.value) - real(a);
    imag_diff = imag(x.value) - imag(a);
    @fact real_diff => roughly(zeros(1), TOL)
    @fact imag_diff => roughly(zeros(1), TOL)
    end

    context("Complex Semidefinite constraint") do
    n = 10
    A = rand(n,n) + im*rand(n,n)
    A = A + A' # now A is hermitian
    x = ComplexVariable(n,n) 
    objective = sumsquares(A - x)
    c1 = x in :SDP
    p = minimize(objective, c1)
    solve!(p)
    # test that X is approximately equal to posA:
    l,v = eig(A)
    posA = v*diagm(max(l,0))*v'

    real_diff = real(x.value) - real(posA);
    imag_diff = imag(x.value) - imag(posA);
    @fact real_diff => roughly(zeros(n,n), TOL)
    @fact imag_diff => roughly(zeros(n,n), TOL)
    end
end