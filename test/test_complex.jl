using Convex
using Test

TOL = 1e-3

@testset "Optimization with complex variables" begin

    @testset "Real Variables with complex equality constraints" begin
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
    @test x1 == x2
  end

  @testset "Complex Variable with complex equality constraints" begin
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

    @test isapprox(real_diff, zeros(10, 1), atol=TOL)
    imag_diff = imag(x1) - xi.value
    @test isapprox(imag_diff, zeros(10, 1), atol=TOL)
    #@fact x1==x2 --> true
  end


  @testset "norm2 atom" begin
    a = 2+4im
    x = ComplexVariable()
    objective = norm2(a-x)
    c1 = real(x)>=0
    p = minimize(objective,c1)
    solve!(p)
    @test isapprox(p.optval, 0, atol=TOL)
    @test isapprox(evaluate(objective), 0, atol=TOL)
    real_diff = real(x.value) - real(a)
    imag_diff = imag(x.value) - imag(a)
    @test isapprox(real_diff, 0, atol=TOL)
    @test isapprox(imag_diff, 0, atol=TOL)
    end

    @testset "sumsquares atom" begin
    a = [2+4im;4+6im]
    x = ComplexVariable(2)
    objective = sumsquares(a-x)
    c1 = real(x)>=0
    p = minimize(objective,c1)
    solve!(p)
    @test isapprox(p.optval, 0, atol=TOL)
    @test isapprox(evaluate(objective), zeros(1, 1), atol=TOL)
    real_diff = real.(x.value) - real.(a)
    imag_diff = imag.(x.value) - imag.(a)
    @test isapprox(real_diff, zeros(2, 1), atol=TOL)
    @test isapprox(imag_diff, zeros(2, 1), atol=TOL)
    end

    @testset "abs atom" begin
    a = [5-4im]
    x = ComplexVariable()
    objective = abs(a-x)
    c1 = real(x)>=0
    p = minimize(objective,c1)
    solve!(p)
    @test isapprox(p.optval, 0, atol=TOL)
    @test isapprox(evaluate(objective), zeros(1), atol=TOL)
    real_diff = real(x.value) - real(a)
    imag_diff = imag(x.value) - imag(a)
    @test isapprox(real_diff, zeros(1), atol=TOL)
    @test isapprox(imag_diff, zeros(1), atol=TOL)
    end

    @testset "Complex Semidefinite constraint" begin
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
    posA = v*diagm(max.(l,0))*v'

    real_diff = real.(x.value) - real.(posA)
    imag_diff = imag.(x.value) - imag.(posA)
    @test isapprox(real_diff, zeros(n, n), atol=TOL)
    @test isapprox(imag_diff, zeros(n, n), atol=TOL)
    end
end
