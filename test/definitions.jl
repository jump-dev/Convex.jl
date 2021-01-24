function _conv(h, x)
    m = length(h)
    n = length(x)
    zero_pad_x(i) = 1 <= i <= n ? x[i] : 0
    [ sum(h[j]*zero_pad_x(i-j+1) for j = 1:m) for i = 1:m+n-1  ]
end

@testset "Conv (issue #364)" begin
    n = 3
    m = 11
    h = rand(m)
    x = rand(n)
    hvar = Variable(m)
    hvar.value = h
    @test evaluate(conv(hvar, x)) ≈ _conv(h, x)
end


@testset "`conj` (issue #416)" begin
    A = [1 1im; -1im 1]
    X = ComplexVariable(2, 2)
    p = minimize(real(tr(conj(X))), [X == A])
    solve!(p, () -> SCS.Optimizer(verbose=1, eps=1e-6))
    @test evaluate(X) ≈ A atol=1e-3
end
