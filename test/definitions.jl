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
