module TestAtoms

using Convex
using Test

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$name" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_RationalNormAtom_complex_matrix()
    x = Variable(2, 2)
    @test_throws(
        ErrorException("[RationalNormAtom] not defined for complex matrices"),
        rationalnorm(im * x, 3 // 2),
    )
    return
end

function test_RationalNormAtom_matrix()
    x = Variable(2, 2)
    atom = rationalnorm(x, 3 // 2)
    x.value = [1.0 2.0; 3.0 4.0]
    @test evaluate(atom) ≈ sum(abs.(x.value) .^ (3 // 2))^(2 // 3)
    return
end

function test_RationalNormAtom_less_than_1()
    x = Variable(3)
    k = 1 // 2
    @test_throws(
        ErrorException(
            "[RationalNormAtom] p-norms not defined for p < 1. Got $k",
        ),
        rationalnorm(x, k),
    )
    return
end

end  # TestAtoms

TestAtoms.runtests()
