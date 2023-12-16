using SparseArrays, LinearAlgebra
using Convex: SparseTape, SparseAffineOperation, add_operation
@testset "SparseTape with type $T" for T in (Float64, Float32, BigFloat)
    optimizer = SCS.Optimizer()
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}()),
            optimizer,
        ),
        T,
    )
    d_in = 5
    variables = MOI.add_variables(model, d_in)
    input = rand(T, d_in)
    A = sprand(T, d_in, d_in, 0.1)
    b = sprand(T, d_in, 0.8)
    A_init = copy(A)
    b_init = copy(b)
    op = SparseAffineOperation(A, b)
    tape = SparseTape(op, variables)
    collapsed_tape = SparseAffineOperation(tape)
    @test collapsed_tape.matrix * input + collapsed_tape.vector ≈ A * input + b

    op2 = SparseAffineOperation(sparse(one(T) * I, d_in, d_in), -b)
    tape = add_operation(tape, op2)
    collapsed_tape2 = SparseAffineOperation(tape)
    @test collapsed_tape2.matrix * input + collapsed_tape2.vector ≈ A * input

    op3 = SparseAffineOperation(ones(T, 1, d_in), [zero(T)])
    tape = add_operation(tape, op3)
    collapsed_tape3 = SparseAffineOperation(tape)
    @test collapsed_tape3.matrix * input + collapsed_tape3.vector ≈
          [sum(A * input)]

    @test A_init ≈ A
    @test b_init ≈ b
end
