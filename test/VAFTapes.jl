using SparseArrays, LinearAlgebra

@testset "SparseVAFTape with type $T" for T in (Float64, Float32, BigFloat)
    optimizer = SCS.Optimizer()
    model = MOIB.full_bridge_optimizer(MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{T}()),
        optimizer), T)
    d_in = 5
    variables = MOI.add_variables(model, d_in)
    input = rand(T, d_in)
    A = sprand(T, d_in, d_in, .1)
    b = sprand(T, d_in, .8)
    A_init = copy(A)
    b_init = copy(b)
    op = Convex.SparseAffineOperation(A,b)
    tape = Convex.SparseVAFTape([op], variables)
    collapsed_tape = Convex.AffineOperation!(tape)
    @test collapsed_tape.matrix * input + collapsed_tape.vector ≈ A*input + b

    op2 = Convex.SparseAffineOperation(sparse(one(T)*I, d_in, d_in),-b)
    Convex.add_operation!(tape, op2)
    collapsed_tape2 = Convex.AffineOperation!(tape)
    @test collapsed_tape2.matrix * input + collapsed_tape2.vector ≈ A*input

    op3 = Convex.SparseAffineOperation(ones(T, 1, d_in), [zero(T)])
    Convex.add_operation!(tape, op3)
    collapsed_tape3 = Convex.AffineOperation!(tape)
    @test collapsed_tape3.matrix * input + collapsed_tape3.vector ≈ [sum(A*input)]
    
    @test A_init ≈ A
    @test b_init ≈ b
end

@testset "VAFTape with type $T" for T in (Float64, Float32, BigFloat)
    optimizer = SCS.Optimizer()
    model = MOIB.full_bridge_optimizer(MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{T}()),
        optimizer), T)
    d_in = 5
    variables = MOI.add_variables(model, d_in)
    input = rand(T, d_in)
    A = rand(T, d_in, d_in)
    b = rand(T, d_in)
    A_init = copy(A)
    b_init = copy(b)
    op = Convex.AffineOperation(A,b)
    tape = Convex.VAFTape(tuple(op), variables)
    collapsed_tape = Convex.AffineOperation!(tape)
    @test collapsed_tape.matrix * input + collapsed_tape.vector ≈ A*input + b


    tape = Convex.VAFTape(tuple(op), variables)
    op2 = Convex.AffineOperation(one(T)*I, -b)
    tape = Convex.add_operation(tape, op2)
    collapsed_tape2 = Convex.AffineOperation!(tape)
    @test collapsed_tape2.matrix * input + collapsed_tape2.vector ≈ A*input

    tape = Convex.VAFTape(tuple(op), variables)
    tape = Convex.add_operation(tape, op2)
    op3 = Convex.AffineOperation(ones(T, 1, d_in), [zero(T)])
    tape = Convex.add_operation(tape, op3)
    collapsed_tape3 = Convex.AffineOperation!(tape)
    @test collapsed_tape3.matrix * input + collapsed_tape3.vector ≈ [sum(A*input)]
    
    @test A_init ≈ A
    @test b_init ≈ b
end

@testset "MOIU.operate with type $T" for T in (Float32, Float64, BigFloat) begin
    d_in = 5

    optimizer = SCS.Optimizer()
    model = MOIB.full_bridge_optimizer(MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{T}()),
        optimizer), T)
    variables = MOI.add_variables(model, d_in)

    A = rand(T, d_in, d_in)
    b = rand(T, d_in)
    A_init = copy(A)
    b_init = copy(b)
    tape = Convex.VAFTape(tuple(Convex.AffineOperation(A,b)), variables)
    sparse_tape = Convex.SparseVAFTape([Convex.SparseAffineOperation(A,b)], variables)


    B = rand(T, d_in, d_in)
    # MOIU.operate(+, T, B, tape)

    # MOIU.operate(+, T, B, tape)

end
