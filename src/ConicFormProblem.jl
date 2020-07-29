using SparseArrays
using AbstractTrees: children
export conic_form_problem_solve

function template(p::Problem, context)
    obj_problem = template(p.objective, context)
    for c in p.constraints
        add_constraints_to_context(c, context)
    end
    return obj_problem
end

function create_context(p::Problem{T}, optimizer) where {T}
    model = MOIB.full_bridge_optimizer(
        MOIU.CachingOptimizer(
            MOIU.UniversalFallback(MOIU.Model{T}()),
            optimizer
        ),
        T
    )

    return (var_id_to_moi_indices = OrderedDict{UInt64, Vector{MOI.VariableIndex}}(), model=model, T=T,
            id_to_variables = OrderedDict{UInt64, AbstractVariable}())
end

function _isempty(v::MOI.VectorOfVariables)
    isempty(v.variables)
end

function template(a::AbstractVariable, context)
    var_inds = get!(context.var_id_to_moi_indices, a.id_hash) do
        add_variables!(context.model, a::Variable)
    end

    context.id_to_variables[a.id_hash] = a
    return MOI.VectorOfVariables(var_inds)
end

function add_variables!(model, var::AbstractVariable)
    var.id_hash == objectid(:constant) && error("Internal error: constant used as variable")
    if sign(var) == ComplexSign()
        error("complex not implemented")
    else
        return MOI.add_variables(model, length(var))
    end
end

function promote_size(values)
    ds = unique( MOI.output_dimension(v) for v in values if v isa MOI.AbstractFunction)
    d = only(ds)
    values2 = ( v isa Number ? fill(v, d) : v for v in values )
    values2, d
end

function template(A::AdditionAtom, context)
    subproblems = template.(children(A), Ref(context))
    objectives, d = promote_size(subproblems)
    obj = MOIU.operate(+, context.T, objectives...)
    return obj
end

function template(A::NegateAtom, context)
    subobj = template(only(children(A)), context)
    obj = MOIU.operate(-, context.T, subobj)
    return obj
end

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix, v::MOI.VectorOfVariables) where {T}
    @assert size(A,2) == length(v.variables)
    terms = MOI.VectorAffineTerm{T}[]
    add_terms!(terms, A, 1:size(A,1), v.variables)
    return MOI.VectorAffineFunction(terms, zeros(T, size(A, 1)))
end


function MOIU.operate(::typeof(+), ::Type{T}, v1::MOI.VectorOfVariables, v2::MOI.ScalarAffineFunction) where {T}
    MOIU.operate(+, T, scalar_fn(v1), v2)
end

function MOIU.operate(::typeof(+), ::Type{T}, v1::MOI.VectorAffineFunction, v2::MOI.ScalarAffineFunction) where {T}
    MOIU.operate(+, T, scalar_fn(v1), v2)
end

function MOIU.operate!(::typeof(+), ::Type{T}, v1::MOI.VectorAffineFunction, v2::MOI.ScalarAffineFunction) where {T}
    MOIU.operate!(+, T, scalar_fn(v1), v2)
end

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix, vaf::MOI.VectorAffineFunction) where {T}
    @assert size(A,2) == MOI.output_dimension(vaf)

    # AB[i,l] = sum(A[i,j] * B[j, l] for j)
    # @assert iszero(vaf.constant) # for now
    new_constant = A * vaf.constants

    vats = MOI.VectorAffineTerm{T}[]
    for v in vaf.terms
        j = v.output_index
        l_ind = v.scalar_term.variable_index
        Bjl = v.scalar_term.coefficient
        for i in 1:size(A,1)
            new_coefficient = A[i,j] * Bjl
            push!(vats, MOI.VectorAffineTerm{T}(i, MOI.ScalarAffineTerm{T}(new_coefficient, l_ind)))
        end
    end
    return MOI.VectorAffineFunction(vats, new_constant)
end

function MOIU.operate(::typeof(*), ::Type{T}, A::SparseMatrixCSC, vaf::MOI.VectorAffineFunction) where {T}
    @assert size(A,2) == MOI.output_dimension(vaf)
    # @show size(A), typeof(A), nnz(A)
    # @assert iszero(vaf.constant) # for now
    new_constant = A * vaf.constants

    vats = MOI.VectorAffineTerm{T}[]

    rows = rowvals(A)
    vals = nonzeros(A)
    # AB[:,k] = sum(A[:,j] *B[j, k] for j)

    for vat in vaf.terms
        j = vat.output_index
        Bjk = vat.scalar_term.coefficient
        k = vat.scalar_term.variable_index

        for n in nzrange(A, j)
            i = rows[n]
            Aij = vals[n]
            new_vat=MOI.VectorAffineTerm{T}(i, MOI.ScalarAffineTerm{T}(Aij*Bjk, k))
            push!(vats, new_vat)
        end
    end

    return MOI.VectorAffineFunction{T}(vats, new_constant)
    for (j, saf) in enumerate(old_rows)
        for n in nzrange(A, j)
            i = rows[n]
            Aij = vals[n]

        end
    end




    old_rows = MOIU.scalarize(vaf)
    new_rows = MOI.ScalarAffineFunction{T}[]

    # algorithm:
    # AB[i,:] = sum(A[i,j] *B[j, :] for j)
    rows = rowvals(A)
    vals = nonzeros(A)
    for (j, saf) in enumerate(old_rows)
        for n in nzrange(A, j)
            i = rows[n]
            Aij = vals[n]
            new_row = MOIU.operate!(*, T, saf, Aij)
            push!(new_rows, new_row)
        end
    end
    return MOIU.vectorize(new_rows)
    # I, J, V = findnz(A)
    # for (i, j, v) in zip(I, J, V)
        # Aij = v
        # new_row = MOIU.operate!(*, T, old_rows[j], Aij)
        # push!(new_rows, new_row)
    # end
    # return MOIU.vectorize(new_rows)
end



function MOIU.operate(::typeof(*), ::Type{T}, A::Diagonal, vaf::MOI.VectorAffineFunction) where {T}
    @assert size(A,2) == MOI.output_dimension(vaf)
    # @show size(A), typeof(A), nnz(A)
    # @assert iszero(vaf.constant) # for now
    v = diag(A)
    new_constant = v .* vaf.constants

    vats = MOI.VectorAffineTerm{T}[]

    for vat in vaf.terms
        j = vat.output_index
        Bjk = vat.scalar_term.coefficient
        k = vat.scalar_term.variable_index
        new_vat=MOI.VectorAffineTerm{T}(j, MOI.ScalarAffineTerm{T}(v[j]*Bjk, k))
        push!(vats, new_vat)
    end

    return MOI.VectorAffineFunction{T}(vats, new_constant)
end


function template(x::MultiplyAtom, context)
    subproblems = template.(children(x), Ref(context))

    @assert length(children(x)) == length(subproblems) == 2

    if children(x)[1].size == (1, 1) || children(x)[2].size == (1, 1)
        obj = MOIU.operate(*, context.T, subproblems...)
    end

    if size(children(x)[1], 2) > 1 && size(children(x)[2], 2) > 1
        A, B = children(x)
        if vexity(A) == ConstVexity()
            left_mult_by_A = kron(sparse(1.0I, x.size[2], x.size[2]), evaluate(A))
            obj = MOIU.operate(*, context.T, left_mult_by_A, subproblems[2])
        elseif vexity(B) == ConstVexity()
            right_mult_by_B = kron(transpose(evaluate(B)), sparse(1.0I, x.size[1], x.size[1])) 
            obj = MOIU.operate(*, context.T, right_mult_by_B, subproblems[1])
        else
            error()
        end
    end

    return obj
end


function template(x::DotMultiplyAtom, context)
    subproblems = template.(children(x), Ref(context))

    @assert length(children(x)) == length(subproblems) == 2
    s1, s2 = subproblems
    if vexity(x.children[1]) != ConstVexity()
        if vexity(x.children[2]) != ConstVexity()
            error("multiplication of two non-constant expressions is not DCP compliant")
        else
            # make sure first child is the one that's constant
            x.children[1], x.children[2] = x.children[2], x.children[1]
            s1, s2 = s2, s1
        end
    end

    # A = s1
    # vaf = s2
    # @assert vaf isa VectorAffineFunction
    # @assert A isa AbstractMatrix

    # # AB[i,j] = A[i,j] * B[i, j]
    # @assert iszero(vaf.constants) # for now


    # vats = MOI.VectorAffineTerm{T}[]
    # for v in vaf.terms
    #     i = v.output_index
    #     j = v.scalar_term.variable_index
    #     Bij = v.scalar_term.coefficient
    #     new_coefficient = A[i,j]*Bij
    #     push!(vats, MOI.VectorAffineTerm{T}(i, MOI.ScalarAffineTerm{T}(new_coefficient, j)))
    # end
    # return MOI.VectorAffineFunction{T}(vats, copy(vaf.constants))


    # promote the size of the coefficient matrix, so eg
    # 3 .* x
    # works regardless of the size of x
    coeff = evaluate(x.children[1]) .* ones(size(x.children[2]))
    # promote the size of the variable
    # we've previously ensured neither x nor y is 1x1
    # and that the sizes are compatible,
    # so if the sizes aren't equal the smaller one is size 1
    var = x.children[2]
    if size(var, 1) < size(coeff, 1)
        var = ones(size(coeff, 1)) * var
    elseif size(var, 2) < size(coeff, 2)
        var = var * ones(1, size(coeff, 1))
    end

    const_multiplier = Diagonal(vec(coeff))

    obj = MOIU.operate(*, context.T, const_multiplier, s2)

    return obj
end


function template(A::EucNormAtom, context)
    obj = template(only(children(A)), context)

    x = only(children(A))
    d = length(x)


    t_obj = template(Variable(), context)

    t_single_var_obj = MOI.SingleVariable(only(t_obj.variables))

    f = MOIU.operate(vcat, context.T, t_single_var_obj, obj)
    set = MOI.SecondOrderCone(d+1)
    MOI.add_constraint(context.model, f, set)

    return t_obj
end

function template(A::ReshapeAtom, context)
    obj = template(only(children(A)), context)

    return obj
end

function MOIU.operate(::typeof(sum), ::Type{T}, vaf::MOI.VectorAffineFunction) where {T}
    return MOI.ScalarAffineFunction( [vat.scalar_term for vat in vaf.terms], sum(vaf.constants))
end

function MOIU.operate(::typeof(sum), ::Type{T}, v::MOI.VectorOfVariables) where {T}
    return MOI.ScalarAffineFunction(  [ MOI.ScalarAffineTerm{T}(one(T), vi) for vi in v.variables ], zero(T))
end

function template(A::SumAtom, context)
    subobj = template(only(children(A)), context)

    obj = MOIU.operate(sum, context.T, subobj)
    return obj
end

function template(A::AbsAtom, context)
    x = only(A.children)
    t = Variable(size(x))

    t_obj = template(t, context)

    # t_vec = fill(t, size(x))
    add_constraints_to_context(t >= x, context)
    add_constraints_to_context(t >= -x, context)
    return t_obj
end

# function template(I::IndexAtom, subproblems, context)
#     expr = only(subproblems)
#     if I.rows !== nothing && I.cols !== nothing

#     elseif I.inds !== nothing

#     else
#         error("indexing error")
#     end

#     T = context.T
#     return ConicFormProblem(MOI.FEASIBILITY_SENSE, T.(real(C.value)), T.(imag(C.value)))
# end

function template(C::Constant, context)
    T = context.T
    return T.(C.value)
end

scalar_fn(x) = only(MOIU.scalarize(x))

# function scalar_fn(v::MOI.VectorOfVariables)
    # return MOI.SingleVariable(only(v.variables))
# end

# function scalar_fn(v::MOI.VectorAffineFunction{T}) where {T}
    # MOI.output_dimension(v) == 1 || error("not scalar")
    # return MOI.ScalarAffineFunction{T}([ vat.scalar_term for vat in  v.terms], only(v.constants))
# end

scalar_fn(v::MOI.AbstractScalarFunction) = v

function conic_form_problem_solve(p::Problem, optimizer)
    context = create_context(p, optimizer)
    cfp = template(p, context)

    obj = scalar_fn(cfp)
    model = context.model
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.set(model, MOI.ObjectiveSense(), p.head == :maximize ? MOI.MAX_SENSE : MOI.MIN_SENSE)

    MOI.optimize!(model)
    p.model = model
    status = MOI.get(model, MOI.TerminationStatus())
    p.status = status

    if p.status != MOI.OPTIMAL
        @warn "Problem wasn't solved optimally" status
    end

    dual_status = MOI.get(model, MOI.DualStatus())
    primal_status = MOI.get(model, MOI.PrimalStatus())

    var_to_indices = context.var_id_to_moi_indices
    id_to_variables = context.id_to_variables

    if primal_status != MOI.NO_SOLUTION
        for (id, var_indices) in var_to_indices
            var = id_to_variables[id]
            vexity(var) == ConstVexity() && continue
            vectorized_value =  MOI.get(model, MOI.VariablePrimal(), var_indices)
            set_value!(var, unpackvec(vectorized_value, size(var), sign(var) == ComplexSign()))
        end
    else
        for (id, var_indices) in var_to_indices
            var = id_to_variables[id]
            vexity(var) == ConstVexity() && continue
            set_value!(var, nothing)
        end
    end
    return context
end


function add_constraints_to_context(lt::GtConstraint, context)
    lhs = template(lt.lhs, context)
    rhs = template(lt.rhs, context)
    T = context.T
    objectives, d = promote_size( (lhs, rhs) )

    f = MOIU.operate(-, T, objectives...)
    MOI.add_constraint(context.model, f, MOI.Nonnegatives(MOI.output_dimension(f)))

    return nothing
end
