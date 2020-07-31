using SparseArrays
using AbstractTrees: children
export conic_form_problem_solve

# This is a proof of concept implementation of a different internal structure for Convex.jl
# Currently, just the minimal amount of methods are implemented to solve `testproblem/testproblem.jl`.
# Probably any other problem will error!

# Some design notes

# `template` basically has the same function as `conic_form!`:
# it recurses through the problem, generating a conic reformulation of it.
# it does so in a different way, however. It takes a `context` parameter
# which acts very similarly to how `unique_conic_forms` currently works
# by holding problem-local state. In this case, it holds a MOI model.
# when constructing the conic reformulation, any Convex AbstractVariables
# are associated to MOI variables added to the model. Any constraints
# are simply immediately added to the model. Then a MOI function is returned.
# In the end, we recover an MOI function which is the objective function
# for the problem.

# Some more notes:
# * We get access to both the Convex.jl object (and it's children), as well
#   as the MOI objects returned from calling `template` on the children.
#   We should think about when to use one versus the other.
# * It's crucial that a redesign allows defining atoms in terms of
#   partially specified problems (https://github.com/jump-dev/Convex.jl/issues/310).
# * I've used the convention of not including the `!` when modifying the context,
#   under the assumption that it is often being modified.
# * It might be quite interesting to multithread problem formulation. To do so,
#   we could use the branching tree structure of the problem. One issue is that
#   adding variables to the MOI model is a global operation we only want to do once.

# This is a variant of MOI.VectorAffineFunction which represents
# the transformation `matrix * variables + vector` lazily.
struct VectorAffineFunctionAsMatrix{M,B,V}
    matrix::M
    vector::B
    variables::V
end

function MOI.output_dimension(v::VectorAffineFunctionAsMatrix)
    return size(v.matrix, 1)
end

# convert to a usual VAF
function to_vaf(vaf_as_matrix::VectorAffineFunctionAsMatrix{<:SparseMatrixCSC}, context)
    T = eltype(vaf_as_matrix.matrix)
    I, J, V = findnz(vaf_as_matrix.matrix)
    vats = MOI.VectorAffineTerm{T}[]
    for n in eachindex(I, J, V)
        i = I[n]
        j = J[n]
        v = V[n]
        push!(vats,
              MOI.VectorAffineTerm{T}(i,
                                      MOI.ScalarAffineTerm{T}(v,
                                                              vaf_as_matrix.variables[j])))
    end

    return MOI.VectorAffineFunction{T}(vats, vaf_as_matrix.vector)
end

function to_vaf(vaf_as_matrix::VectorAffineFunctionAsMatrix{<:AbstractMatrix})
    T = eltype(vaf_as_matrix.matrix)
    vats = MOI.VectorAffineTerm{T}[]
    M = vaf_as_matrix.matrix
    for i in 1:size(M, 1)
        for j in 1:size(M, 2)
            # this check turns out to be key for performance on the test problem
            # which means I suspect something is being densely represented when possibly it should be sparse
            iszero(M[i, j]) && continue
            push!(vats,
                  MOI.VectorAffineTerm{T}(i,
                                          MOI.ScalarAffineTerm{T}(M[i, j],
                                                                  vaf_as_matrix.variables[j])))
        end
    end
    return MOI.VectorAffineFunction{T}(vats, vaf_as_matrix.vector)
end

# method for adding constraints and coverting to standard VAFs as needed
MOI_add_constraint(model, f, set) = MOI.add_constraint(model, f, set)

function MOI_add_constraint(model, f::VectorAffineFunctionAsMatrix, set)
    return MOI.add_constraint(model, to_vaf(f), set)
end

# `MOIU.operate` methods for `VectorAffineFunctionAsMatrix`

function MOIU.operate(::typeof(-), ::Type{T},
                      vafasmatrix::VectorAffineFunctionAsMatrix) where {T}
    return VectorAffineFunctionAsMatrix(-vafasmatrix.matrix, -vafasmatrix.vector,
                                        vafasmatrix.variables)
end

function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector,
                      vafasmatrix::VectorAffineFunctionAsMatrix) where {T}
    return VectorAffineFunctionAsMatrix(vafasmatrix.matrix, vafasmatrix.vector + v,
                                        vafasmatrix.variables)
end

function MOIU.operate(::typeof(-), ::Type{T}, vafasmatrix::VectorAffineFunctionAsMatrix,
                      v::AbstractVector) where {T}
    return VectorAffineFunctionAsMatrix(vafasmatrix.matrix, vafasmatrix.vector - v,
                                        vafasmatrix.variables)
end

# A simple type representing a vector of zeros. Maybe should include the size or use FillArrays or similar.
struct Zero end

Base.:(+)(a, z::Zero) = a
Base.:(+)(z::Zero, a) = a
Base.:(-)(z::Zero, a) = -a
Base.:(-)(a, z::Zero) = a
Base.:(-)(z::Zero) = z
Base.:(*)(A, z::Zero) = z

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
                      v::MOI.VectorOfVariables) where {T}
    @assert size(A, 2) == length(v.variables)
    return VectorAffineFunctionAsMatrix(A, Zero(), v.variables)
end

function MOIU.operate(::typeof(+), ::Type{T}, v1::MOI.VectorOfVariables,
                      v2::MOI.ScalarAffineFunction) where {T}
    return MOIU.operate(+, T, scalar_fn(v1), v2)
end

function MOIU.operate(::typeof(+), ::Type{T}, v1::MOI.VectorAffineFunction,
                      v2::MOI.ScalarAffineFunction) where {T}
    return MOIU.operate(+, T, scalar_fn(v1), v2)
end

function MOIU.operate!(::typeof(+), ::Type{T}, v1::MOI.VectorAffineFunction,
                       v2::MOI.ScalarAffineFunction) where {T}
    return MOIU.operate!(+, T, scalar_fn(v1), v2)
end

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
                      vafasmatrix::VectorAffineFunctionAsMatrix) where {T}
    return VectorAffineFunctionAsMatrix(A * vafasmatrix.matrix, A * vafasmatrix.vector,
                                        vafasmatrix.variables)
end

# Other `MOIU.operate` methods

function MOIU.operate(::typeof(sum), ::Type{T}, vaf::MOI.VectorAffineFunction) where {T}
    return MOI.ScalarAffineFunction([vat.scalar_term for vat in vaf.terms],
                                    sum(vaf.constants))
end

function MOIU.operate(::typeof(sum), ::Type{T}, v::MOI.VectorOfVariables) where {T}
    return MOI.ScalarAffineFunction([MOI.ScalarAffineTerm{T}(one(T), vi)
                                     for vi in v.variables], zero(T))
end

function template(p::Problem, context)
    obj_problem = template(p.objective, context)
    for c in p.constraints
        add_constraints_to_context(c, context)
    end
    return obj_problem
end

function create_context(p::Problem{T}, optimizer) where {T}
    model = MOIB.full_bridge_optimizer(MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{T}()),
                                                             optimizer), T)

    return (var_id_to_moi_indices=OrderedDict{UInt64,Vector{MOI.VariableIndex}}(),
            model=model, T=T, id_to_variables=OrderedDict{UInt64,AbstractVariable}())
end

function _isempty(v::MOI.VectorOfVariables)
    return isempty(v.variables)
end

function template(a::AbstractVariable, context)
    var_inds = get!(context.var_id_to_moi_indices, a.id_hash) do
        return add_variables!(context.model, a::Variable)
    end

    context.id_to_variables[a.id_hash] = a
    return MOI.VectorOfVariables(var_inds)
end

function add_variables!(model, var::AbstractVariable)
    var.id_hash == objectid(:constant) && error("Internal error: constant used as variable")
    return if sign(var) == ComplexSign()
        error("complex not implemented")
    else
        return MOI.add_variables(model, length(var))
    end
end

# this is kind of awful; currently needed since sometimes Convex treats numbers as vectors or matrices
# need to think of a better way, possibly via breaking changes to Convex's syntax
function promote_size(values)
    ds = unique(MOI.output_dimension(v)
                for v in values
                if v isa MOI.AbstractFunction || v isa VectorAffineFunctionAsMatrix)
    d = only(ds)
    values2 = (v isa Number ? fill(v, d) : v for v in values)
    return values2, d
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
            right_mult_by_B = kron(transpose(evaluate(B)),
                                   sparse(1.0I, x.size[1], x.size[1]))
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

    # we're going to convert early since we haven't defined `vcat`...
    if obj isa VectorAffineFunctionAsMatrix
        obj = to_vaf(obj)
    end

    f = MOIU.operate(vcat, context.T, t_single_var_obj, obj)
    set = MOI.SecondOrderCone(d + 1)
    MOI_add_constraint(context.model, f, set)

    return t_obj
end

function template(A::ReshapeAtom, context)
    obj = template(only(children(A)), context)

    return obj
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

function template(C::Constant, context)
    T = context.T
    return T.(C.value)
end

scalar_fn(x) = only(MOIU.scalarize(x))

scalar_fn(v::MOI.AbstractScalarFunction) = v

function conic_form_problem_solve(p::Problem, optimizer)
    context = create_context(p, optimizer)
    cfp = template(p, context)

    obj = scalar_fn(cfp)
    model = context.model
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.set(model, MOI.ObjectiveSense(),
            p.head == :maximize ? MOI.MAX_SENSE : MOI.MIN_SENSE)

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
            vectorized_value = MOI.get(model, MOI.VariablePrimal(), var_indices)
            set_value!(var,
                       unpackvec(vectorized_value, size(var), sign(var) == ComplexSign()))
        end
    else
        for (id, var_indices) in var_to_indices
            var = id_to_variables[id]
            vexity(var) == ConstVexity() && continue
            set_value!(var, nothing)
        end
    end
    return nothing
end

function add_constraints_to_context(lt::GtConstraint, context)
    lhs = template(lt.lhs, context)
    rhs = template(lt.rhs, context)
    T = context.T
    objectives, d = promote_size((lhs, rhs))
    f = MOIU.operate(-, T, objectives...)
    MOI_add_constraint(context.model, f, MOI.Nonnegatives(MOI.output_dimension(f)))
    return nothing
end
