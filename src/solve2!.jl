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
    d = only(unique(MOI.output_dimension(v)
                for v in values
                if v isa MOI.AbstractFunction || v isa VectorAffineFunctionAsMatrix || v isa VAFTape || v isa SparseVAFTape|| v isa SparseVAFTape2 || v isa SparseVAFTape3))
    return (v isa Number ? fill(v, d) : v for v in values)
end

scalar_fn(x::Number) = x # for `satisfy` problems? Not sure...
scalar_fn(x) = only(MOIU.scalarize(x))
scalar_fn(x::VAFTapes) = scalar_fn(to_vaf(x))
scalar_fn(v::MOI.AbstractScalarFunction) = v

function solve2!(problem::Problem{T}, optimizer; kwargs...) where {T}
    if Base.applicable(optimizer)
        return solve2!(problem, optimizer(); kwargs...)
    else
        throw(ArgumentError("MathProgBase solvers like `solve!(problem, SCSSolver())` are no longer supported. Use instead e.g. `solve!(problem, SCS.Optimizer)`."))
    end
end

function solve2!(p::Problem{T}, optimizer::MOI.ModelLike) where {T}
    global context = Context{T}(optimizer)
    cfp = template(p, context)

    model = context.model

    if p.head == :satisfy
        MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    else
        obj = scalar_fn(cfp)
        MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
        MOI.set(model, MOI.ObjectiveSense(),
                p.head == :maximize ? MOI.MAX_SENSE : MOI.MIN_SENSE)
    end
    

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
    

    if dual_status != MOI.NO_SOLUTION
        for (constr, MOI_constr_indices) in pairs(context.constr_to_moi_inds)
            populate_dual!(model, constr, MOI_constr_indices)
        end
    else
        populate_dual!(model, constr, nothing)

    end


    return nothing
end

function populate_dual!(model::MOI.ModelLike, constr::Constraint, ::Nothing)
    constr.dual = nothing
    return nothing
end
