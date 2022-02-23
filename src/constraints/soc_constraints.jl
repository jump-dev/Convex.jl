# TODO: Document this. How is this different from SOCElemConstraint? Why do we need both. How does
# conic form work for SOC constraints.
struct SOCConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    children::Tuple
    dual::ValueOrNothing

    function SOCConstraint(args::AbstractExpr...)
        children = tuple(args...)
        id_hash = hash((children, :soc))
        return new(:soc, id_hash, children, nothing)
    end
end

function conic_form!(c::SOCConstraint, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, c)
        objectives = Vector{ConicObj}(undef, length(c.children))
        @inbounds for iobj in 1:length(c.children)
            objectives[iobj] = conic_form!(c.children[iobj], unique_conic_forms)
        end
        cache_conic_form!(
            unique_conic_forms,
            c,
            ConicConstr(objectives, :SOC, [length(x) for x in c.children]),
        )
    end
    return get_conic_form(unique_conic_forms, c)
end

# For debugging created this constraint
socp(args::AbstractExpr...) = SOCConstraint(args::AbstractExpr...)

struct SOCElemConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    children::Tuple
    dual::ValueOrNothing

    function SOCElemConstraint(args::AbstractExpr...)
        children = tuple(args...)
        id_hash = hash((children, :soc_elem))
        return new(:soc_elem, id_hash, children, nothing)
    end
end

function conic_form!(c::SOCElemConstraint, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, c)
        num_constrs = length(c.children[1])
        num_children = length(c.children)
        conic_constrs = Vector{ConicConstr}(undef, num_constrs)
        objectives = Vector{ConicObj}(undef, num_children)
        @inbounds for iobj in 1:num_children
            objectives[iobj] = conic_form!(c.children[iobj], unique_conic_forms)
        end
        @inbounds for row in 1:num_constrs
            objectives_by_row = [get_row(obj, row) for obj in objectives]
            conic_constrs[row] = ConicConstr(objectives_by_row, :SOC, [1, 1, 1])
        end
        cache_conic_form!(unique_conic_forms, c, conic_constrs)
    end
    return get_conic_form(unique_conic_forms, c)
end
