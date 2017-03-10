import Base.sign
export partialtrace
export sign, curvature, monotonicity, evaluate

type PartialTraceAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    sys::Int
    dims::Vector

    function PartialTraceAtom(x::AbstractExpr, sys::Int, dims::Vector)
        if x.size[1] ≠ x.size[2]
            error("Only square matrices are supported")
        end
        if ! (1 ≤ sys ≤ length(dims))
            error("Invalid system, should between 1 and ", length(dims), " got ", sys)
        end
        if x.size[1] ≠ prod(dims)
            error("Dimension of system doesn't correspond to dimension of subsystems")
        end
        children = (x, )
        newsize = (round(Int, x.size[2]/dims[sys]), round(Int, x.size[1]/dims[sys]))
        return new(:partialtrace, hash(children), children, newsize, sys, dims)
    end
end

function sign(x::PartialTraceAtom)
     return sign(x.children[1])
 end

function curvature(x::PartialTraceAtom)
  return ConstVexity()
end

function monotonicity(x::PartialTraceAtom)
    return (Nondecreasing(),)
end

function evaluate(x::PartialTraceAtom)
    ρ = evaluate(x.children[1])

    subsystem = function(sys)
        function term(ρ, j::Int)
            a = speye(1)
            b = speye(1)
            i_sys = 1
            for dim in dims
                if i_sys == sys
                # create a vector that is only 1 at its jth component
                v = spzeros(dim, 1);
                v[j] = 1;
                a = kron(a, v')
                b = kron(b, v)
                else
                    a = kron(a, speye(dim))
                    b = kron(b, speye(dim))
                end
                i_sys += 1
            end
            return a * ρ * b
        end
        return sum([term(ρ, j) for j in 1:dims[sys]])
    end
    sub_systems = [subsystem(i) for i in 1:length(dims)]
end


function conic_form!(x::PartialTraceAtom, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, x)
        sys = x.sys
        dims = x.dims

        # compute the partial trace of x by summing over
        # (I ⊗ <j| ⊗ I) x (I ⊗ |j> ⊗ I) for all j's
        # in the system we want to trace out
        # This function returns every term in the sum
        function term(ρ, j::Int)
            a = speye(1)
            b = speye(1)
            i_sys = 1
            for dim in dims
                if i_sys == sys
                    # create a vector that is only 1 at its jth component
                    v = spzeros(dim, 1);
                    v[j] = 1;
                    a = kron(a, v')
                    b = kron(b, v)
                else
                    a = kron(a, speye(dim))
                    b = kron(b, speye(dim))
                end
                i_sys += 1
            end
            return a * ρ * b
        end

        # sum all terms described above for all j's
        objective = conic_form!(sum([term(x.children[1], j) for j in 1:dims[sys]]), unique_conic_forms)
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

partialtrace(x::AbstractExpr, sys::Int, dim::Vector) = PartialTraceAtom(x, sys, dim)