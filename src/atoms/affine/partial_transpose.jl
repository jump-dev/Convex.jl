import Base.sign
export partialtranspose
export sign, curvature, monotonicity, evaluate

type PartialTransposeAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    sys::Int
    dims::Vector

    function PartialTransposeAtom(x::AbstractExpr, sys::Int, dims::Vector)
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
        #newsize = (round(Int, x.size[2]/dims[sys]), round(Int, x.size[1]/dims[sys]))
        return new(:partialtranspose, hash(children), children, x.size, sys, dims)
    end
end

function sign(x::PartialTransposeAtom)
     return sign(x.children[1])
 end

function curvature(x::PartialTransposeAtom)
  return ConstVexity()
end

function monotonicity(x::PartialTransposeAtom)
    return (Nondecreasing(),)
end

function evaluate(x::PartialTransposeAtom)
    ρ = evaluate(x.children[1])
    dims = x.dims

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
    a = eye(1)
    for i in 1:length(dims)
        if i == x.sys
            a = kron(a,transpose(sub_systems[i]))
        else
            a = kron(a,sub_systems[i])
        end
    end
    return a
end


function conic_form!(x::PartialTransposeAtom, unique_conic_forms::UniqueConicForms)
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
        objective = conic_form!(transpose(sum([term(x.children[1], j) for j in 1:dims[sys]])), unique_conic_forms)
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

partialtranspose(x::AbstractExpr, sys::Int, dim::Vector) = PartialTransposeAtom(x, sys, dim)