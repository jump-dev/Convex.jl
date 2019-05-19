import LinearAlgebra.isposdef
import Base.in
export SDPConstraint, isposdef, in, ⪰, ⪯

### Positive semidefinite cone constraint

# TODO: Terrible documentation. Please fix.
struct SDPConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    child::AbstractExpr
    size::Tuple{Int, Int}
    dual::ValueOrNothing

    function SDPConstraint(child::AbstractExpr)
        sz = child.size
        if sz[1] != sz[2]
            error("positive semidefinite expressions must be square")
        end
        id_hash = hash((child, :sdp))
        return new(:sdp, id_hash, child, sz, nothing)
    end
end

function vexity(c::SDPConstraint)
    vex = vexity(c.child)
    if vex == AffineVexity() || vex == ConstVexity()
        return AffineVexity()
    else
        return NotDcp()
    end
end

# users specify SDPs as `A in :SDP` where A is an n x n square matrix
# solvers (Mosek and SCS) specify the *lower triangular part of A* is in the SDP cone
# so we need the lower triangular part (n*(n+1)/2 entries) of c.child to be in the SDP cone
# and we need the corresponding upper elements to match the lower elements,
# which we enforce via equality constraints
function conic_form!(c::SDPConstraint, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    # TODO:
        # 1) propagate dual values
        # 2) presolve to eliminate (many) variables --- variables on upper triangular part often don't matter at all
    if !has_conic_form(unique_conic_forms, c)
        n = c.size[1]
        # construct linear indices to pick out the lower triangular part (including diagonal),
        # the upper triangular part (not including diagonal)
        # and the corresponding entries in the lower triangular part, so
        # symmetry => c.child[upperpart]
        # scale off-diagonal elements by sqrt(2)
        rescale = sqrt(2)*tril(ones(n,n))
        rescale[diagind(n, n)] .= 1.0
        diagandlowerpart = findall(!iszero, vec(rescale))
        lowerpart = Vector{Int}(undef, div(n*(n-1),2))
        upperpart = Vector{Int}(undef, div(n*(n-1),2))
        klower = 0
        # diagandlowerpart in column-major order:
        # ie the (1,1), (2,1), ..., (n,1), (2,2), (3,2), ...
        # consider using  and find(triu(ones(3,3)))
        @inbounds for j = 1:n
            for i = j+1:n
                klower += 1
                upperpart[klower] = n*(i-1) + j # (j,i)th element
                lowerpart[klower] = n*(j-1) + i # (i,j)th element
            end
        end
        objective = conic_form!((broadcast(*,rescale,c.child))[diagandlowerpart], unique_conic_forms)
        sdp_constraint = ConicConstr([objective], :SDP, [div(n*(n+1),2)])
        cache_conic_form!(unique_conic_forms, c, sdp_constraint)
        # make sure upper and lower triangular part match in the solution
        # note that this introduces all-zero rows into the constraint matrix
        # if the matrix we require to be PSD has already been constructed to be symmetric
        equality_constraint = conic_form!(c.child[lowerpart] == c.child[upperpart], unique_conic_forms)

        # for computing duals --- won't work b/c
            # 1) sizes are wrong (will be n(n+1)/2, expects n^2), and
            # 2) the "correct" dual is the sum of the sdp dual and the equality constraint dual (so yes, dual of sdp constraint is *not* symmetric)
        # conic_constr_to_constr[sdp_constraint] = c

    end
    return get_conic_form(unique_conic_forms, c)
end

# TODO: Remove isposdef, change tests to use in. Update documentation and notebooks
function isposdef(x::AbstractExpr)
    if sign(x) == ComplexSign()
        SDPConstraint([real(x) -imag(x);imag(x) real(x)])
    else
        SDPConstraint(x)
    end
end

# TODO: Throw error if symbol is invalid.
function in(x::AbstractExpr, y::Symbol)
    if y == :semidefinite || y == :SDP
        if sign(x) == ComplexSign()
            SDPConstraint([real(x) -imag(x);imag(x) real(x)])
        else
            SDPConstraint(x)
        end
    end
end

function ⪰(x::AbstractExpr, y::AbstractExpr)
    if sign(x) == ComplexSign() || sign(y) == ComplexSign()
        SDPConstraint([real(x-y) -imag(x-y);imag(x-y) real(x-y)])
    else
        SDPConstraint(x-y)
    end
end

function ⪯(x::AbstractExpr, y::AbstractExpr)
    if sign(x) == ComplexSign() || sign(y) == ComplexSign()
        SDPConstraint([real(y-x) -imag(y-x);imag(y-x) real(y-x)])
    else
        SDPConstraint(y-x)
    end
end

function ⪰(x::AbstractExpr, y::Value)
    if sign(x) == ComplexSign() || !isreal(y)
        all(y .== 0) ? SDPConstraint([real(x) -imag(x);imag(x) real(x)]) : SDPConstraint([real(x-Constant(y)) -imag(x-Constant(y));imag(x-Constant(y)) real(x-Constant(y))])
    else
        all(y .== 0) ? SDPConstraint(x) : SDPConstraint(x - Constant(y))
    end

end

function ⪰(x::Value, y::AbstractExpr)
    if sign(y) == ComplexSign() || !isreal(x)
        all(x .== 0) ? SDPConstraint([real(-y) -imag(-y);imag(-y) real(-y)]) : SDPConstraint([real(Constant(x)-y) -imag(Constant(x)-y);imag(Constant(x)-y) real(Constant(x)-y)])
    else
        all(x .== 0) ? SDPConstraint(-y) : SDPConstraint(Constant(x) - y)
    end
end

function ⪯(x::Value, y::AbstractExpr)
    if sign(y) == ComplexSign() || !isreal(x)
        all(x .== 0) ? SDPConstraint([real(y) -imag(y);imag(y) real(y)]) : SDPConstraint([real(y-Constant(x)) -imag(y-Constant(x));imag(y-Constant(x)) real(y-Constant(x))])
    else
        all(x .== 0) ? SDPConstraint(y) : SDPConstraint(y - Constant(x))
    end
end

function ⪯(x::AbstractExpr, y::Value)
    if sign(x) == ComplexSign() || !isreal(y)
        all(y .== 0) ? SDPConstraint([real(-x) -imag(-x);imag(-x) real(-x)]) : SDPConstraint([real(Constant(y)-x) -imag(Constant(y)-x);imag(Constant(y)-x) real(Constant(y)-x)])
    else
        all(y .== 0) ? SDPConstraint(-x) : SDPConstraint(Constant(y) - x)
    end
end
