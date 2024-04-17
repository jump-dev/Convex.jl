function set_with_size(::Type{MOI.PositiveSemidefiniteConeSquare}, sz::Tuple{Int,Int})
    if sz[1] != sz[2]
        error("Positive semidefinite expressions must be square")
    end
    return MOI.PositiveSemidefiniteConeSquare(sz[1])
end

head(io::IO, ::MOI.PositiveSemidefiniteConeSquare) = print(io, "sdp")

function vexity(vex, ::MOI.PositiveSemidefiniteConeSquare)
    if !(vex in (AffineVexity(), ConstVexity()))
        return NotDcp()
    end
    return AffineVexity()
end

function is_feasible(x, ::MOI.PositiveSemidefiniteConeSquare, tol)
    if !(x ≈ transpose(x))
        @warn "constant SDP constraint is violated"
        return false
    elseif evaluate(LinearAlgebra.eigmin(c.child)) < -tol
        @warn "constant SDP constraint is violated"
        return false
    end
    return true
end

function LinearAlgebra.isposdef(x::AbstractExpr)
    if iscomplex(x)
        return GenericConstraint{MOI.PositiveSemidefiniteConeSquare}(
            [real(x) -imag(x); imag(x) real(x)],
        )
    end
    return GenericConstraint{MOI.PositiveSemidefiniteConeSquare}(x)
end

⪰(x::AbstractExpr, y::AbstractExpr) = isposdef(x - y)

function ⪰(x::AbstractExpr, y::Value)
    if all(y .== 0)
        return isposdef(x)
    end
    return isposdef(x - constant(y))
end

function ⪰(x::Value, y::AbstractExpr)
    if all(x .== 0)
        return isposdef(-y)
    end
    return isposdef(constant(x) - y)
end

⪯(x::AbstractExpr, y::AbstractExpr) = ⪰(y, x)
⪯(x::Value, y::AbstractExpr) = ⪰(y, x)
⪯(x::AbstractExpr, y::Value) = ⪰(y, x)
