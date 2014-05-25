export CanonicalConstr, CvxConstr, ==, >=, <=, <, >, +

# Holds a constraint in canonical_form
# We have an array of `coeffs` where these coefficients corresponds to `vars`
# is_eq: the constraint is either == or <=
# is_conic: specifies whether or not is conic
type CanonicalConstr
  coeffs::VecOrMatOrSparse
  vars::Array{Int64, 1}
  constant::Value
  is_eq::Bool
  is_conic::Bool

  function CanonicalConstr(coeffs::VecOrMat, vars::Array{Int64, 1}, constant::Value, is_eq::Bool, is_conic::Bool)
    return new(coeffs, vars, constant, is_eq, is_conic)
  end

  function CanonicalConstr(coeffs::Number, vars::Array{Int64, 1}, constant::Value, is_eq::Bool, is_conic::Bool)
    return new([coeffs], vars, constant, is_eq, is_conic)
  end

  function CanonicalConstr(coeffs::VecOrMat, vars::Int64, constant::Value, is_eq::Bool, is_conic::Bool)
    return new(coeffs, [vars], constant, is_eq, is_conic)
  end

  function CanonicalConstr(coeffs::Number, vars::Int64, constant::Value, is_eq::Bool, is_conic::Bool)
    return new([coeffs], [vars], constant, is_eq, is_conic)
  end
end

# A constraint between two `AbstractCvxExpr`
type CvxConstr
  head::Symbol
  lhs::AbstractCvxExpr
  rhs::AbstractCvxExpr
  vexity::Symbol
  dual_value
  canon_form::Function

  function CvxConstr(head::Symbol, lhs::AbstractCvxExpr, rhs::AbstractCvxExpr)
    # Check vexity
    if head == :(==)
      if lhs.vexity in (:linear, :constant)  && rhs.vexity in (:linear, :constant)
        vexity = :linear
      else
        error("equality constraints between nonlinear expressions are not DCP compliant")
      end
    elseif head == :(<=)
      if lhs.vexity in (:linear, :constant, :convex) && rhs.vexity in (:linear, :constant, :concave)
        vexity = :convex
      else
        error("constraint is not DCP compliant")
      end
    elseif head == :(>=)
      error(">= should have been transformed to <=")
    else
      error("unrecognized comparison $head")
    end

    canon_form = ()->
      begin
        # If both lhs and rhs are constants, we check if the constraint is feasible
        # If it is not feasible, we throw an error. If it is feasible, the
        # canonical_form will be an empty array
        if lhs.vexity == :constant && rhs.vexity == :constant
          if head == :(==)
            if lhs.value != rhs.value
              error("Infeasible problem. $(lhs.value) != $(rhs.value)")
            end
          elseif head == :(<=)
            if !(lhs.value <= rhs.value)
              error("Infeasible problem. $(lhs.value) is not <= $(rhs.value)")
            end
          end
          return CanonicalConstr[]
        # In this case, we have lhs <= constant
        # We promote the size of the rhs if it is 1 x 1
        # If the lhs is 1 x 1, we write this as ones(...) * lhs <= constant
        # to make comparisons of the same size
        elseif rhs.vexity == :constant
          if rhs.size == (1, 1) && lhs.size != (1, 1)
            rhs = Constant(rhs.value * ones(lhs.size...), rhs.sign)
            coeffs = VecOrMatOrSparse[speye(get_vectorized_size(lhs))]
          elseif rhs.size != (1, 1) && lhs.size == (1, 1)
            coeffs = VecOrMatOrSparse[ones(get_vectorized_size(rhs), 1)]
          elseif lhs.size != rhs.size
            error("Can't compare expressions of size $(x.size) and $(y.size)")
          else
            coeffs = VecOrMatOrSparse[speye(get_vectorized_size(lhs))]
          end

          constant = vec(rhs.value)
          canon_constr = CanonicalConstr(coeffs, unique_id(lhs), constant, (head == :(==)), false)
          canon_constr_array = lhs.canon_form()
          push!(canon_constr_array, canon_constr)

        # If we have lhs <= rhs, the canonical form is [I -I] * [lhs rhs]' <= 0
        # Again, we multiply by ones(...) if the size of lhs or rhs is (1, 1)
        else
          if lhs.size == (1, 1) && rhs.size != (1, 1)
            sz = get_vectorized_size(rhs.size)
            coeffs = VecOrMatOrSparse[ones(sz, 1), -speye(sz)]
          elseif lhs.size != (1, 1) && rhs.size == (1, 1)
            sz = get_vectorized_size(lhs.size)
            coeffs = VecOrMatOrSparse[speye(sz), -ones(sz, 1)]
          elseif lhs.size != rhs.size
            error("Can't compare expressions of size $(x.size) and $(y.size)")
          else
            sz = get_vectorized_size(rhs.size)
            coeffs = VecOrMatOrSparse[speye(sz), -speye(sz)]
          end
          vars = [unique_id(lhs); unique_id(rhs)]
          constant = zeros(sz)

          canon_constr = CanonicalConstr(coeffs, vars, constant, (head == :(==)), false)
          canon_constr_array = lhs.canon_form()
          append!(canon_constr_array, rhs.canon_form())
          push!(canon_constr_array, canon_constr)
        end
        return canon_constr_array
      end

    return new(head, lhs, rhs, vexity, nothing, canon_form)
  end
end

==(x::AbstractCvxExpr, y::AbstractCvxExpr) = CvxConstr(:(==), x, y)
>=(x::AbstractCvxExpr, y::AbstractCvxExpr) = CvxConstr(:(<=), y, x)
<=(x::AbstractCvxExpr, y::AbstractCvxExpr) = CvxConstr(:(<=), x, y)
>(x::AbstractCvxExpr, y::AbstractCvxExpr) = >=(x, y)
<(x::AbstractCvxExpr, y::AbstractCvxExpr) = <=(x, y)

==(x::Constant, y::AbstractCvxExpr) = CvxConstr(:(==), y, x)
<=(x::Constant, y::AbstractCvxExpr) = CvxConstr(:(<=), -y, -x)
>=(x::Constant, y::AbstractCvxExpr) = CvxConstr(:(<=), y, x)
>(x::Constant, y::AbstractCvxExpr) = <=(y, x)
<(x::Constant, y::AbstractCvxExpr) = >=(y, x)

==(x::Value, y::AbstractCvxExpr) = CvxConstr(:(==), y, convert(CvxExpr, x))
>=(x::Value, y::AbstractCvxExpr) = CvxConstr(:(<=), y, convert(CvxExpr, x))
<=(x::Value, y::AbstractCvxExpr) = CvxConstr(:(<=), -y, -convert(CvxExpr, x))
>(x::Value, y::AbstractCvxExpr) = <=(y, x)
<(x::Value, y::AbstractCvxExpr) = >=(y, x)

==(x::AbstractCvxExpr, y::Value)= CvxConstr(:(==), x, convert(CvxExpr, y))
>=(x::AbstractCvxExpr, y::Value) = CvxConstr(:(<=), -x, -convert(CvxExpr, y))
<=(x::AbstractCvxExpr, y::Value) = CvxConstr(:(<=), x, convert(CvxExpr, y))
>(x::AbstractCvxExpr, y::Value) = >=(x, y)
<(x::AbstractCvxExpr, y::Value) = <=(x, y)

+(constraints::Array{CvxConstr}, new_constraint::CvxConstr) =
  push!(constraints, new_constraint)
+(constraints::Array{CvxConstr}, new_constraints::Array{CvxConstr}) =
  append!(constraints, new_constraints)
