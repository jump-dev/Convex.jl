#############################################################################
# rational_norm.jl
# Handles the k-norms for k > 1 where k is rational. Reduces the k-norm
# constraint to O(dlog(n+m)) SOC constraints where d is the dimension of
# the problem and k = n // m. See Duchi & Namkoong (2014)'s paper
# for details.
#############################################################################

export rational_norm

### k-norm for rational k

type RationalNormAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::(AbstractExpr,)
  size::(Int, Int)
  k::Rational

  function RationalNormAtom(x::AbstractExpr, k::Rational)
    children = (x,)
    k >= 1 || error("p-norms not defined for p < 1")
    return new(:rational_norm, hash(children), children, (1,1), k)
  end
end

function sign(x::RationalNormAtom)
  return Positive()
end

# The monotonicity
function monotonicity(x::RationalNormAtom)
    return (sign(x.children[1]) * Nondecreasing(),)
end

function curvature(x::RationalNormAtom)
  return ConvexVexity()
end

function evaluate(x::RationalNormAtom)
  return sum(abs(evaluate(x.children[1])).^x.k)^(1/x.k);
end

rational_norm(x::AbstractExpr, k::Rational) = RationalNormAtom(x, k::Rational)


# First formulates the constraint sum(abs(x.^k))^(1/k) <= t
# into constraints of the form
# v[ii] >= abs(x[ii]), t >= 0, s[ii] >= 0,
# sum(s) <= t, v[ii]^n <= s[ii]^m*t^(n-m)
# where k = n // m.
# Then, we use the reduction delineated in Duchi & Namkoong (2014)
# to reduce the last constraint into SOC constraints.

function conic_form!(x::RationalNormAtom, unique_conic_forms)
  if !has_conic_form(unique_conic_forms, x)
    # add extra variables and constraints
    d = length(x.children[1]);
    v = Variable(d, 1);
    s = Variable(d, 1);
    t = Variable(1);

    obj = conic_form!(t, unique_conic_forms);
    conic_form!(v >= abs(x.children[1]), unique_conic_forms)
    conic_form!(t >= 0, unique_conic_forms)
    conic_form!(s >= 0, unique_conic_forms)
    conic_form!(sum(s) <= t, unique_conic_forms)

    # reduce to SOC constraints
    numerator = num(x.k);
    denominator = den(x.k);
    (ineq_list, var_list) = psocp.ProductToSimpleInequalities(numerator, numerator - denominator);

    u = Variable(length(var_list), d);
    # v is the first variable
    conic_form!(v == u[var_list[1], :]', unique_conic_forms)
    # s[jj] is the second variable
    conic_form!(s == u[var_list[2], :]', unique_conic_forms)

    for jj = 1:d
        # t is the third variable
        conic_form!(t == u[var_list[3], jj], unique_conic_forms)
        for ii = 1:length(ineq_list)
            temp1 = u[ineq_list[ii].left_var_ind, jj];
            temp2 = u[ineq_list[ii].t_var_ind, jj];
            temp3 = u[ineq_list[ii].s_var_ind, jj];
            conic_form!(SOCElemConstraint(temp2 + temp3, temp2 - temp3, 2 * temp1), unique_conic_forms)
        end
    end
    cache_conic_form!(unique_conic_forms, x, obj)
  end
  return get_conic_form(unique_conic_forms, x)
end
