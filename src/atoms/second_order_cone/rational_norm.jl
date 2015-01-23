######

######

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
    k >= 1 || error("pnorms not implemented for p < 1")
    return new(:rational_norm, hash(children), children, (1,1), k)
  end
end

function sign(x::RationalNormAtom)
  return Positive()
end

# The monotonicity
function monotonicity(x::RationalNormAtom)
  # check if arg is positive or negative, otherwise
  return (NoMonotonicity(),)
end

function curvature(x::RationalNormAtom)
  return ConvexVexity()
end

function evaluate(x::RationalNormAtom)
  return sum(evaluate(x.children[1]).^x.k)^(1/x.k);
end

rational_norm(x::AbstractExpr) = RationalNormAtom(x)

function conic_form!(x::RationalNormAtom, unique_conic_forms)
  if !has_conic_form(unique_conic_forms, x)

    d = x.children[1].size[1];
    v = Variables(d, 1);
    s = Variables(d, 1);
    t = Variables(1);

    p = minimize(t);
    p.constraints += (v >= abs(x.children[1]));
    p.constraints += (t >= 0);
    p.constraints += (s >= 0);
    p.constraints += (sum(s) <= t);

    n = num(k);
    m = den(k);

    [ineq_list, var_list] = PSOCP.ProductToSimpleInequalities(n,m);
    u = Variables(length(var_list), d);
      for jj = 1:d
        for ii = 1:length(ineq_list)
            x = u[ineq_list[ii].left_var_ind, jj];
            t = u[ineq_list[ii].t_var_ind, jj];
            s = u[ineq_list[ii].s_var_ind, jj];
            conic_form!(SOCElemConstraint(t - s, t + s, 2 * x)), unique_conic_forms)
        end
    end
    cache_conic_form!(unique_conic_forms, x, p)
  end
  return get_conic_form(unique_conic_forms, x)
end
