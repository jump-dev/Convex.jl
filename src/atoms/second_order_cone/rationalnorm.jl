#############################################################################
# rationalnorm.jl
#
# Handles the k-norms for k > 1 where k is rational. Reduces the
# k-norm constraint to at most 2d ceil(log2(n + m)) + O(d) second
# order cone constraints. Here d is the dimension of the problem and k
# is converted to a rational with value k = n / m. The procedure for
# this is based on the paper "Second-order cone programming" by
# F. Alizadeh and D. Goldfarb, Mathematical Programming, Series B,
# 95:3-51, 2001, and is documented in the pdf available at
# https://github.com/JuliaOpt/Convex.jl/raw/master/docs/supplementary/rational_to_socp.pdf
#############################################################################

export rationalnorm

### k-norm for rational k

struct RationalNormAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    k::Rational{Int64}

    function RationalNormAtom(x::AbstractExpr, k::Rational{Int})
        children = (x,)
        k >= 1 || error("p-norms not defined for p < 1")
        return new(:rationalnorm, hash(children), children, (1,1), k)
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
    return sum(abs.(evaluate(x.children[1])).^x.k)^(1/x.k);
end



# conic_form!(x::RationalNormAtom, unique_conic_forms)
#
# Formulate the conic constraint
#
# sum(abs(x.^k))^(1/k) <= t
#
# into constraints of the form
#
# v[j] >= abs(x[j]) for all j, s >= 0, sum(s) <= t,
# v[j]^n <= s[j]^m * t^(n-m)
#
# where k is the rational n // m. Then, we use the reduction
# delineated in Duchi & Namkoong (2014) to reduce the last constraint
# into SOC constraints.
function conic_form!(x::RationalNormAtom, unique_conic_forms)
    if !has_conic_form(unique_conic_forms, x)
        # Add extra variables and constraints
        d = length(x.children[1]);
        v = Variable(d, 1);
        s = Variable(d, 1);
        t = Variable(1);

        # Adding non-negativity constraints for all variables to problem
        obj = conic_form!(t, unique_conic_forms);
        conic_form!(v >= abs(x.children[1]), unique_conic_forms)
        # conic_form!(t >= 0, unique_conic_forms)
        conic_form!(s >= 0, unique_conic_forms)
        conic_form!(sum(s) <= t, unique_conic_forms)

        # Reduce to SOC constraints (get powers for each element)
        num = Int(numerator(x.k));
        denom = Int(denominator(x.k));
        # Construct list of inequalities of form u^2 <= vw, where the list
        # is given by triples in ineq_list.
        (ineq_list,
         var_list) = psocp.ProductToSimpleInequalities(denom, num - denom);
        if (length(ineq_list) > 10)
            @warn """
                  Rational norm generating $(length(ineq_list)) intermediate constraints.
                  Increasing :max_iters or decreasing solver tolerance may give more accurate solutions.
                  """
        end
        # u corresponds to "introduced" variables; make a matrix of them
        # and then add equality constraints for the first and second
        # variable in the matrix (i.e. first and second columns), as they
        # correspond to s and v.
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
                # Want: t1^2 <= t2 * t3, so add constraints of form
                # norm([2 * t1, t2 - t3], 2) <= t2 + t3.
                conic_form!(SOCElemConstraint(temp2 + temp3,
                                              temp2 - temp3, 2 * temp1),
                            unique_conic_forms);
            end
        end
        cache_conic_form!(unique_conic_forms, x, obj)
    end
    return get_conic_form(unique_conic_forms, x)
end
function rationalnorm(x::AbstractExpr, k::Rational{Int})
    if sign(x) == ComplexSign()
        row,col = size(x)
        if row == 1 || col == 1
            return RationalNormAtom(abs(x),k)
        end
    else
        return RationalNormAtom(x,k)
    end
end
