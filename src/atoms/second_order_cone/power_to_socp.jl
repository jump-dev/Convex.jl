# power_to_socp
#
# Functions that recursively simplify inequalities of the form
#
# x^n <= t^{n - m} s^m, t >= 0, s >= 0
#
# into a series of equations of the form x^2 <= st, which is
# second-order cone representable as
#
# norm([2x, t - s]) <= t + s, t >= 0, s >= 0.
#
# This reduction is documented in the pdf available at
# https://github.com/JuliaOpt/Convex.jl/raw/master/docs/supplementary/rational_to_socp.pdf

module psocp

mutable struct InequalityExpression
    # Represents an inequality of the form
    #
    # x^n <= t^{p_1} s^{p_2} u^{p_3},
    #
    # where we always have p_1 + p_2 + p_3 == n, and we must have p_3
    # belonging to {0, 1, 2}. This allows a particular form of recursive
    # enumeration of inequalities.

    # The powers p_1, p_2, p_3
    power1::Int
    power2::Int
    power3::Int

    # Index for variable represented by x in the above equation
    left_var_ind::Int
    # Index for variable represented by t
    t_var_ind::Int
    # Index for variable represented by s
    s_var_ind::Int
    # Index for variable represented by u (may be -1 if power3 == 0).
    u_var_ind::Int
end

struct SimpleInequalityExpression
    # Represents an inequality of the form
    #
    # x^2 <= t s,
    #
    # which is representable as the SOCP norm([2 x, t - s]) <= t + s. A
    # list of these is the first of the tuple returned by the method
    # ProductToSimpleInequalities.

    # Index for variable represented by x in the above equation
    left_var_ind::Int
    # Index for variable represented by t
    t_var_ind::Int
    # Index for variable represented by s
    s_var_ind::Int
end

# ProductToSimpleInequalities(first_power, second_power)
#
# Computes a sequence of inequalities of the form x^2 <= st, beginning
# from an initial inequality of the form
#
# x^n <= t^{p1} s^{p2},
#
# where n = p1 + p2 and p1 = first_power, p2 = second_power.
#
# Returns the tuple (inequalities, vars). The first of the tuple (the
# list inequalities) of which is an array of
# SimpleInequalityExpressions that consists of all the variable
# inequalities, the second of the tuple (the array vars) is an array
# of all the variable indices in the inequalities.
function ProductToSimpleInequalities(first_power::Int, second_power::Int)
    # Construct the first InequalityExpression, which is an inequality
    # of the form x^n <= t^p1 s^p2.
    @assert first_power > 0 && second_power > 0;
    n = first_power + second_power;
    var_list = [1, 2, 3];
    init_inequality = InequalityExpression(first_power, second_power, 0,
                                                                                 1, 2, 3, -1);
    return ReducePowers(init_inequality, Array{SimpleInequalityExpression}(undef, 0),
                                            var_list);
end

# ReducePowers(curr_inequality::InequalityExpression,
#              ineq_array::Array{SimpleInequalityExpression, 1},
#              variable_array::Array{Int64, 1})
#
# Takes an initial inequality (curr_inequality) representing an
# expression of the form
#
# x^n <= t^{p1} s^{p2} u^{p3},
#
# where p3 is one of 0, 1, or 2, and n = p1 + p2 + p3, and recursively
# reduces it (see the file method-notes.tex) until the inequality
# array is fully populated by inequalities of the form x^2 <= t s,
# which are SOCP representable.
#
# The inequality array is populated (should be initially called with
# an empty array of SimpleInequalityExpressions), as is the
# variable_array, and returns the pair (ineq_array, variable_array)
# once both are fully populated.
function ReducePowers(curr_inequality::InequalityExpression,
                      ineq_array::Array{SimpleInequalityExpression, 1},
                      variable_array::Array{Int, 1})
    ReArrangePowers!(curr_inequality);
    p1 = curr_inequality.power1;
    p2 = curr_inequality.power2;
    p3 = curr_inequality.power3;
    n = (p1 + p2 + p3);
    @assert p1 >= 1;  # , "Must have at least 1 on first power");
    @assert p2 >= 1;  # , "Must have at least 1 on second power");
    # Evaluate cases for variables
    if (p3 == 0)
        # Double check if we have p1 == 1, p2 == 1
        if (p1 == 1 && p2 == 1)
            new_ineq = SimpleInequalityExpression(curr_inequality.left_var_ind,
                                                  curr_inequality.t_var_ind,
                                                  curr_inequality.s_var_ind);
            return (vcat(ineq_array, new_ineq), variable_array);
        end
        return ReduceThirdZero(curr_inequality,
                               ineq_array, variable_array);
    elseif (p3 == 1)
        return ReduceThirdOne(curr_inequality, ineq_array, variable_array);
    elseif (p3 == 2)
        return ReduceThirdTwo(curr_inequality, ineq_array, variable_array);
    else
        error("third power in inequality (", curr_inequality.power3, ") must be 0, 1, or 2");
    end
end

# ReArrangePowers!(current_inequality::InequalityExpression)
#
# Function that takes in the current inequality, then re-orders
# variables so that the variable with smallest power is stored as
# u_var_ind. Also divides each power by their gcd.
function ReArrangePowers!(inequality::InequalityExpression)
    p1 = inequality.power1;
    p2 = inequality.power2;
    p3 = inequality.power3;
    divisor = gcd(p1, p2, p3);
    if (divisor > 1)
        inequality.power1 /= divisor;
        inequality.power2 /= divisor;
        inequality.power3 /= divisor;
        p1 = inequality.power1;
        p2 = inequality.power2;
        p3 = inequality.power3;
    end
    if (p1 <= min(p2, p3))
        # Swap t and u
        inequality.power3 = p1;
        inequality.power1 = p3;
        ind1 = inequality.t_var_ind;
        ind3 = inequality.u_var_ind;
        inequality.t_var_ind = ind3;
        inequality.u_var_ind = ind1;
    elseif (p2 <= min(p1, p3))
        # Swap s and u
        inequality.power3 = p2;
        inequality.power2 = p3;
        ind2 = inequality.s_var_ind;
        ind3 = inequality.u_var_ind;
        inequality.s_var_ind = ind3;
        inequality.u_var_ind = ind2;
    end  # In this case, p3 <= min(p1, p2), so no swap necessary
end

# ReduceThirdZero(curr_inequality::InequalityExpression,
#                 ineq_array::Array{SimpleInequalityExpression, 1},
#                 var_array::Array{Int64, 1})
#
# Recursively reduces powers when the third power (u) in
# curr_inequality is equal to 1.
function ReduceThirdZero(curr_inequality::InequalityExpression,
                         ineq_array::Array{SimpleInequalityExpression, 1},
                         var_array::Array{Int, 1})
    p1 = curr_inequality.power1;
    p2 = curr_inequality.power2;
    p3 = curr_inequality.power3;
    n = (p1 + p2);
    @assert p3 == 0;
    @assert n >= 2;
    if (mod(n, 2) == 0)
        # n is even, so check even-ness of power1, power2
        if (p1 == p2)
            # We are seriously done
            new_ineq = SimpleInequalityExpression(curr_inequality.left_var_ind,
                                                  curr_inequality.t_var_ind,
                                                  curr_inequality.s_var_ind);
            return (vcat(ineq_array, new_ineq), var_array);
        end
        if (mod(p1, 2) == 0)
            # Both are even, reduce and repeat
            curr_inequality.power1 /= 2;
            curr_inequality.power2 /= 2;
            return ReducePowers(curr_inequality, ineq_array, var_array);
        else
            # p1 is odd, p2 is odd, and they are unequal. Go to two
            # equations, which are
            #
            # x^(n / 2) <= t^((p1 - 1)/2) s^((p2 - 1)/2) u,
            # u^2 <= ts
            #
            # If n/2 == 2, then do not recurse, but simply catenate the
            # inequality on and return
            new_ind = var_array[end] + 1;
            new_ineq = SimpleInequalityExpression(new_ind,
                                                  curr_inequality.t_var_ind,
                                                  curr_inequality.s_var_ind);
            curr_inequality.power1 = (p1 - 1) / 2;
            curr_inequality.power2 = (p2 - 1) / 2;
            curr_inequality.power3 = 1;
            curr_inequality.u_var_ind = new_ind;
            # Recurse on the reduced inequality
            return ReducePowers(curr_inequality,
                                vcat(ineq_array, new_ineq),
                                vcat(var_array, new_ind));
        end  # if (mod(p1, 2) == 0)
    else    # n is odd, and we must recurse differently.
        if (mod(p1, 2) == 1)
            new_ind = var_array[end] + 1;
            new_ineq = SimpleInequalityExpression(new_ind,
                                                  curr_inequality.t_var_ind,
                                                  curr_inequality.left_var_ind);
            curr_inequality.power1 = (p1 - 1) / 2;
            curr_inequality.power2 = p2 / 2;
            curr_inequality.power3 = 1;
            curr_inequality.u_var_ind = new_ind;
            return ReducePowers(curr_inequality,
                                vcat(ineq_array, new_ineq),
                                vcat(var_array, new_ind));
        else    # mod(p2, 2) == 1
            new_ind = var_array[end] + 1;
            new_ineq = SimpleInequalityExpression(new_ind,
                                                  curr_inequality.s_var_ind,
                                                  curr_inequality.left_var_ind);
            curr_inequality.power1 = p1 / 2;
            curr_inequality.power2 = (p2 - 1) / 2;
            curr_inequality.power3 = 1;
            curr_inequality.u_var_ind = new_ind;
            return ReducePowers(curr_inequality,
                                vcat(ineq_array, new_ineq),
                                vcat(var_array, new_ind));
        end
    end  # if mod(n, 2) == 0
end

# ReduceThirdOne(curr_inequality::InequalityExpression,
#                ineq_array::Array{SimpleInequalityExpression, 1},
#                var_array::Array{Int64, 1})
#
# Recursively reduces powers when the third power (u) in
# curr_inequality is equal to 1.
function ReduceThirdOne(curr_inequality::InequalityExpression,
                        ineq_array::Array{SimpleInequalityExpression, 1},
                        var_array::Array{Int, 1})
    p1 = curr_inequality.power1;
    p2 = curr_inequality.power2;
    p3 = curr_inequality.power3;
    n = (p1 + p2 + p3);
    @assert p3 == 1; # , "Must have third power 1");
    @assert n >= 3; #  "Must have power n >= 3");
    if (mod(n, 2) == 0)
        # Exactly one of p1, p2 is odd, find it and reduce
        if (mod(p1, 2) == 1)
            new_ind = var_array[end] + 1;
            new_ineq = SimpleInequalityExpression(new_ind,
                                                  curr_inequality.t_var_ind,
                                                  curr_inequality.u_var_ind);
            curr_inequality.power1 = (p1 - 1)/2;
            curr_inequality.power2 = (p2 / 2);
            curr_inequality.power3 = 1;
            curr_inequality.u_var_ind = new_ind;
            return ReducePowers(curr_inequality,
                                vcat(ineq_array, new_ineq),
                                vcat(var_array, new_ind));
        else    # mod(p2, 2) == 1
            new_ind = var_array[end] + 1;
            new_ineq = SimpleInequalityExpression(new_ind,
                                                  curr_inequality.s_var_ind,
                                                  curr_inequality.u_var_ind);
            curr_inequality.power1 = p1 / 2;
            curr_inequality.power2 = (p2 - 1) / 2;
            curr_inequality.power3 = 1;
            curr_inequality.u_var_ind = new_ind;
            return ReducePowers(curr_inequality,
                                vcat(ineq_array, new_ineq),
                                vcat(var_array, new_ind));
        end
    else    # mod(n, 2) == 1
        # Either both of p1, p2 are odd or both are even
        if (mod(p1, 2) == 0)
            new_ind = var_array[end] + 1;
            new_ineq = SimpleInequalityExpression(new_ind,
                                                  curr_inequality.left_var_ind,
                                                  curr_inequality.u_var_ind);
            curr_inequality.power1 = p1 / 2;
            curr_inequality.power2 = p2 / 2;
            curr_inequality.power3 = 1;
            curr_inequality.u_var_ind = new_ind;
            return ReducePowers(curr_inequality,
                                vcat(ineq_array, new_ineq),
                                vcat(var_array, new_ind));
        else    # mod(p1, 2) == 1, so mod(p2, 2) == 1
            if (p1 == 1 && p2 == 1)
                # There is a simpler reduction to return: we replace
                #
                # x^3 <= tsu
                #
                # with the three inequalities
                #
                # x^4 <= w^2 v^2,  w^2 <= ts,  v^2 <= ux,
                #
                # which is equivalent to x^2 <= wv, w^2 <= ts, v^2 <= ux
                new_w_ind = var_array[end] + 1;
                new_v_ind = var_array[end] + 2;
                new_x_ineq = SimpleInequalityExpression(curr_inequality.left_var_ind,
                                                        new_w_ind, new_v_ind);
                new_w_ineq = SimpleInequalityExpression(new_w_ind,
                                                        curr_inequality.t_var_ind,
                                                        curr_inequality.s_var_ind);
                new_v_ineq = SimpleInequalityExpression(new_v_ind,
                                                        curr_inequality.u_var_ind,
                                                        curr_inequality.left_var_ind);
                return (vcat(ineq_array, [new_x_ineq, new_w_ineq, new_v_ineq]),
                        vcat(var_array, [new_w_ind, new_v_ind]));
            else
                # Do the full reduction
                new_w_ind = var_array[end] + 1;
                new_v_ind = var_array[end] + 2;
                new_y_ind = var_array[end] + 3;
                new_y_ineq = SimpleInequalityExpression(new_y_ind,
                                                        new_w_ind, new_v_ind);
                new_w_ineq = SimpleInequalityExpression(new_w_ind,
                                                        curr_inequality.t_var_ind,
                                                        curr_inequality.s_var_ind);
                new_v_ineq = SimpleInequalityExpression(new_v_ind,
                                                        curr_inequality.u_var_ind,
                                                        curr_inequality.left_var_ind);
                curr_inequality.power1 = (p1 - 1) / 2;
                curr_inequality.power2 = (p2 - 1) / 2;
                curr_inequality.power3 = 2;
                curr_inequality.u_var_ind = new_y_ind;
                return ReducePowers(curr_inequality,
                                    vcat(ineq_array, [new_y_ineq, new_w_ineq, new_v_ineq]),
                                    vcat(var_array, [new_w_ind, new_v_ind, new_y_ind]));
            end
        end
    end
end

# ReduceThirdTwo(curr_inequality::InequalityExpression,
#                ineq_array::Array{SimpleInequalityExpression, 1},
#                var_array::Array{Int64, 1})
#
# Recursively reduces powers when the third power (u) in
# curr_inequality is equal to 2.
function ReduceThirdTwo(curr_inequality::InequalityExpression,
                        ineq_array::Array{SimpleInequalityExpression, 1},
                        var_array::Array{Int, 1})
    p1 = curr_inequality.power1;
    p2 = curr_inequality.power2;
    p3 = curr_inequality.power3;
    n = (p1 + p2 + p3);
    @assert p3 == 2;  # , "Must have third power 2");
    @assert n >= 6;  # , "Must have power n >= 6");
    @assert p1 >= 2 && p2 >= 2;  # , "Did not rearrange powers properly");
    if (mod(n, 2) == 0)
        # Either both p1, p2 are even or both are odd
        if (mod(p1, 2) == 0)
            # Simply reduce powers on everything by a factor of two
            curr_inequality.power1 /= 2;
            curr_inequality.power2 /= 2;
            curr_inequality.power3 /= 2;
            return ReducePowers(curr_inequality, ineq_array, var_array);
        else
            # Both are odd, so introduce two new variables
            new_y_ind = var_array[end] + 1;
            new_w_ind = var_array[end] + 2;
            new_y_ineq = SimpleInequalityExpression(new_y_ind,
                                                    curr_inequality.u_var_ind,
                                                    new_w_ind);
            new_w_ineq = SimpleInequalityExpression(new_w_ind,
                                                    curr_inequality.s_var_ind,
                                                    curr_inequality.t_var_ind);
            curr_inequality.power1 = (p1 - 1) / 2;
            curr_inequality.power2 = (p2 - 1) / 2;
            curr_inequality.power3 = 2;
            curr_inequality.u_var_ind = new_y_ind;
            return ReducePowers(curr_inequality,
                                vcat(ineq_array, [new_y_ineq, new_w_ineq]),
                                vcat(var_array, [new_y_ind, new_w_ind]));
        end
    else    # mod(n, 2) == 1, so exactly one of p1 and p2 is odd
        if (mod(p1, 2) == 1)
            new_y_ind = var_array[end] + 1;
            new_w_ind = var_array[end] + 2;
            new_y_ineq = SimpleInequalityExpression(new_y_ind,
                                                    curr_inequality.u_var_ind,
                                                    new_w_ind);
            new_w_ineq = SimpleInequalityExpression(new_w_ind,
                                                    curr_inequality.left_var_ind,
                                                    curr_inequality.t_var_ind);
            curr_inequality.power1 = (p1 - 1) / 2;
            curr_inequality.power2 = p2 / 2;
            curr_inequality.power3 = 2;
            curr_inequality.u_var_ind = new_y_ind;
            return ReducePowers(curr_inequality,
                                vcat(ineq_array, [new_y_ineq, new_w_ineq]),
                                vcat(var_array, [new_y_ind, new_w_ind]));
        else    # mod(p2, 2) == 1
            new_y_ind = var_array[end] + 1;
            new_w_ind = var_array[end] + 2;
            new_y_ineq = SimpleInequalityExpression(new_y_ind,
                                                    curr_inequality.u_var_ind,
                                                    new_w_ind);
            new_w_ineq = SimpleInequalityExpression(new_w_ind,
                                                    curr_inequality.left_var_ind,
                                                    curr_inequality.s_var_ind);
            curr_inequality.power1 = p1 / 2;
            curr_inequality.power2 = (p2 - 1) / 2;
            curr_inequality.power3 = 2;
            curr_inequality.u_var_ind = new_y_ind;
            return ReducePowers(curr_inequality,
                                vcat(ineq_array, [new_y_ineq, new_w_ineq]),
                                vcat(var_array, [new_y_ind, new_w_ind]));
        end
    end
end

end  # module PSOCP
