# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

function get_counts(model)
    n_constraints = 0
    total_len_constraints = 0
    n_variables = MOI.get(model, MOI.NumberOfVariables())
    n_constants = 0
    for (F, S) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
        n_constraints += MOI.get(model, MOI.NumberOfConstraints{F,S}())

        for idx in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
            set = MOI.get(model, MOI.ConstraintSet(), idx)
            total_len_constraints += MOI.dimension(set)
            f = MOI.get(model, MOI.ConstraintFunction(), idx)
            n_constants += count_size(f)
        end
    end
    return (; n_variables, n_constraints, total_len_constraints, n_constants)
end

# If we use other MOI functions, we will need to add methods.
function count_size(f::MOI.VectorAffineFunction)
    return length(f.terms) + length(f.constants)
end

function show_moi_counts(io::IO, model)
    c = get_counts(model)
    maybe_s(n) = n == 1 ? "" : "s"
    print(io, c.n_variables, " scalar variable", maybe_s(c.n_variables))
    print(
        io,
        ", ",
        underscorise(c.n_constraints),
        " constraint",
        maybe_s(c.n_constraints),
    )
    if c.n_constraints > 0
        print(
            io,
            " (",
            underscorise(c.total_len_constraints),
            " element",
            maybe_s(c.total_len_constraints),
            ")",
        )
    end
    print(
        io,
        ", ",
        underscorise(c.n_constants),
        " sparse constant",
        maybe_s(c.n_constants),
    )
    return nothing
end
