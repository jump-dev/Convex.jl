# Interface from CVX.jl canonical form to MathProgBase canonical form
using ECOS
using MathProgBase

function solve_mpb_ecos(canonical_cons)
    # First task is to build a mapping from CVX's multidimensional
    # variable objects to a flat list of variables
    num_var = 0
    var_indices = Dict{Int,Int}()
    for con in canonical_cons
        #println("Con")
        for (j, var) in enumerate(con.vars)
            #dump(j)
            #dump(var)
            # If haven't seen this variable before...
            if !haskey(var_indices, var)
                # ... give it next available index
                var_indices[var] = num_var + 1
                # Track how many variables we just added
                num_var += size(con.coeffs[j], 2)
            end
        end
    end
    free_cones = num_var

    # We now know how many variables we've got, allocate storage
    # for the coefficient matrix and the rhs
    A = zeros(0, num_var)
    b = zeros(0)
    for con in canonical_cons
        #dump(con)
        if con.is_eq
            # Easy case - constraint is a simple equality
            num_row, num_col = size(con.coeffs[1])
            new_rows = zeros(num_row, num_var)
            for (var_pos, var) in enumerate(con.vars)
                base_idx = var_indices[var]
                for i = 1:num_row
                    for j = 1:num_col
                        new_rows[i,base_idx+j-1] = con.coeffs[var_pos][i,j]
                    end
                end
            end
            for i = 1:num_row
                push!(b, con.constant[i])
            end
            A = vcat(A, new_rows)
        else
            # Inequality - need to introduce slacks as well
            num_row, num_col = size(con.coeffs[1])
            num_var += num_row
            A = hcat(A, zeros(size(A,1),num_row))
            new_rows = zeros(num_row, num_var)
            for (var_pos, var) in enumerate(con.vars)
                base_idx = var_indices[var]
                for i = 1:num_row
                    for j = 1:num_col
                        new_rows[i,base_idx+j-1] = con.coeffs[var_pos][i,j]
                    end
                end
            end
            for i = 1:num_row
                push!(b, con.constant[i])
                new_rows[i,num_var-num_row+i] = 1.0
            end
            A = vcat(A, new_rows)
        end
    end

    c = zeros(num_var)
    c[free_cones] = 1.0
    cones = [(:Free,1:free_cones),(:NonNeg,(free_cones+1):num_var)]

    s = ECOSSolver()
    m = MathProgBase.model(s)
    MathProgBase.loadconicproblem!(m, c, A, b, cones)
    MathProgBase.optimize!(m)
    
    return Solution(MathProgBase.getsolution(m), 
                    MathProgBase.getsolution(m),  # Dual equality
                    MathProgBase.getsolution(m),  # Dual inequality
                    0)  # Really should just use status(m)
end