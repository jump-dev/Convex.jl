export 
    ConicProblem, 
    IneqConicProblem,
    ECOSConicProblem

# The ConicProblem type stores the problem data
# a ConicProblem instance corresponds to the problem

    # minimize     c'*x
    # subject to   Ax == b 
    #              x \in cones

# The parameter c is the objective vector, the parameter A is the constraint matrix (typically sparse), the parameter b is the vector of right-hand side values, and cones is a list of (Symbol,vars) tuples, where Symbol is one of the above recognized cones and vars is a list of indices of variables which belong to this cone (may be given as a Range). All variables must be listed in exactly one cone, and the indices given must correspond to the order of the columns in in the constraint matrix A. Cones may be listed in any order, and cones of the same class may appear multiple times. For the semidefinite cone, the number of variables present must be a square integer n corresponding to a sqrt(n) x sqrt(n) matrix; 
# variables should be listed in column-major, or by symmetry, row-major order.
type ConicProblem
    c::Array{Float64, 1}
    A::Array{Float64, 1}
    b::Array{Float64, 1}
    cones::Array # Array of (Symbol,Range) or (Symbol,Array{Integer}) tuples
end

# The IneqConicProblem type stores the problem data in conic inequality form
# an IneqConicProblem instance corresponds to the problem

    # minimize     c'*x
    # subject to   Ax == b 
    #              Gx \leq_cones h

# (this corresponds to the form accepted by the solver ECOS)
type IneqConicProblem
    c::Array{Float64, 1}
    A::Array{Float64, 1}
    b::Array{Float64, 1}
    G::Array{Float64, 1}
    h::Array{Float64, 1}
    cones::Array # Array of (Symbol,Range) or (Symbol,Array{Integer}) tuples
end

type ECOSConicProblem
    n
    m
    p
    l
    ncones
    q
    G
    c
    h
    A
    b
end
ECOSConicProblem(;n=n, m=m, p=p, l=l, ncones=ncones, q=q, G=G, c=c, h=h, A=A, b=b) = 
    ECOSConicProblem(n, m, p, l, ncones, q, G, c, h, A, b)

# CAUTION: For now, we assume we are solving a linear program.
# Loops over the objective and constraints to get the canonical constraints array.
# It then calls create_ecos_matrices which will create the inequality/equality matrices
# and their corresponding coefficients, which are then passed to ecos_solve
function ECOSConicProblem(problem::Problem)
    canonical_constraints_array = canonical_constraints(problem)
    m, n, p, l, ncones, q, G, h, A, b, variable_index, eq_constr_index, ineq_constr_index =
        create_ecos_matrices(canonical_constraints_array)

    # Now, all we need to is create c
    c = zeros(n, 1)
    objective = problem.objective
    if objective.vexity != :constant
        uid = objective.uid
        c[variable_index[uid] : variable_index[uid] + objective.size[1] - 1] = 1
    end  

    return ECOSConicProblem(n=n, m=m, p=p, l=l, ncones=ncones, q=q, G=G, c=c, h=h, A=A, b=b), variable_index, eq_constr_index, ineq_constr_index
end

function IneqConicProblem(p::ConicProblem)
    return p
end

function IneqConicProblem(p::ECOSConicProblem)
    cones = [(:NonNeg,p.q)]
    for dim in p.l
        lastidx = cones[-1][-1]
        append!(cones,(:SOC,lastidx+1:lastidx+dim))
    end
    return p
end

function ConicProblem(ip::IneqConicProblem)
    return ip
end

function ECOSConicProblem(ip::IneqConicProblem)
    q = 0; l = [];
    for (cone,idx) in ip.cones
        if cone == :free
            continue
        elseif cone == :NonNeg
            q += length(idx)
        elseif cone == :SOC
            l += 1
        else
            error("ECOS does not support cone $cone")
        end
    end
    return ECOSConicProblem(n=n, m=m, p=p, l=l, ncones=ncones, q=q, G=G, c=c, h=h, A=A, b=b)
end