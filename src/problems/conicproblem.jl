import Base.show

export 
    ConicProblem, 
    IneqConicProblem,
    ConicProblemSolution, 
    IneqConicProblemSolution,
    ECOSConicProblem,
    cones,
    problemtype,
    show

cones_in = {:LP=>Set({:NonNeg}),
            :SOC=>Set({:NonNeg, :SOC}),
            :SDP=>Set({:NonNeg, :SOC, :SDP})}

import Base.convert
convert(::Type{Range}, i::Integer) = i:i 

# The ConicProblem type stores the problem data
# a ConicProblem instance corresponds to the problem

    # minimize     c'*x
    # subject to   Ax == b 
    #              x \in cones

# The parameter c is the objective vector, the parameter A is the constraint matrix (typically sparse), the parameter b is the vector of right-hand side values, and cones is a list of (Symbol,vars) tuples, where Symbol is one of the above recognized cones and vars is a list of indices of variables which belong to this cone (may be given as a Range). All variables must be listed in exactly one cone, and the indices given must correspond to the order of the columns in in the constraint matrix A. Cones may be listed in any order, and cones of the same class may appear multiple times. For the semidefinite cone, the number of variables present must be a square integer n corresponding to a sqrt(n) x sqrt(n) matrix; 
# variables should be listed in column-major, or by symmetry, row-major order.
type ConicProblem{T}
    c::Array{T, 2}
    A::AbstractArray{T, 2}
    b::Array{T, 2}
    cones::Array # Array of (Symbol,Iterable{Integer}) tuples

    function ConicProblem(c::Array{T, 2}, A::AbstractArray{T, 2}, b::Array{T, 2}, cones::Array)
        # verify all variables are in exactly one cone
        free_vars = Set(1:length(c))
        set_cones = Array((Symbol,Set), length(cones)+1)
        for i= 1:length(cones)
            cone_type, idxs = cones[i]
            set_cones[i] = cone_type, Set(idxs)
            free_vars = setdiff(free_vars, set_cones[i][2])
        end
        if length(free_vars) > 0
            set_cones[end] = (:Free, free_vars)
            cones = set_cones
        end
        return new(c, A, b, cones)
    end
end
ConicProblem{T<:Real}(c::Array{T, 2}, A::AbstractArray{T, 2}, b::Array{T, 2}, cones::Array) = ConicProblem{T}(c, A, b, cones)

function show(io::IO, p::ConicProblem)
  print(io, "ConicProblem: minimize c^T x subject to A x <= b, x in K\n
            c:\n$(p.c)
            A:\n$(p.A)
            b:\n$(p.b)
            K:\n$(p.cones)")
end
cones(cp::ConicProblem) = Set(keys(cp.cones))
function problemtype(cp::ConicProblem)
    mycones = cones(cp)
    for problem_type in [:LP, :SOC, :SDP]
        if isempty(mycones - cones_in(problem_type))
            return problem_type
        end
    end
    error("No solvers found to solve problems with cones $mycones")
end

# The IneqConicProblem type stores the problem data in conic inequality form
# an IneqConicProblem instance corresponds to the problem

    # minimize     c'*x
    # subject to   Ax == b 
    #              Gx \leq_cones h

# (this corresponds to the form accepted by the solver ECOS)
type IneqConicProblem{T<:Number}
    c::Array{T, 2}
    A::AbstractArray{T, 2}
    b::Array{T, 2}
    G::AbstractArray{T, 2}
    h::Array{T, 2}
    cones::Array # eg Array{(Symbol,Int64),1} # Array of (Symbol,Range) or (Symbol,Array{Integer}) tuples
end
function show(io::IO, p::IneqConicProblem)
  print(io, "IneqConicProblem: minimize c^T x subject to A x <= b, G x <=_K h\n
            c:\n$(p.c)
            A:\n$(p.A)
            b:\n$(p.b)
            G:\n$(p.G)
            h:\n$(p.h)
            K:\n$(p.cones)")
end

cones(p::IneqConicProblem) = Set(keys(p.cones))

type ECOSConicProblem
    n          # number of variables
    m          # number of inequality constraints, ie, size(G,1)
    p          # number of equality constraints
    l          # dimension of the positive orthant constraints
    ncones     # number of SOC cones
    q          # array of integers of length ncones, where each element defines the dimension of the cone
    G      # inequality constraint Gx \leq_K h, m,n = size(G)
    c      # objective function
    h      # rhs of inequality constraints
    A      # equality constraints, p,n = size(A)
    b      # rhs of equality constraints
end
ECOSConicProblem(;n=n, m=m, p=p, l=l, ncones=ncones, q=q, G=G, c=c, h=h, A=A, b=b) = 
    ECOSConicProblem(n, m, p, l, ncones, q, G, c, h, A, b)

function cones(p::ECOSConicProblem)
    mycones = Set()
    if sum(p.l) > 0
        push!(mycones, :SOC)
    end 
    if p.q > 0
        push!(mycones, :NonNeg)
    end
    return mycones
end

function IneqConicProblem(p::ConicProblem)
    m,n = size(p.A)
    G = -speye(n)
    h = zeros(n,1)
    return IneqConicProblem(p.c,p.A,p.b,G,h,p.cones)
end

function IneqConicProblem(p::ECOSConicProblem)
    if p.l > 0
        cones = [(:NonNeg,1:p.l)]
    else 
        cones = Any[]
    end
    lastidx = p.l
    for dim in p.q
        push!(cones,(:SOC,lastidx+1:lastidx+dim))
        lastidx += dim
    end
    n = length(p.c); T = eltype(p.c)
    G = (p.G == nothing ? spzeros(T,0,n) : p.G )
    A = (p.A == nothing ? spzeros(T,0,n) : p.A )
    b = (p.b == nothing ? Array(T,(0,1)) : p.b )
    h = (p.h == nothing ? Array(T,(0,1)) : p.h )
    return IneqConicProblem(p.c,A,b,G,h,cones)
end

function ConicProblem(ip::IneqConicProblem)
    nslacks = length(ip.h)
    neq, nvars = size(cp.A)
    c = [ip.c; zeros(nslacks)]
    nslacks, n = size(cp.G)
    A = blkdiag(ip.A, speye(nslacks)) 
    A[neq+1:end, 1:n] = cp.G
    b = [ip.b; ip.h]
    cones = [(cone, idx_range + nvars) for (cone, idx_range) in ip.cones]
    push!(cones, (:free, 1:nvars)) # XXX check this
    return ConicProblem(c,A,b,cones)
end

function ECOSConicProblem(ip::IneqConicProblem)
    # TODO: preserve sparsity
    l = 0; q = Int64[];
    nonneg_indices = Range[]
    soc_indices = Range[]
    ncones = 0
    for (cone,idx) in ip.cones
        if cone == :free
            continue
        elseif cone == :NonNeg
            l += length(idx)
            push!(nonneg_indices, idx)
        elseif cone == :SOC
            push!(q, length(idx))
            push!(soc_indices, idx)
            ncones += 1
        else
            error("ECOS does not support cone $cone")
        end
    end
    # rearrange rows so nonneg cone comes before SOC
    nonneg_indices,soc_indices = map(l->vcat(l...), (nonneg_indices,soc_indices))
    G = vcat(ip.G[nonneg_indices,:], ip.G[soc_indices,:])
    h = vcat(ip.h[nonneg_indices], ip.h[soc_indices])
    m,n = size(G)
    p,n = size(ip.A)
    return ECOSConicProblem(n=n, m=m, p=p, l=l, ncones=ncones, q=q, G=G, c=ip.c, h=h, A=ip.A, b=ip.b)
end

ECOSConicProblem(p::ConicProblem) = ECOSConicProblem(IneqConicProblem(p))
ConicProblem(p::ECOSConicProblem) = ConicProblem(IneqConicProblem(p))