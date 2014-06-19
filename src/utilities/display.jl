import Base.show
export show 

function show(io::IO, x::Constant)
    print(io,"$(x.value)")
end

function show(io::IO, x::Variable)
    print(io,"Variable$(x.size)")
end

function show(io::IO,c::CvxConstr)
    print(io,"CvxConstr: $(c.lhs) $(c.head) $(c.rhs)")
end

function show(io::IO,e::CvxExpr)
    print(io,"CvxExpr($(e.size),$(e.vexity))")
end

function show(io::IO,p::Problem)
    print(io,"Problem: $(p.head) $(p.objective) \n \t subject to\n\t\t")
    print_joined(io,p.constraints,"\n\t\t")
end