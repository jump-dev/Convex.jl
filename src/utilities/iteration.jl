import Base.iterate
export iterate

function iterate(x::Variable, (el, s)=(x[1], 0))
    return s >= length(x) ? nothing : (el, (x[s+1], s+1))
end
