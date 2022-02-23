import Base.iterate

function iterate(x::AbstractVariable, s = 0)
    return s >= length(x) ? nothing : (x[s+1], s + 1)
end
