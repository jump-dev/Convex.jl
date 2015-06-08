import Base: start, next, done
export start, next, done

start(x::Variable) = 1
next(x::Variable, s::Int) = x[s], s+1
done(x::Variable, s::Int) = (s > length(x))