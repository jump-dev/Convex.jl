module Convex

### modeling framework
include("dcp.jl")
include("expressions.jl")
include("variable.jl")
include("constant.jl")
include("conic_form.jl")
include("constraints.jl")
include("problems.jl")
include("solution.jl")

### affine atoms
include("add_subtract.jl")
include("multiply_divide.jl")
include("sum.jl")
include("abs.jl")
include("transpose.jl")
include("index.jl")
include("diag.jl")
include("stack.jl")
include("dot.jl")

### elementwise atoms
include("maximum.jl")
include("minimum.jl")

### SOC atoms
include("norm2.jl")

### SDP atoms
include("matrix_norm.jl")

### exponential atoms

### utilities
include("utilities/show.jl")
end
