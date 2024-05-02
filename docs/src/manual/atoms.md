# Atoms

Convex.jl supports the following functions. These functions may be composed
according to the [DCP](http://dcp.stanford.edu) composition rules to form new
convex, concave, or affine expressions.
## `*`
```@docs
Base.:*(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `+`
```@docs
Base.:+(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `-`
```@docs
Base.:-(::Convex.AbstractExpr)
Base.:-(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `/`
```@docs
Base.:/(::Convex.AbstractExpr, ::Convex.Value)
```

## `.*`
```@docs
Base.Broadcast.broadcasted(::typeof(*), ::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `./`
```@docs
Base.Broadcast.broadcasted(::typeof(/), ::Convex.AbstractExpr, ::Convex.Value)
```

## `.^`
```@docs
Base.Broadcast.broadcasted(::typeof(^), ::Convex.AbstractExpr, ::Int)
```

## `abs`
```@docs
Base.abs(::Convex.AbstractExpr)
```

## `abs2`
```@docs
Base.abs2(::Convex.AbstractExpr)
```

## `adjoint`
```@docs
Convex.adjoint(::Convex.AbstractExpr)
```

## `conj`
```@docs
Base.conj(::Convex.AbstractExpr)
```

## `conv`
```@docs
Convex.conv(::Convex.Value, ::Convex.AbstractExpr)
```

## `diag`
```@docs
Convex.diag(::Convex.AbstractExpr, ::Int = 0)
```

## `diagm`
```@docs
Convex.diagm(::Convex.AbstractExpr)
```

## `dot`
```@docs
Convex.dot(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `dotsort`
```@docs
Convex.dotsort(::Convex.AbstractExpr, ::Convex.Value)
```

## `eigmax`
```@docs
Convex.eigmax(::Convex.AbstractExpr)
```

## `eigmin`
```@docs
Convex.eigmin(::Convex.AbstractExpr)
```

## `entropy`
```@docs
Convex.entropy(::Convex.AbstractExpr)
```

## `entropy_elementwise`
```@docs
Convex.entropy_elementwise(::Convex.AbstractExpr)
```

## `exp`
```@docs
Base.exp(::Convex.AbstractExpr)
```

## `geomean`
```@docs
Convex.geomean(::Union{Convex.AbstractExpr,Convex.Value}...)
```

## `hcat`
```@docs
Base.hcat(::Convex.AbstractExpr...)
```

## `hinge_loss`
```@docs
Convex.hinge_loss(::Convex.AbstractExpr)
```

## `huber`
```@docs
Convex.huber(::Convex.AbstractExpr, ::Real = 1.0)
```

## `hvcat`
```@docs
Base.hvcat(
    ::Tuple{Vararg{Int}},
    ::Union{Convex.AbstractExpr,Convex.Value}...,
)
```

## `imag`
```@docs
Base.imag(::Convex.AbstractExpr)
```

## `inner_product`
```@docs
Convex.inner_product(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `invpos`
```@docs
Convex.invpos(::Convex.AbstractExpr)
```

## `kron`
```@docs
Base.kron(::Convex.Value, ::Convex.AbstractExpr)
```

## `lieb_ando`
```@docs
Convex.lieb_ando(
    ::Convex.AbstractExpr,
    ::Convex.AbstractExpr,
    ::Union{AbstractMatrix,Convex.Constant},
    ::Rational,
)
```

## `log`
```@docs
Base.log(::Convex.AbstractExpr)
```

## `log_perspective`
```@docs
Convex.log_perspective(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `logdet`
```@docs
Convex.logdet(::Convex.AbstractExpr)
```

## `logisticloss`
```@docs
Convex.logisticloss(::Convex.AbstractExpr)
```

## `logsumexp`
```@docs
Convex.logsumexp(::Convex.AbstractExpr)
```

## `matrixfrac`
```@docs
Convex.matrixfrac(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `max`
```@docs
Base.max(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `maximum`
```@docs
Base.maximum(::Convex.AbstractExpr)
```

## `min`
```@docs
Base.min(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `minimum`
```@docs
Base.minimum(::Convex.AbstractExpr)
```

## `neg`
```@docs
Convex.neg(::Convex.AbstractExpr)
```

## `norm`
```@docs
Convex.norm(::Convex.AbstractExpr, ::Real = 2)
```

## `norm2`
```@docs
Convex.norm2(::Convex.AbstractExpr)
```

## `nuclearnorm`
```@docs
Convex.nuclearnorm(::Convex.AbstractExpr)
```

## `opnorm`
```@docs
Convex.opnorm(::Convex.AbstractExpr, ::Real = 2)
```

## `partialtrace`
```@docs
Convex.partialtrace(x, ::Int, ::Vector)
```

## `partialtranspose`
```@docs
Convex.partialtranspose(
    ::Union{AbstractMatrix,Convex.AbstractExpr},
    ::Int,
    ::Vector,
)
```

## `pos`
```@docs
Convex.pos(::Convex.AbstractExpr)
```

## `qol_elementwise`
```@docs
Convex.qol_elementwise(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `quadform`
```@docs
Convex.quadform(::Convex.Value, ::Convex.AbstractExpr; kwargs...)
```

## `quadoverlin`
```@docs
Convex.quadoverlin(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `quantum_entropy`
```@docs
Convex.quantum_entropy(::Convex.AbstractExpr, ::Integer = 3, ::Integer = 3)
```

## `quantum_relative_entropy`
```@docs
Convex.quantum_relative_entropy(
    ::Convex.AbstractExpr,
    ::Convex.AbstractExpr,
    ::Integer = 3,
    ::Integer = 3,
)
```

## `rationalnorm`
```@docs
Convex.rationalnorm(::Convex.AbstractExpr, ::Rational{Int})
```

## `real`
```@docs
Base.real(::Convex.AbstractExpr)
```

## `relative_entropy`
```@docs
Convex.relative_entropy(::Convex.AbstractExpr, ::Convex.AbstractExpr)
```

## `reshape`
```@docs
Base.reshape(::Convex.AbstractExpr, ::Int, ::Int)
```

## `rootdet`
```@docs
Convex.rootdet(::Convex.AbstractExpr)
```

## `sigmamax`
```@docs
Convex.sigmamax(::Convex.AbstractExpr)
```

## `sqrt`
```@docs
Base.sqrt(::Convex.AbstractExpr)
```

## `square`
```@docs
Convex.square(::Convex.AbstractExpr)
```

## `sum`
```@docs
Base.sum(::Convex.AbstractExpr; dims = :)
```

## `sumlargest`
```@docs
Convex.sumlargest(::Convex.AbstractExpr, ::Int)
```

## `sumlargesteigs`
```@docs
Convex.sumlargesteigs(::Convex.AbstractExpr, ::Int)
```

## `sumsmallest`
```@docs
Convex.sumsmallest(::Convex.AbstractExpr, ::Int)
```

## `sumsquares`
```@docs
Convex.sumsquares(::Convex.AbstractExpr)
```

## `tr`
```@docs
Convex.tr(::Convex.AbstractExpr)
```

## `trace_logm`
```@docs
Convex.trace_logm(
    ::Convex.AbstractExpr,
    ::Union{AbstractMatrix,Convex.Constant},
    ::Integer = 3,
    ::Integer = 3,
)
```

## `trace_mpower`
```@docs
Convex.trace_mpower(
    ::Convex.AbstractExpr,
    ::Rational,
    ::Union{AbstractMatrix,Convex.Constant},
)
```

## `transpose`
```@docs
Convex.transpose(::Convex.AbstractExpr)
```

## `vcat`
```@docs
Base.vcat(::Convex.AbstractExpr...)
```

## `vec`
```@docs
Base.vec(::Convex.AbstractExpr)
```

