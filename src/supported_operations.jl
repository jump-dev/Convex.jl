# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

"""
    Base.:+(x::Convex.AbstractExpr, y::Convex.AbstractExpr)
    Base.:+(x::Convex.Value, y::Convex.AbstractExpr)
    Base.:+(x::Convex.AbstractExpr, y::Convex.Value)

The addition operator \$x + y\$.

## Examples

Applies to scalar expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> x + 1
+ (affine; real)
├─ real variable (id: 110…477)
└─ [1;;]
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = [1, 2, 3];

julia> atom = x + y
+ (affine; real)
├─ 3-element real variable (id: 458…482)
└─ [1; 2; 3;;]

julia> size(atom)
(3, 1)
```
"""
Base.:+(x::AbstractExpr, y::AbstractExpr) = AdditionAtom(x, y)

Base.:+(x::Value, y::AbstractExpr) = constant(x) + y

Base.:+(x::AbstractExpr, y::Value) = x + constant(y)

"""
    Base.:-(x::Convex.AbstractExpr, y::Convex.AbstractExpr)
    Base.:-(x::Convex.Value, y::Convex.AbstractExpr)
    Base.:-(x::Convex.AbstractExpr, y::Convex.Value)

The subtraction operator \$x - y\$.

## Examples

Applies to scalar expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> x - 1
+ (affine; real)
├─ real variable (id: 161…677)
└─ [-1;;]
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = [1, 2, 3];

julia> atom = y - x
+ (affine; real)
├─ [1; 2; 3;;]
└─ Convex.NegateAtom (affine; real)
   └─ 3-element real variable (id: 242…661)
```
"""
Base.:-(x::AbstractExpr, y::AbstractExpr) = x + (-y)

Base.:-(x::Value, y::AbstractExpr) = constant(x) + (-y)

Base.:-(x::AbstractExpr, y::Value) = x + constant(-y)

"""
    Base.:-(x::Convex.AbstractExpr)

The univariate negation operator \$-x\$.

## Examples

Applies to scalar expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> -x
Convex.NegateAtom (affine; real)
├─ real variable (id: 161…677)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = -x
Convex.NegateAtom (affine; real)
└─ 3-element real variable (id: 137…541)

julia> size(atom)
(3, 1)
```
"""
Base.:-(x::AbstractExpr) = NegateAtom(x)

Base.:-(x::Union{Constant,ComplexConstant}) = constant(-evaluate(x))

"""
    Base.:*(x::Convex.AbstractExpr, y::Convex.AbstractExpr)

The binary multiplication operator \$x \\times y\$.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
ulia> x = Variable();

julia> 2 * x
* (affine; real)
├─ [2;;]
└─ real variable (id: 709…007)
```

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = [1, 2, 3];

julia> x' * y
* (affine; real)
├─ reshape (affine; real)
│  └─ * (affine; real)
│     ├─ 3×3 SparseArrays.SparseMatrixCSC{Int64, Int64} with 3 stored entries
│     └─ reshape (affine; real)
│        └─ …
└─ [1; 2; 3;;]
```
"""
function Base.:*(x::AbstractExpr, y::AbstractExpr)
    if isequal(x, y) && x.size == (1, 1)
        return square(x)
    end
    return MultiplyAtom(x, y)
end

Base.:*(x::Value, y::AbstractExpr) = MultiplyAtom(constant(x), y)

Base.:*(x::AbstractExpr, y::Value) = MultiplyAtom(x, constant(y))

"""
    Base.:/(x::Convex.AbstractExpr, y::Convex.Value)

The binary division operator \$\\frac{x}{y}\$.

## Examples

Applies to a scalar expression:
```jldoctest; filter=r"id: [0-9]+…[0-9]+"
ulia> x = Variable();

julia> x / 2
```
and element-wise to a matrix:
```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = x / 2
* (affine; real)
├─ 3-element real variable (id: 129…611)
└─ [0.5;;]

julia> size(atom)
(3, 1)
```
"""
Base.:/(x::AbstractExpr, y::Value) = MultiplyAtom(x, constant(1 ./ y))

function Base.:/(x::Value, y::AbstractExpr)
    if size(y) != (1, 1)
        error("cannot divide by a variable of size $(size(y))")
    end
    return MultiplyAtom(constant(x), invpos(y))
end

"""
    x::Convex.AbstractExpr .* y::Convex.AbstractExpr

Element-wise multiplication between matrices `x` and `y`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> atom = x .* 2
* (affine; real)
├─ 2-element real variable (id: 197…044)
└─ [2;;]

julia> atom = x .* [2, 4]
.* (affine; real)
├─ 2-element real variable (id: 197…044)
└─ [2; 4;;]

julia> size(atom)
(2, 1)
```
"""
function Base.Broadcast.broadcasted(
    ::typeof(*),
    x::AbstractExpr,
    y::AbstractExpr,
)
    if isequal(x, y)
        return square(x)
    elseif x.size == (1, 1) || y.size == (1, 1)
        return x * y
    end
    return BroadcastMultiplyAtom(x, y)
end

function Base.Broadcast.broadcasted(::typeof(*), x::Value, y::AbstractExpr)
    return constant(x) .* y
end

function Base.Broadcast.broadcasted(::typeof(*), x::AbstractExpr, y::Value)
    return x .* constant(y)
end

"""
    x::Convex.AbstractExpr ./ y::Convex.AbstractExpr

Element-wise division between matrices `x` and `y`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> atom = x ./ 2
* (affine; real)
├─ 2-element real variable (id: 875…859)
└─ [0.5;;]

julia> atom = x ./ [2, 4]
.* (affine; real)
├─ 2-element real variable (id: 875…859)
└─ [0.5; 0.25;;]

julia> size(atom)
(2, 1)
```
"""
function Base.Broadcast.broadcasted(::typeof(/), x::AbstractExpr, y::Value)
    return x .* constant(1 ./ y)
end

function Base.Broadcast.broadcasted(::typeof(/), x::Value, y::AbstractExpr)
    return constant(x) .* invpos(y)
end

"""
    x::Convex.AbstractExpr .^ k::Int

Element-wise exponentiation of `x` to the power of `k`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> atom = x .^ 2
qol_elem (convex; positive)
├─ 2-element real variable (id: 131…737)
└─ [1.0; 1.0;;]

julia> size(atom)
(2, 1)
```
"""
function Base.Broadcast.broadcasted(::typeof(^), x::AbstractExpr, k::Int)
    if k != 2
        error("raising variables to powers other than 2 is not implemented")
    end
    return square(x)
end

function Base.Broadcast.broadcasted(
    ::typeof(Base.literal_pow),
    ::typeof(^),
    x::AbstractExpr,
    ::Val{k},
) where {k}
    return Base.Broadcast.broadcasted(^, x, k)
end

"""
    Base.abs(x::Convex.AbstractExpr)

The epigraph of \$|x|\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> abs(x)
abs (convex; positive)
└─ real variable (id: 103…720)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = abs(x)
abs (convex; positive)
└─ 3-element real variable (id: 389…882)

julia> size(atom)
(3, 1)
```
"""
Base.abs(x::AbstractExpr) = AbsAtom(x)

"""
    Base.abs2(x::Convex.AbstractExpr)

The epigraph of \$|x|^2\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> abs2(x)
qol_elem (convex; positive)
├─ abs (convex; positive)
│  └─ real variable (id: 319…413)
└─ [1.0;;]
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = abs2(x)
qol_elem (convex; positive)
├─ abs (convex; positive)
│  └─ 3-element real variable (id: 123…996)
└─ [1.0; 1.0; 1.0;;]

julia> size(atom)
(3, 1)
```
"""
Base.abs2(x::AbstractExpr) = square(abs(x))

"""
    LinearAlgebra.adjoint(x::AbstractExpr)

The transpose of the conjugated matrix `x`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = ComplexVariable(2, 2);

julia> atom = adjoint(x)
reshape (affine; complex)
└─ * (affine; complex)
   ├─ 4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 4 stored entries
   └─ reshape (affine; complex)
      └─ conj (affine; complex)
         └─ …

julia> size(atom)
(2, 2)
```
"""
LinearAlgebra.adjoint(x::AbstractExpr) = transpose(conj(x))

function LinearAlgebra.adjoint(x::Union{Constant,ComplexConstant})
    return constant(adjoint(evaluate(x)))
end

"""
    Base.conj(x::Convex.AbstractExpr)

The complex conjugate of `x`.

If `x` is real, this function returns `x`.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = ComplexVariable();

julia> conj(x)
conj (affine; complex)
└─ complex variable (id: 180…137)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
conj (affine; complex)
└─ complex variable (id: 180…137)

julia> x = ComplexVariable(3);

julia> atom = conj(x)
conj (affine; complex)
└─ 3-element complex variable (id: 104…031)

julia> size(atom)
(3, 1)
```
"""
function Base.conj(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return ConjugateAtom(x)
    end
    return x
end

Base.conj(x::Constant) = x

Base.conj(x::ComplexConstant) = ComplexConstant(real(x), -imag(x))

"""
    LinearAlgebra.diag(x::Convex.AbstractExpr, k::Int = 0)

Return the `k`-th diagonnal of the matrix `X` as a column vector.

## Examples

Applies to a single square matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = diag(x, 0)
diag (affine; real)
└─ 2×2 real variable (id: 724…318)

julia> size(atom)
(2, 1)

julia> atom = diag(x, 1)
diag (affine; real)
└─ 2×2 real variable (id: 147…856)

julia> size(atom)
(1, 1)
```
"""
LinearAlgebra.diag(x::AbstractExpr, k::Int = 0) = DiagAtom(x, k)

"""
    LinearAlgebra.diagm(x::Convex.AbstractExpr)

Create a diagonal matrix out of the vector `x`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> atom = diagm(x)
diagm (affine; real)
└─ 2-element real variable (id: 541…968)

julia> size(atom)
(2, 2)
```
"""
LinearAlgebra.diagm(x::AbstractExpr) = DiagMatrixAtom(x)

function LinearAlgebra.diagm((d, x)::Pair{<:Integer,<:AbstractExpr})
    if d != 0
        msg = "[DiagMatrixAtom] only the main diagonal is supported. Got `d=$d`"
        throw(ArgumentError(msg))
    end
    return diagm(x)
end

LinearAlgebra.Diagonal(x::AbstractExpr) = diagm(x)

"""
    dotsort(x::Convex.AbstractExpr, y::Convex.Value)
    dotsort(x::Convex.Value, y::Convex.AbstractExpr)

Computes `dot(sort(x), sort(y))`, where `x` or `y` is constant.

For example, if `x = Variable(6)` and `y = [1 1 1 0 0 0]`, this atom computes
the sum of the three largest elements of `x`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(4);

julia> atom = dotsort(x, [1, 0, 0, 1])
dotsort (convex; real)
└─ 4-element real variable (id: 128…367)

julia> size(atom)
(1, 1)
```
"""
dotsort(a::AbstractExpr, b::Value) = DotSortAtom(a, b)
dotsort(b::Value, a::AbstractExpr) = DotSortAtom(a, b)

"""
    LinearAlgebra.eigmax(X::Convex.AbstractExpr)

The epigraph of the maximum eigen value of \$X\$.

## Examples

Applies to a single square matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = eigmax(x)
eigmin (convex; real)
└─ 2×2 real variable (id: 428…695)

julia> size(atom)
(1, 1)
```
"""
LinearAlgebra.eigmax(x::AbstractExpr) = EigMaxAtom(x)

"""
    LinearAlgebra.eigmin(X::Convex.AbstractExpr)

The hypograph of the minimum eigen value of \$X\$.

## Examples

Applies to a single square matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = eigmin(x)
eigmin (concave; real)
└─ 2×2 real variable (id: 428…695)

julia> size(atom)
(1, 1)
```
"""
LinearAlgebra.eigmin(x::AbstractExpr) = EigMinAtom(x)

"""
    entropy(x::Convex.AbstractExpr)

The hypograph of \$\\sum_i -x_i \\log x_i\$.

## Examples

Applies to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = entropy(x)
sum (concave; real)
└─ entropy (concave; real)
   └─ 3-element real variable (id: 901…778)

julia> size(atom)
(1, 1)
```
"""
entropy(x::AbstractExpr) = sum(EntropyAtom(x))

"""
    entropy_elementwise(x::Convex.AbstractExpr)

The hypograph of \$-x \\log x\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> entropy_elementwise(x)
entropy (concave; real)
└─ real variable (id: 172…395)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = entropy_elementwise(x)
entropy (concave; real)
└─ 3-element real variable (id: 140…126)

julia> size(atom)
(3, 1)
```
"""
entropy_elementwise(x::AbstractExpr) = EntropyAtom(x)

"""
    Base.exp(x::Convex.AbstractExpr)

The epigraph of \$e^x\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> exp(x)
exp (convex; positive)
└─ real variable (id: 103…720)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = exp(x)
exp (convex; positive)
└─ 3-element real variable (id: 389…882)

julia> size(atom)
(3, 1)
```
"""
Base.exp(x::AbstractExpr) = ExpAtom(x)

"""
    geomean(x::Convex.AbstractExpr...)

The hypograph of the geometric mean \$\\sqrt[n]{x_1 \\cdot x_2 \\cdot \\ldots x_n}\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> y = Variable();

julia> geomean(x, y)
geomean (concave; positive)
├─ real variable (id: 163…519)
└─ real variable (id: 107…393)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = Variable(3);

julia> atom = geomean(x, y)
geomean (concave; positive)
├─ 3-element real variable (id: 177…782)
└─ 3-element real variable (id: 307…913)

julia> size(atom)
(3, 1)
```
"""
geomean(args::Union{AbstractExpr,Value}...) = GeoMeanAtom(args...)

"""
    Base.hcat(args::AbstractExpr...)

Horizontally concatenate `args`.

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = hcat(x, x)
hcat (affine; real)
├─ 2×2 real variable (id: 111…376)
└─ 2×2 real variable (id: 111…376)

julia> size(atom)
(2, 4)
```

You can also use the Julia `[x x]` syntax:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = [x x]
hcat (affine; real)
├─ 2×2 real variable (id: 111…376)
└─ 2×2 real variable (id: 111…376)

julia> size(atom)
(2, 4)
```
"""
Base.hcat(args::AbstractExpr...) = HcatAtom(args...)

function Base.hcat(args::Union{AbstractExpr,Value}...)
    if all(Base.Fix2(isa, Value), args)
        return Base.cat(args..., dims = Val(2))
    end
    return HcatAtom(args...)
end

"""
    Base.hvcat(
        rows::Tuple{Vararg{Int}},
        args::Union{AbstractExpr,Value}...,
    )

Horizontally and vertically concatenate `args` in single call.

`rows` is the number of arguments to vertically concatenate into each column.

## Examples

Applies to a matrix:

To make the matrix:
```
a    b[1] b[2]
c[1] c[2] c[3]
```
do:
```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> a = Variable();

julia> b = Variable(1, 2);

julia> c = Variable(1, 3);

julia> atom = hvcat((2, 1), a, b, c)
vcat (affine; real)
├─ hcat (affine; real)
│  ├─ real variable (id: 429…021)
│  └─ 1×2 real variable (id: 120…326)
└─ hcat (affine; real)
   └─ 1×3 real variable (id: 124…615)

julia> size(atom)
(2, 3)
```
"""
function Base.hvcat(
    rows::Tuple{Vararg{Int}},
    args::Union{AbstractExpr,Value}...,
)
    output_rows = Vector{HcatAtom}(undef, length(rows))
    offset = 0
    for (i, n) in enumerate(rows)
        output_rows[i] = HcatAtom(args[offset.+(1:n)]...)
        offset += n
    end
    return vcat(output_rows...)
end

"""
    hinge_loss(x::Convex.AbstractExpr)

The epigraph of \$\\max(1 - x, 0)\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> hinge_loss(x)
max (convex; positive)
├─ + (affine; real)
│  ├─ [1;;]
│  └─ Convex.NegateAtom (affine; real)
│     └─ real variable (id: 129…000)
└─ [0;;]
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = hinge_loss(x)
max (convex; positive)
├─ + (affine; real)
│  ├─ * (constant; positive)
│  │  ├─ [1;;]
│  │  └─ [1.0; 1.0; 1.0;;]
│  └─ Convex.NegateAtom (affine; real)
│     └─ 3-element real variable (id: 125…591)
└─ [0;;]

julia> size(atom)
(3, 1)
```
"""
hinge_loss(x::AbstractExpr) = pos(1 - x)

"""
    huber(x::Convex.AbstractExpr, M::Real = 1.0)

The epigraph of the Huber loss function:
```math
\\begin{cases}
    x^2         & |x| \\le M \\\\
    2M|x| - M^2 & |x| > M
\\end{cases}
```
where \$M \\ge 1\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> huber(x, 2.5)
huber (convex; positive)
└─ real variable (id: 973…369)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = huber(x)
huber (convex; positive)
└─ 3-element real variable (id: 896…728)

julia> size(atom)
(3, 1)
```
"""
huber(x::AbstractExpr, M::Real = 1.0) = HuberAtom(x, M)

"""
    Base.imag(x::Convex.AbstractExpr)

Return the imaginary component of `x`.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = ComplexVariable();

julia> imag(x)
imag (affine; real)
└─ complex variable (id: 407…692)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = ComplexVariable(3);

julia> atom = imag(x)
imag (affine; real)
└─ 3-element complex variable (id: 435…057)

julia> size(atom)
(3, 1)
```
"""
Base.imag(x::AbstractExpr) = ImaginaryAtom(x)

Base.imag(x::ComplexConstant) = x.imag_constant

Base.imag(x::Constant) = Constant(zero(x.value))

"""
    inner_product(x::AbstractExpr, y::AbstractExpr)

The inner product \$tr(x^\\top y)\$ where `x` and `y` are square matrices.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> y = [1 3; 2 4];

julia> atom = inner_product(x, y)
real (affine; real)
└─ sum (affine; real)
   └─ diag (affine; real)
      └─ * (affine; real)
         ├─ …
         └─ …

julia> size(atom)
(1, 1)
```
"""
function inner_product(x::AbstractExpr, y::AbstractExpr)
    if !(x.size == y.size && x.size[1] == x.size[2])
        error("arguments must be square matrices of the same dimension")
    end
    return real(LinearAlgebra.tr(x' * y))
end

inner_product(x::Value, y::AbstractExpr) = inner_product(constant(x), y)

inner_product(x::AbstractExpr, y::Value) = inner_product(x, constant(y))

"""
    invpos(x::Convex.AbstractExpr)

The epigraph of \$\\frac{1}{x}\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> invpos(x)
qol_elem (convex; positive)
├─ [1.0;;]
└─ real variable (id: 139…839)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = invpos(x)
qol_elem (convex; positive)
├─ [1.0; 1.0; 1.0;;]
└─ 3-element real variable (id: 133…285)

julia> size(atom)
(3, 1)
```
"""
invpos(x::AbstractExpr) = QolElemAtom(constant(ones(x.size)), x)

"""
    Base.log(x::Convex.AbstractExpr)

The hypograph of \$\\log(x)\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> log(x)
log (concave; real)
└─ real variable (id: 103…720)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = log(x)
log (concave; real)
└─ 3-element real variable (id: 161…499)

julia> size(atom)
(3, 1)
```
"""
Base.log(x::AbstractExpr) = LogAtom(x)

"""
    log_perspective(x::Convex.AbstractExpr, y::Convex.AbstractExpr)

The hypograph the perspective of of the log function:
\$\\sum y_i*\\log \\frac{x_i}{y_i}\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> y = Variable();

julia> log_perspective(x, y)
Convex.NegateAtom (concave; real)
└─ relative_entropy (convex; real)
   ├─ real variable (id: 136…971)
   └─ real variable (id: 131…344)
```

And to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = Variable(3);

julia> atom = log_perspective(x, y)
Convex.NegateAtom (concave; real)
└─ relative_entropy (convex; real)
   ├─ 3-element real variable (id: 854…248)
   └─ 3-element real variable (id: 111…174)

julia> size(atom)
(1, 1)
```
"""
log_perspective(x::AbstractExpr, y::AbstractExpr) = -relative_entropy(y, x)

"""
    LinearAlgebra.logdet(X::Convex.AbstractExpr)

The hypograph of \$\\log(\\det(X))\$.

## Examples

Applies to a single matrix expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> X = Variable(2, 2);

julia> atom = logdet(X)
logdet (concave; real)
└─ 2×2 real variable (id: 159…883)

julia> size(atom)
(1, 1)
```
"""
LinearAlgebra.logdet(x::AbstractExpr) = LogDetAtom(x)

"""
    logisticloss(x::Convex.AbstractExpr)

Reformulation for epigraph of the logistic loss: \$\\sum_i \\log(e^x_i + 1)\$.

This reformulation uses [`logsumexp`](@ref).

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> logisticloss(x)
logsumexp (convex; real)
└─ vcat (affine; real)
   ├─ real variable (id: 444…892)
   └─ [0;;]
```

And to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = logisticloss(x)
+ (convex; real)
├─ logsumexp (convex; real)
│  └─ vcat (affine; real)
│     ├─ index (affine; real)
│     │  └─ …
│     └─ [0;;]
├─ logsumexp (convex; real)
│  └─ vcat (affine; real)
│     ├─ index (affine; real)
│     │  └─ …
│     └─ [0;;]
└─ logsumexp (convex; real)
   └─ vcat (affine; real)
      ├─ index (affine; real)
      │  └─ …
      └─ [0;;]

julia> size(atom)
(1, 1)
```
"""
function logisticloss(e::AbstractExpr)
    if length(e) == 1
        return logsumexp([e; 0])
    end
    return sum(logsumexp(hcat(vec(e), zeros(length(e))); dims = 2))
end

"""
    logsumexp(x::Convex.AbstractExpr)

The epigraph of \$\\log\\left(\\sum_i e^{x_i}\\right)\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = logsumexp(x)
logsumexp (convex; real)
└─ 3-element real variable (id: 102…634)

julia> size(atom)
(1, 1)
```
"""
logsumexp(x::AbstractExpr; dims = Colon()) = LogSumExpAtom(x, dims)

"""
    matrixfrac(x::AbstractExpr, P::AbstractExpr)

The epigraph of \$x^\\top P^{-1} x\$.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> P = Variable(2, 2);

julia> atom = matrixfrac(x, P)
matrixfrac (convex; positive)
├─ 2-element real variable (id: 139…388)
└─ 2×2 real variable (id: 126…414)

julia> size(atom)
(1, 1)
```
"""
matrixfrac(x::AbstractExpr, P::AbstractExpr) = MatrixFracAtom(x, P)

matrixfrac(x::Value, P::AbstractExpr) = MatrixFracAtom(constant(x), P)

matrixfrac(x::AbstractExpr, P::Value) = MatrixFracAtom(x, constant(P))

"""
    Base.max(x::Convex.AbstractExpr, y::Convex.AbstractExpr)
    Base.max(x::Convex.AbstractExpr, y::Convex.Value)
    Base.max(x::Convex.Value, y::Convex.AbstractExpr)

The hypograph of \$max(x, y)\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> max(x, 1)
max (convex; real)
├─ real variable (id: 183…974)
└─ [1;;]
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = [1, 2, 3];

julia> atom = max(x, y)
max (convex; real)
├─ 3-element real variable (id: 153…965)
└─ [1; 2; 3;;]

julia> size(atom)
(3, 1)
```
"""
Base.max(x::AbstractExpr, y::AbstractExpr) = MaxAtom(x, y)

Base.max(x::AbstractExpr, y::Value) = max(x, constant(y))

Base.max(x::Value, y::AbstractExpr) = max(constant(x), y)

"""
    Base.maximum(x::Convex.AbstractExpr)

The hypograph of \$max(x...)\$.

## Examples

Applies to a matrix expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = maximum(x)
maximum (convex; real)
└─ 3-element real variable (id: 159…219)

julia> size(atom)
(1, 1)
```
"""
Base.maximum(x::AbstractExpr) = MaximumAtom(x)

"""
    Base.min(x::Convex.AbstractExpr, y::Convex.AbstractExpr)
    Base.min(x::Convex.Value, y::Convex.AbstractExpr)
    Base.min(x::Convex.AbstractExpr, y::Convex.Value)

The epigraph of \$min(x, y)\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> min(x, 1)
min (concave; real)
├─ real variable (id: 183…974)
└─ [1;;]
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = [1, 2, 3];

julia> atom = min(x, y)
min (concave; real)
├─ 3-element real variable (id: 153…965)
└─ [1; 2; 3;;]

julia> size(atom)
(3, 1)
```
"""
Base.min(x::AbstractExpr, y::AbstractExpr) = MinAtom(x, y)

Base.min(x::AbstractExpr, y::Value) = min(x, constant(y))

Base.min(x::Value, y::AbstractExpr) = min(constant(x), y)

"""
    Base.minimum(x::Convex.AbstractExpr)

The epigraph of \$min(x...)\$.

## Examples

Applies to a matrix expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = minimum(x)
minimum (convex; real)
└─ 3-element real variable (id: 159…219)

julia> size(atom)
(1, 1)
```
"""
Base.minimum(x::AbstractExpr) = MinimumAtom(x)

"""
    neg(x::Convex.AbstractExpr)

The epigraph of \$\\max(-x, 0)\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> neg(x)
max (convex; positive)
├─ Convex.NegateAtom (affine; real)
│  └─ real variable (id: 467…111)
└─ [0;;]
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = neg(x)
max (convex; positive)
├─ Convex.NegateAtom (affine; real)
│  └─ 3-element real variable (id: 224…439)
└─ [0;;]

julia> size(atom)
(3, 1)
```
"""
neg(x::AbstractExpr) = max(-x, 0)

"""
    LinearAlgebra.norm2(x::Convex.AbstractExpr)

The epigraph of the 2-norm \$||x||_2\$.

## Examples

Applies to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = norm2(x)
norm2 (convex; positive)
└─ 3-element real variable (id: 162…975)

julia> size(atom)
(3, 1)
```

And to a complex:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> y = ComplexVariable(3);

julia> atom = norm2(y)
norm2 (convex; positive)
└─ vcat (affine; real)
   ├─ real (affine; real)
   │  └─ 3-element complex variable (id: 120…942)
   └─ imag (affine; real)
      └─ 3-element complex variable (id: 120…942)

julia> size(atom)
(1, 1)
```
"""
function LinearAlgebra.norm2(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return EuclideanNormAtom([real(x); imag(x)])
    end
    return EuclideanNormAtom(x)
end

"""
    nuclearnorm(x::Convex.AbstractExpr)

The epigraph of the nuclear norm \$||X||_*\$, which is the sum of the singular
values of \$X\$.

## Examples

Applies to a real-valued matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = nuclearnorm(x)
nuclearnorm (convex; positive)
└─ 2×2 real variable (id: 106…758)

julia> size(atom)
(1, 1)

julia> y = ComplexVariable(2, 2);

julia> atom = nuclearnorm(y)
nuclearnorm (convex; positive)
└─ 2×2 complex variable (id: 577…313)

julia> size(atom)
(1, 1)
```
"""
nuclearnorm(x::AbstractExpr) = NuclearNormAtom(x)

"""
    LinearAlgebra.opnorm(x::Convex.AbstractExpr, p::Real = 2)

The epigraph of the matrix norm \$||X||_p\$.

## Examples

Applies to a real- or complex-valued matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = LinearAlgebra.opnorm(x, 1)
maximum (convex; positive)
└─ * (convex; positive)
   ├─ [1.0 1.0]
   └─ abs (convex; positive)
      └─ 2×2 real variable (id: 106…758)

julia> atom = LinearAlgebra.opnorm(x, 2)
opnorm (convex; positive)
└─ 2×2 real variable (id: 106…758)

julia> atom = LinearAlgebra.opnorm(x, Inf)
maximum (convex; positive)
└─ * (convex; positive)
    ├─ abs (convex; positive)
    │  └─ 2×2 real variable (id: 106…758)
    └─ [1.0; 1.0;;]


julia> y = ComplexVariable(2, 2);

julia> atom = maximum (convex; positive)
└─ * (convex; positive)
   ├─ abs (convex; positive)
   │  └─ 2×2 complex variable (id: 116…943)
   └─ [1.0; 1.0;;]

julia> size(atom)
(1, 1)
```
"""
function LinearAlgebra.opnorm(x::AbstractExpr, p::Real = 2)
    if length(size(x)) <= 1 || minimum(size(x)) == 1
        throw(ArgumentError("argument to `opnorm` must be a matrix"))
    end
    if p == 1
        return maximum(sum(abs(x), dims = 1))
    elseif p == 2
        return OperatorNormAtom(x)
    elseif p == Inf
        return maximum(sum(abs(x), dims = 2))
    else
        throw(
            ArgumentError("matrix p-norms only defined for p = 1, 2, and Inf"),
        )
    end
end

"""
    pos(x::Convex.AbstractExpr)

The epigraph of \$\\max(x, 0)\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> pos(x)
max (convex; positive)
├─ real variable (id: 467…111)
└─ [0;;]
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = pos(x)
max (convex; positive)
├─ 3-element real variable (id: 154…809)
└─ [0;;]

julia> size(atom)
(3, 1)
```
"""
pos(x::AbstractExpr) = max(x, 0)

"""
    qol_elementwise(x::AbstractExpr, y::AbstractExpr)

The elementwise epigraph of \$\\frac{x^2}{y}\$.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = Variable(3, Positive());

julia> atom = qol_elementwise(x, y)
qol_elem (convex; positive)
├─ 3-element real variable (id: 155…648)
└─ 3-element positive variable (id: 227…080)

julia> size(atom)
(3, 1)
```
"""
qol_elementwise(x::AbstractExpr, y::AbstractExpr) = QolElemAtom(x, y)

"""
    quadoverlin(x::AbstractExpr, y::AbstractExpr)

The epigraph of \$\\frac{||x||_2^2}{y}\$.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = Variable(Positive());

julia> atom = quadoverlin(x, y)
qol (convex; positive)
├─ 3-element real variable (id: 868…883)
└─ positive variable (id: 991…712)

julia> size(atom)
(1, 1)
```
"""
quadoverlin(x::AbstractExpr, y::AbstractExpr) = QuadOverLinAtom(x, y)

"""
    quantum_entropy(X::AbstractExpr, m::Integer, k::Integer)

`quantum_entropy` returns `-LinearAlgebra.tr(X*log(X))` where `X` is a positive
semidefinite.

Note this function uses logarithm base e, not base 2, so return value is in
units of nats, not bits.

Quantum entropy is concave. This function implements the semidefinite
programming approximation given in the reference below.  Parameters `m` and `k`
control the accuracy of this approximation: `m` is the number of quadrature
nodes to use and `k` the number of square-roots to take. See reference for more
details.

The implementation uses the expression
```math
H(X) = -tr(D_{op}(X||I))
```
where \$D_{op}\$ is the operator relative entropy:
```math
D_{op}(X||Y) = X^{1/2}*logm(X^{1/2} Y^{-1} X^{1/2})*X^{1/2}
```

## Reference

Ported from CVXQUAD which is based on the paper: "Lieb's concavity theorem,
matrix geometric means and semidefinite optimization" by Hamza Fawzi and James
Saunderson (arXiv:1512.03401)

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> X = Variable(2, 2);

julia> atom = quantum_entropy(X)
quantum_entropy (concave; positive)
└─ 2×2 real variable (id: 700…694)

julia> size(atom)
(1, 1)
```
"""
function quantum_entropy(X::AbstractExpr, m::Integer = 3, k::Integer = 3)
    return QuantumEntropyAtom(X, m, k)
end

function quantum_entropy(
    X::Union{AbstractMatrix,Constant},
    m::Integer = 0,
    k::Integer = 0,
)
    return -quantum_relative_entropy(X, Matrix(1.0 * LinearAlgebra.I, size(X)))
end

"""
    quantum_relative_entropy(
        A::AbstractExpr,
        B::AbstractExpr,
        m::Integer,
        k::Integer,
    )

`quantum_relative_entropy` returns `LinearAlgebra.tr(A*(log(A)-log(B)))` where
`A` and `B` are positive semidefinite matrices.

Note this function uses logarithm base e, not base 2, so return value is in
units of nats, not bits.

Quantum relative entropy is convex (jointly) in (A, B). This function implements
the semidefinite programming approximation given in the reference below.
Parameters m and k control the accuracy of this approximation: `m` is the number
of quadrature nodes to use and `k` the number of square-roots to take. See
reference for more details.

Implementation uses the expression
```math
D(A||B) = e'*D_{op} (A \\otimes I || I \\otimes B) )*e
```
where \$D_{op}\$ is the operator relative entropy and `e = vec(Matrix(I, n, n))`.

## Reference

Ported from CVXQUAD which is based on the paper: "Lieb's concavity theorem,
matrix geometric means and semidefinite optimization" by Hamza Fawzi and James
Saunderson (arXiv:1512.03401)

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> A = Variable(2, 2);

julia> B = Variable(2, 2);

julia> atom = quantum_relative_entropy(A, B)
quantum_relative_entropy (convex; positive)
├─ 2×2 real variable (id: 144…849)
└─ 2×2 real variable (id: 969…693)

julia> size(atom)
(1, 1)
```
"""
function quantum_relative_entropy(
    A::AbstractExpr,
    B::AbstractExpr,
    m::Integer = 3,
    k::Integer = 3,
)
    return QuantumRelativeEntropy1Atom(A, B, m, k)
end

function quantum_relative_entropy(
    A::AbstractExpr,
    B::Union{AbstractMatrix,Constant},
    m::Integer = 3,
    k::Integer = 3,
    nullspace_tol::Real = 1e-6,
)
    # This `evaluate` is safe since it is not a `fix!`ed variable
    # (it must be a constant or matrix)
    return QuantumRelativeEntropy2Atom(A, evaluate(B), m, k, nullspace_tol)
end

function quantum_relative_entropy(
    A::Union{AbstractMatrix,Constant},
    B::AbstractExpr,
    m::Integer = 3,
    k::Integer = 3,
)
    # This `evaluate` is safe since it is not a `fix!`ed variable
    # (it must be a constant or matrix)
    A = evaluate(A)
    return -quantum_entropy(A, m, k) - trace_logm(B, A, m, k)
end

function quantum_relative_entropy(
    A::Union{AbstractMatrix,Constant},
    B::Union{AbstractMatrix,Constant},
    m::Integer = 0,
    k::Integer = 0,
    nullspace_tol::Real = 1e-6,
)
    # These `evaluate`s are safe since it is not a `fix!`ed variable
    # (it must be a constant or matrix)
    A, B = evaluate(A), evaluate(B)
    if size(A) != size(B)
        throw(DimensionMismatch("A and B must be the same size"))
    elseif size(A) != (size(A)[1], size(A)[1])
        throw(DimensionMismatch("A and B must be square"))
    elseif norm(A - A') > nullspace_tol
        throw(DomainError(A, "A must be Hermitian"))
    elseif norm(B - B') > nullspace_tol
        throw(DomainError(B, "B must be Hermitian"))
    end
    v, U = LinearAlgebra.eigen(LinearAlgebra.Hermitian(A))
    if any(v .< -nullspace_tol)
        throw(DomainError(A, "A must be positive semidefinite"))
    end
    if any(LinearAlgebra.eigvals(LinearAlgebra.Hermitian(B)) .< -nullspace_tol)
        throw(DomainError(B, "B must be positive semidefinite"))
    end
    J = U'[v.>nullspace_tol, :]
    Ap = LinearAlgebra.Hermitian(J * A * J')
    Bp = LinearAlgebra.Hermitian(J * B * J')
    if any(LinearAlgebra.eigvals(Bp) .< nullspace_tol)
        return Inf
    end
    return real(LinearAlgebra.tr(Ap * (log(Ap) - log(Bp))))
end

"""
    rationalnorm(x::AbstractExpr, k::Rational{Int})

The epigraph of `||x||_k`.

## Examples

Applies to a single matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> atom = rationalnorm(x, 3 // 2)
rationalnorm (convex; positive)
└─ 2-element real variable (id: 182…293)

julia> size(atom)
(1, 1)
```
"""
function rationalnorm(x::AbstractExpr, k::Rational{Int})
    if sign(x) == ComplexSign()
        row, col = size(x)
        if !(row == 1 || col == 1)
            error("[RationalNormAtom] not defined for complex matrices")
        end
        return RationalNormAtom(abs(x), k)
    end
    return RationalNormAtom(x, k)
end

"""
    Base.real(x::Convex.AbstractExpr)

Return the real component of `x`.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = ComplexVariable();

julia> real(x)
real (affine; real)
└─ complex variable (id: 407…692)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = ComplexVariable(3);

julia> atom = real(x)
real (affine; real)
└─ 3-element complex variable (id: 435…057)

julia> size(atom)
(3, 1)
```
"""
Base.real(x::AbstractExpr) = RealAtom(x)

Base.real(x::ComplexConstant) = x.real_constant

Base.real(x::Constant) = x

"""
    relative_entropy(x::Convex.AbstractExpr, y::Convex.AbstractExpr)

The epigraph of \$\\sum y_i*\\log \\frac{x_i}{y_i}\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> y = Variable();

julia> relative_entropy(x, y)
relative_entropy (convex; real)
├─ real variable (id: 124…372)
└─ real variable (id: 409…346)
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> y = Variable(3);

julia> atom = relative_entropy(x, y)
julia> atom = relative_entropy(x, y)
relative_entropy (convex; real)
├─ 3-element real variable (id: 906…671)
└─ 3-element real variable (id: 118…912)

julia> size(atom)
(1, 1)
```
"""
relative_entropy(x::AbstractExpr, y::AbstractExpr) = RelativeEntropyAtom(x, y)

"""
    Base.reshape(x::AbstractExpr, m::Int, n::Int)

Reshapes the expression `x` into a matrix with `m` rows and `n` columns.

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(6, 1);

julia> size(x)
(6, 1)

julia> atom = reshape(x, 2, 3)
reshape (affine; real)
└─ 6-element real variable (id: 103…813)

julia> size(atom)
(2, 3)
```
"""
Base.reshape(x::AbstractExpr, m::Int, n::Int) = ReshapeAtom(x, m, n)

"""
    Convex.rootdet(X::Convex.AbstractExpr)

The hypograph of \$\\det(X)^{\\frac{1}{n}}\$, where \$n\$ is the side-dimension
of the square matrix \$X\$.

## Examples

Applies to a single matrix expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> X = Variable(2, 2);

julia> atom = rootdet(X)
rootdet (concave; real)
└─ 2×2 real variable (id: 159…883)

julia> size(atom)
(1, 1)
```
"""
rootdet(x::AbstractExpr) = RootDetAtom(x)

"""
    sigmamax(x::Convex.AbstractExpr)

The epigraph of the spectral norm \$||X||_2\$, which is the maximum of the
singular values of \$X\$.

## Examples

Applies to a real- or complex-valued matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = sigmamax(x)
opnorm (convex; positive)
└─ 2×2 real variable (id: 106…758)

julia> size(atom)
(1, 1)

julia> y = ComplexVariable(2, 2);

julia> atom = sigmamax(y)
opnorm (convex; positive)
└─ 2×2 complex variable (id: 577…313)

julia> size(atom)
(1, 1)
```
"""
sigmamax(x::AbstractExpr) = OperatorNormAtom(x)

"""
    Base.sqrt(x::Convex.AbstractExpr)

The hypograph of \$\\sqrt x\$.

## Examples

Applies to a single expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable();

julia> sqrt(x)
geomean (concave; positive)
├─ real variable (id: 576…546)
└─ [1.0;;]
```

And element-wise to a matrix of expressions:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = sqrt(x)
geomean (concave; positive)
├─ 3-element real variable (id: 181…583)
└─ [1.0; 1.0; 1.0;;]

julia> size(atom)
(3, 1)
```
"""
function Base.sqrt(x::AbstractExpr)
    return GeoMeanAtom(x, constant(ones(x.size[1], x.size[2])))
end

"""
    Base.sum(x::Convex.AbstractExpr; dims = :)

Sum `x`, optionally along a dimension `dims`.

## Examples

Sum all elements in an expression:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = sum(x)
sum (affine; real)
└─ 2×2 real variable (id: 263…449)

julia> size(atom)
(1, 1)
```

Sum along the first dimension, creating a row vector:
```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = sum(x; dims = 1)
* (affine; real)
├─ [1.0 1.0]
└─ 2×2 real variable (id: 143…826)

julia> size(atom)
(1, 2)
```

Sum along the second dimension, creating a columnn vector:
```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> atom = sum(x; dims = 2)
* (affine; real)
├─ 2×2 real variable (id: 143…826)
└─ [1.0; 1.0;;]

julia> size(atom)
(2, 1)
```
"""
Base.sum(x::AbstractExpr; dims = :) = _sum(x, dims)

_sum(x::AbstractExpr, ::Colon) = SumAtom(x)

function _sum(x::AbstractExpr, dimension::Integer)
    if dimension == 1
        return Constant(ones(1, x.size[1]), Positive()) * x
    elseif dimension == 2
        return x * Constant(ones(x.size[2], 1), Positive())
    else
        error("[SumAtom] sum not implemented for `dims=$dimension`")
    end
end

"""
    sumlargesteigs(x::Convex.AbstractExpr, k::Int)

Sum the `k` largest eigen values of `x`.

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3, 3);

julia> atom = sumlargesteigs(x, 2)
sumlargesteigs (convex; real)
├─ 3×3 real variable (id: 833…482)
└─ [2;;]

julia> size(atom)
(1, 1)
```
"""
function sumlargesteigs(x::AbstractExpr, k::Int)
    return k == 0 ? Constant(0) : SumLargestEigsAtom(x, Constant(k))
end

"""
    sumlargest(x::Convex.AbstractExpr, k::Int)

Sum the `k` largest values of `x`.

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3, 3);

julia> atom = sumlargest(x, 2)
sumlargest (convex; real)
├─ 3×3 real variable (id: 833…482)
└─ [2;;]

julia> size(atom)
(1, 1)
```
"""
function sumlargest(x::AbstractExpr, k::Int)
    return k == 0 ? Constant(0) : SumLargestAtom(x, k)
end

"""
    sumsmallest(x::Convex.AbstractExpr, k::Int)

Sum the `k` smallest values of `x`.

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3, 3);

julia> atom = sumsmallest(x, 2)
Convex.NegateAtom (concave; real)
└─ sumlargest (convex; real)
   └─ Convex.NegateAtom (affine; real)
      └─ 3×3 real variable (id: 723…082)

julia> size(atom)
(1, 1)
```
"""
function sumsmallest(x::AbstractExpr, k::Int)
    return k == 0 ? Constant(0) : -SumLargestAtom(-x, k)
end

"""
    square(x::AbstractExpr)

The epigraph of \$x^2\$.

## Examples

Applies elementwise to a matrix

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = square(x)
qol_elem (convex; positive)
├─ 3-element real variable (id: 438…681)
└─ [1.0; 1.0; 1.0;;]

julia> size(atom)
(3, 1)
```
"""
function square(x::AbstractExpr)
    if sign(x) == ComplexSign()
        error(
            "square of a complex number is not DCP. Did you mean square_modulus?",
        )
    end
    return QolElemAtom(x, constant(ones(x.size)))
end

"""
    sumsquares(x::AbstractExpr)

The epigraph of \$||x||_2^2\$.

## Examples

Applies to a single matrix

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(3);

julia> atom = sumsquares(x)
qol (convex; positive)
├─ 3-element real variable (id: 125…181)
└─ [1;;]

julia> size(atom)
(1, 1)
```
"""
function sumsquares(x::AbstractExpr)
    if size(x, 2) != 1
        return QuadOverLinAtom(reshape(x, length(x), 1), constant(1))
    end
    return QuadOverLinAtom(x, constant(1))
end

"""
    LinearAlgebra.tr(x::AbstractExpr)

The trace of the matrix `x`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = tr(x)
sum (affine; real)
└─ diag (affine; real)
   └─ 2×2 real variable (id: 844…180)

julia> size(atom)
(1, 1)
```
"""
LinearAlgebra.tr(e::AbstractExpr) = sum(LinearAlgebra.diag(e))

"""
    trace_logm(
        X::Convex.AbstractExpr,
        C::AbstractMatrix,
        m::Integer = 3,
        k::Integer = 3,
    )

`trace_logm(X, C)` returns `LinearAlgebra.tr(C*logm(X))` where `X` and `C` are
positive definite matrices and `C` is constant.

`trace_logm` is concave in `X`.

This function implements the semidefinite programming approximation given in the
reference below. Parameters `m` and `k` control the accuracy of the
approximation: `m` is the number of quadrature nodes to use and `k` is the
number of square-roots to take. See reference for more details.

Implementation uses the expression

```math
tr(C \\times logm(X)) = -tr(C \\times D_{op}(I||X))
```
where `D_{op}` is the operator relative entropy:
```math
D_{op}(X||Y) = X^{1/2}*logm(X^{1/2} Y^{-1} X^{1/2})*X^{1/2}
```

## Reference

Ported from CVXQUAD which is based on the paper: "Lieb's concavity theorem,
matrix geometric means and semidefinite optimization" by Hamza Fawzi and James
Saunderson (arXiv:1512.03401)

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> X = Variable(2, 2);

julia> C = [1 0; 0 1];

julia> atom = trace_logm(X, C)
trace_logm (concave; real)
└─ 2×2 real variable (id: 608…362)

julia> size(atom)
(1, 1)
```
"""
function trace_logm(
    X::AbstractExpr,
    C::Union{AbstractMatrix,Constant},
    m::Integer = 3,
    k::Integer = 3,
)
    # This `evaluate` is safe since it is not a `fix!`ed variable
    # (it must be a constant or matrix)
    return TraceLogmAtom(X, evaluate(C), m, k)
end

function trace_logm(
    X::Union{AbstractMatrix,Constant},
    C::Union{AbstractMatrix,Constant},
    m::Integer = 3,
    k::Integer = 3,
)
    I = Matrix(1.0 * LinearAlgebra.I(size(X, 1)))
    return -quantum_relative_entropy(C, X) + quantum_relative_entropy(C, I)
end

"""
    trace_mpower(A::Convex.AbstractExpr, t::Rational, C::AbstractMatrix)

`trace_mpower(A, t, C)` returns `LinearAlgebra.tr(C*A^t)` where `A` and `C` are
positive definite matrices, `C` is constant and `t` is a rational in `[-1, 2]`.

When `t` is in `[0, 1]`, `trace_mpower(A, t, C)` is concave in `A` (for fixed
positive semidefinite matrix `C`) and convex for `t` in `[-1, 0)` or `(1, 2]`.

## Reference

Ported from CVXQUAD which is based on the paper: "Lieb's concavity theorem,
matrix geometric means and semidefinite optimization" by Hamza Fawzi and James
Saunderson (arXiv:1512.03401)

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> A = Variable(2, 2);

julia> C = [1 0; 0 1];

julia> atom = trace_mpower(A, 1 // 2, C)
trace_mpower (concave; real)
└─ 2×2 real variable (id: 150…626)

julia> size(atom)
(1, 1)
```
"""
function trace_mpower(
    A::AbstractExpr,
    t::Rational,
    C::Union{AbstractMatrix,Constant},
)
    # This `evaluate` is safe since it is not a `fix!`ed variable (it must be a
    # constant or a matrix).
    return TraceMpowerAtom(A, t, evaluate(C))
end

function trace_mpower(
    A::Union{AbstractMatrix,Constant},
    t::Rational,
    C::Union{AbstractMatrix,Constant},
)
    return LinearAlgebra.tr(C * A^t)
end

function trace_mpower(
    A::Union{AbstractExpr,Value},
    t::Integer,
    C::Union{AbstractMatrix,Constant},
)
    return trace_mpower(A, t // 1, C)
end

"""
    LinearAlgebra.transpose(x::AbstractExpr)

The transpose of the matrix `x`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = transpose(x)
reshape (affine; real)
└─ * (affine; real)
   ├─ 4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 4 stored entries
   └─ reshape (affine; real)
      └─ 2×2 real variable (id: 151…193)

julia> size(atom)
(2, 2)
```
"""
function LinearAlgebra.transpose(x::AbstractExpr)
    P = permutedims_matrix(size(x), (2, 1))
    return reshape(P * vec(x), size(x, 2), size(x, 1))
end

# These `evaluate` calls are safe since they are not a `fix!`ed variable:
function LinearAlgebra.transpose(x::Union{Constant,ComplexConstant})
    return constant(transpose(evaluate(x)))
end

"""
    Base.vcat(args::AbstractExpr...)

Vertically concatenate `args`.

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = vcat(x, x)
vcat (affine; real)
├─ 2×2 real variable (id: 111…376)
└─ 2×2 real variable (id: 111…376)

julia> size(atom)
(4, 2)
```

You can also use the Julia `[x; x]` syntax:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = [x; x]
vcat (affine; real)
├─ 2×2 real variable (id: 111…376)
└─ 2×2 real variable (id: 111…376)

julia> size(atom)
(4, 2)
```
"""
Base.vcat(args::AbstractExpr...) = VcatAtom(args...)

function Base.vcat(args::Union{AbstractExpr,Value}...)
    if all(Base.Fix2(isa, Value), args)
        return Base.cat(args..., dims = Val(1))
    end
    return VcatAtom(args...)
end

"""
    Base.vec(x::AbstractExpr)

Reshapes the expression `x` into a column vector.

## Examples

Applies to a matrix:

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2, 2);

julia> atom = vec(x)
reshape (affine; real)
└─ 2×2 real variable (id: 115…295)

julia> size(atom)
(4, 1)
```
"""
Base.vec(x::AbstractExpr) = reshape(x, length(x), 1)
