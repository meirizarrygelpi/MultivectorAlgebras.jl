"""
    Multivector2{T <: Real} <: AbstractMultivector{T}

An immutable pair of pairs of real numbers that represents a member
of a 2-dimensional multivector algebra.

Each `Multivector2` has the form
```math
    a+bW+cX+dWX
```
where ``a``, ``b``, ``c``, and ``d`` are real (and of the same type),
and ``W * W = 0``, ``X * X = 0``, and ``WX = W * X = -X * W``.
Here ``*`` is the wedge product.
"""
struct Multivector2{T <: Real} <: AbstractMultivector{T}
    l::Multivector1{T}
    r::Multivector1{T}

    Multivector2{U}(l::Multivector1{U}, r::Multivector1{U}) where {U <: Real} = new(l, r)
end

function Multivector2(x::Multivector1{T}, y::Multivector1{T}) where {T <: Real}
    Multivector2{T}(x, y)
end

function Multivector2(z::Multivector1)
    Multivector2(z, zero(z))
end

function Multivector2(a::T, b::T, c::T, d::T) where {T <: Real}
    Multivector2{T}(Multivector1(a, b), Multivector1(c, d))
end

function Multivector2(a::T, b::T, c::T) where {T <: Real}
    Multivector2{T}(Multivector1(a, b), Multivector1(c, zero(T)))
end

function Multivector2(a::T, b::T) where {T <: Real}
    Multivector2{T}(Multivector1(a, b), zero(Multivector1{T}))
end

function Multivector2(a::T) where {T <: Real}
    Multivector2{T}(Multivector1(a, zero(T)), zero(Multivector1{T}))
end

function Multivector2(a::Real, b::Real, c::Real, d::Real)
    Multivector2(promote(a, b, c, d)...)
end

function Multivector2(a::Real, b::Real, c::Real)
    Multivector2(promote(a, b, c)...)
end

function Multivector2(a::Real, b::Real)
    Multivector2(promote(a, b)...)
end

function Multivector2(x::Multivector1, y::Multivector1)
    Multivector2(promote(x.l, x.r, y.l, y.r)...)
end

function show(io::IO, z::Multivector2)
    print(io, "[1: ")
    print(io, z.l.l)
    print(io, ", W: ")
    print(io, z.l.r)
    print(io, ", X: ")
    print(io, z.r.l)
    print(io, ", WX: ")
    print(io, z.r.r)
    print(io, "]")
end

function real(z::Multivector2)
    real(z.l)
end

function unreal(z::Multivector2)
    vcat(unreal(z.l), asarray(z.r))
end

function zero(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(zero(z.l), zero(z.r))
end

function zero(::Type{Multivector2{T}}) where {T <: Real}
    Multivector2{T}(zero(T), zero(T), zero(T), zero(T))
end

function one(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(one(z.l), zero(z.r))
end

function one(::Type{Multivector2{T}}) where {T <: Real}
    Multivector2{T}(one(T), zero(T), zero(T), zero(T))
end

"""
    conj{T <: Real}(z::Multivector2{T})

The `Multivector2` conjugate. If ``z=a+bW+cX+dWX``, then `conj(z)` gives
```math
    a-bW-cX-dWX
```
This operation is an involution.
"""
function conj(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(conj(z.l), -z.r)
end

"""
    cloak{T <: Real}(z::Multivector2{T})

The cloak conjugate changes the sign of even blades. If ``z=a+bW+cX+dWX``, then `cloak(z)` gives
```math
    -a+bW+cX-dWX
```
This operation is equivalent to `-dagger(z)` and thus is also an involution.
"""
function cloak(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(cloak(z.l), dagger(z.r))
end

"""
    dagger{T <: Real}(z::Multivector2{T})

The dagger conjugate changes the sign of odd blades. If ``z=a+bW+cX+dWX``, then `dagger(z)` gives
```math
    a-bW-cX+dWX
```
This operation is an involution.
"""
function dagger(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(dagger(z.l), cloak(z.r))
end

"""
    star{T <: Real}(z::Multivector2{T})

Returns the Hodge star conjugate. If ``z=a+bW+cX+dWX``, then `star(z)` gives
```math
    d-cW+bX+aWX
```
This operation is an involution.
"""
function star(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(dagger(star(z.r)), star(z.l))
end

"""
    selfstar{T <: Real}(z::Multivector2{T})

The self-star-conjugate part. If ``z=a+bW+cX+dWX``, then `selfstar(z)` gives
```math
    \\frac{1}{2}(a+d)(1+WX)
```
This operation is not idempotent.
"""
function selfstar(z::Multivector2{T}) where {T <: Real}
    (z + star(z)) / 2
end

"""
    antiselfstar{T <: Real}(z::Multivector2{T})

The anti-self-star-conjugate part. If ``z=a+bW+cX+dWX``, then `antiselfstar(z)` gives
```math
    \\frac{1}{2}(a-d)(1-WX)
```
This operation is not idempotent.
"""
function antiselfstar(z::Multivector2{T}) where {T <: Real}
    AntiSelfStar2{T}((z.l.l - z.r.r) / 2)
end

(+)(z::Multivector2{T}) where {T <: Real} = z

function (+)(x::Multivector2, y::Multivector2)
    Multivector2(x.l + y.l, x.r + y.r)
end

function (+)(x::Multivector2, y::Multivector1)
    Multivector2(x.l + y, x.r)
end

function (+)(x::Multivector1, y::Multivector2)
    Multivector2(x + y.l, y.r)
end

function (+)(z::Multivector2, a::Real)
    Multivector2(z.l + a, z.r)
end

function (+)(a::Real, z::Multivector2)
    Multivector2(z.l + a, z.r)
end

function (-)(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(-z.l, -z.r)
end

function (-)(x::Multivector2, y::Multivector2)
    Multivector2(x.l - y.l, x.r - y.r)
end

function (-)(x::Multivector2, y::Multivector1)
    Multivector2(x.l - y, x.r)
end

function (-)(x::Multivector1, y::Multivector2)
    Multivector2(x - y.l, -y.r)
end

function (-)(z::Multivector2, a::Real)
    Multivector2(z.l - a, z.r)
end

function (-)(a::Real, z::Multivector2)
    Multivector2(a - z.l, -z.r)
end

"""
    (*)(x::Multivector2, y::Multivector2)

Wedge product of two 2-dimensional multivectors.
"""

function (*)(x::Multivector2, y::Multivector2)
    Multivector2(x.l * y.l,  (y.r * x.l) + (x.r * conj(y.l)))
end

function (*)(x::Multivector2, y::Multivector1)
    Multivector2(x.l * y, x.r * conj(y))
end

function (*)(x::Multivector1, y::Multivector2)
    Multivector2(x * y.l,  y.r * x)
end

"""
    (*)(z::Multivector2, a::Real)
    (*)(a::Real, z::Multivector2)

Scaling and/or reflection of a `Multivector2` by a real number.
"""
function (*)(z::Multivector2, a::Real)
    Multivector2(z.l * a, z.r * a)
end

function (*)(a::Real, z::Multivector2)
    Multivector2(z.l * a, z.r * a)
end

function commutator(x::AbstractMultivector, y::AbstractMultivector)
    (x * y) - (y * x)
end

"""
    abs(z::Multivector2)
"""
abs(z::Multivector2) = abs(z.l)

"""
    abs2(z::Multivector2)
"""
abs2(z::Multivector2) = abs2(z.l)

"""
    iszerodivisor(z::Multivector2)

Returns true if `z` is of the form ``bW+cX+dWX`` such that ``z * conj(z) = 0``.
"""
iszerodivisor(z::Multivector2{T}) where {T <: Real} = iszerodivisor(z.l)

"""
    inv(z::Multivector2)

Inverse of a `Multivector2`. An error occurs if argument is a zero divisor.
"""
function inv(z::Multivector2)
    if iszerodivisor(z)
        error(ZeroDivisorInverse)
    end

    conj(z) / abs2(z)
end

function (/)(x::Multivector2, y::Multivector2)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    x * inv(y)
end

function (/)(a::Real, z::Multivector2)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    a * inv(z)
end

function (/)(z::Multivector2, a::Real)
    if iszero(a)
        error(ZeroDenominator)
    end

    Multivector2(z.l / a, z.r / a)
end

function (\)(y::Multivector2, x::Multivector2)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    inv(y) * x
end

function (\)(z::Multivector2, a::Real)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    a * inv(z)
end

function (\)(a::Real, z::Multivector2)
    if iszero(a)
        error(ZeroDenominator)
    end
    
    Multivector2(z.l / a, z.r / a)
end

"""
    crossratioL(w::AbstractMultivector,
               x::AbstractMultivector,
               y::AbstractMultivector,
               z::AbstractMultivector)

The left cross-ratio:
```julia
    ((x-y) \\ (w-y)) * ((w-z) \\ (x-z))
```
The left cross-ratio is invariant under right Möbius transformations.
"""
function crossratioL(w::AbstractMultivector,
                     x::AbstractMultivector,
                     y::AbstractMultivector,
                     z::AbstractMultivector)
    (inv(x - y) * (w - y)) * (inv(w - z) * (x - z))
end

"""
    crossratioR(w::AbstractMultivector,
                x::AbstractMultivector,
                y::AbstractMultivector,
                z::AbstractMultivector)

The right cross-ratio:
```julia
    ((w-y) / (x-y)) * ((x-z) / (w-z))
```
The right cross-ratio is invariant under left Möbius transformations.
"""
function crossratioR(w::AbstractMultivector,
                     x::AbstractMultivector,
                     y::AbstractMultivector,
                     z::AbstractMultivector)
    ((w - y) * inv(x - y)) * ((x - z) * inv(w - z))
end

"""
    möbiusL(z::AbstractMultivector,
           a::AbstractMultivector,
           b::AbstractMultivector,
           c::AbstractMultivector,
           d::AbstractMultivector)
    möbiusL(z::AbstractMultivector, a::Real, b::Real, c::Real, d::Real)

The left Möbius transformation:
```julia
    ((z * c) + d) \\ ((z * a) + b)
```
This transformation is also know as a fractional linear transformation.
"""
function möbiusL(z::AbstractMultivector,
                a::AbstractMultivector,
                b::AbstractMultivector,
                c::AbstractMultivector,
                d::AbstractMultivector)
    inv((z * c) + d) * ((z * a) + b)
end

function möbiusL(z::AbstractMultivector, a::Real, b::Real, c::Real, d::Real)
    inv((c * z) + d) * ((a * z) + b)
end

"""
    möbiusR(z::AbstractMultivector,
           a::AbstractMultivector,
           b::AbstractMultivector,
           c::AbstractMultivector,
           d::AbstractMultivector)
    möbiusR(z::AbstractMultivector, a::Real, b::Real, c::Real, d::Real)

The right Möbius transformation:
```julia
    ((a * z) + b) / ((c * z) + d)
```
This transformation is also know as a fractional linear transformation.
"""
function möbiusR(z::AbstractMultivector,
                a::AbstractMultivector,
                b::AbstractMultivector,
                c::AbstractMultivector,
                d::AbstractMultivector)
    ((a * z) + b) * inv((c * z) + d)
end

function möbiusR(z::AbstractMultivector, a::Real, b::Real, c::Real, d::Real)
    ((a * z) + b) * inv((c * z) + d)
end
