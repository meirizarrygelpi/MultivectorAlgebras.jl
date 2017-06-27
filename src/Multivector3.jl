"""
    AbstractMultivector3{T <: Real} <: AbstractMultivector{T}

An abstract 3-dimensional multivector.
"""
abstract type AbstractMultivector3{T <: Real} <: AbstractMultivector{T} end

"""
    Multivector3{T <: Real} <: AbstractMultivector3{T}

An immutable pair of pairs of real numbers that represents a member
of a 3-dimensional multivector algebra.

Each `Multivector3` has the form
```math
    a+bW+cX+dWX+fY+gWY+hXY+j(WX)Y
```
where ``a``, ``b``, ``c``, ``d``, ``f``, ``g``, ``h``, and ``j`` are real
(and of the same type), and ``W ∧ W = 0``, ``X ∧ X = 0``, `` Y ∧ Y = 0``;
``WX = W ∧ X = -X ∧ W``, ``WY = W ∧ Y = -Y ∧ W``, and ``XY = X ∧ Y = -Y ∧ X``.
Here ``∧`` is the wedge product.
"""
struct Multivector3{T <: Real} <: AbstractMultivector3{T}
    l::Multivector2{T}
    r::Multivector2{T}
    
    Multivector3{U}(l::Multivector2{U}, r::Multivector2{U}) where {U <: Real} = new(l, r)
end

function Multivector3(x::Multivector2{T}, y::Multivector2{T}) where {T <: Real}
    Multivector3{T}(x, y)
end

function Multivector3(a::T, b::T, c::T, d::T, f::T, g::T, h::T, j::T) where {T <: Real}
    Multivector3{T}(Multivector2(a, b, c, d), Multivector2(f, g, h, j))
end

function Multivector3(a::T, b::T, c::T, d::T, f::T, g::T, h::T) where {T <: Real}
    Multivector3{T}(Multivector2(a, b, c, d), Multivector2(f, g, h))
end

function Multivector3(a::T, b::T, c::T, d::T, f::T, g::T) where {T <: Real}
    Multivector3{T}(Multivector2(a, b, c, d), Multivector2(f, g))
end

function Multivector3(a::T, b::T, c::T, d::T, f::T) where {T <: Real}
    Multivector3{T}(Multivector2(a, b, c, d), Multivector2(f))
end

function Multivector3(a::T, b::T, c::T, d::T) where {T <: Real}
    Multivector3{T}(Multivector2(a, b, c, d), zero(Multivector2{T}))
end

function Multivector3(a::T, b::T, c::T) where {T <: Real}
    Multivector3{T}(Multivector2(a, b, c), zero(Multivector2{T}))
end

function Multivector3(a::T, b::T) where {T <: Real}
    Multivector3{T}(Multivector2(a, b), zero(Multivector2{T}))
end

function Multivector3(a::T) where {T <: Real}
    Multivector3{T}(Multivector2(a), zero(Multivector2{T}))
end

function Multivector3(a::Real, b::Real, c::Real, d::Real, f::Real, g::Real, h::Real, j::Real)
    Multivector3(promote(a, b, c, d, f, g, h, j)...)
end

function Multivector3(a::Real, b::Real, c::Real, d::Real, f::Real, g::Real, h::Real)
    Multivector3(promote(a, b, c, d, f, g, h)...)
end

function Multivector3(a::Real, b::Real, c::Real, d::Real, f::Real, g::Real)
    Multivector3(promote(a, b, c, d, f, g)...)
end

function Multivector3(a::Real, b::Real, c::Real, d::Real, f::Real)
    Multivector3(promote(a, b, c, d, f)...)
end

function Multivector3(a::Real, b::Real, c::Real, d::Real)
    Multivector3(promote(a, b, c, d)...)
end

function Multivector3(a::Real, b::Real, c::Real)
    Multivector3(promote(a, b, c)...)
end

function Multivector3(a::Real, b::Real)
    Multivector3(promote(a, b)...)
end

function Multivector3(x::Multivector2, y::Multivector2)
    Multivector3(promote(asarray(x)..., asarray(y)...)...)
end

function Multivector3(x::Multivector2)
    Multivector3(x, zero(x))
end

function Multivector3(r::Multivector1, s::Multivector1, t::Multivector1, u::Multivector1)
    Multivector3(Multivector2(r, s), Multivector2(t, u))
end

function Multivector3(r::Multivector1, s::Multivector1, t::Multivector1)
    Multivector3(Multivector2(r, s), Multivector2(t))
end

function Multivector3(r::Multivector1, s::Multivector1)
    Multivector3(Multivector2(r, s), Multivector2(zero(s)))
end

function Multivector3(r::Multivector1)
    Multivector3(Multivector2(r))
end


function show(io::IO, z::Multivector3)
    print(io, "[1: ")
    print(io, z.l.l.l)
    print(io, ", W: ")
    print(io, z.l.l.r)
    print(io, ", X: ")
    print(io, z.l.r.l)
    print(io, ", WX: ")
    print(io, z.l.r.r)
    print(io, ", Y: ")
    print(io, z.r.l.l)
    print(io, ", WY: ")
    print(io, z.r.l.r)
    print(io, ", XY: ")
    print(io, z.r.r.l)
    print(io, ", (WX)Y: ")
    print(io, z.r.r.r)
    print(io, "]")
end

function real(z::Multivector3)
    real(z.l)
end

function unreal(z::Multivector3)
    vcat(unreal(z.l), asarray(z.r))
end

function zero(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(zero(z.l), zero(z.r))
end

function zero(::Type{Multivector3{T}}) where {T <: Real}
    Multivector3{T}(zero(T))
end

function one(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(one(z.l))
end

function one(::Type{Multivector3{T}}) where {T <: Real}
    Multivector3{T}(one(T))
end

"""
    conj{T <: Real}(z::Multivector3{T})

The `Multivector3` conjugate. If ``z=a+bW+cX+dWX+fY+gWY+hXY+j(WX)Y``,
then `conj(z)` gives
```math
    a-bW-cX-dWX-fY-gWY-hXY-j(WX)Y
```
This operation is an involution.
"""
function conj(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(conj(z.l), -z.r)
end

"""
    cloak{T <: Real}(z::Multivector3{T})

The cloak conjugate changes the sign of even blades. If ``z=a+bW+cX+dWX+fY+gWY+hXY+j(WX)Y``,
then `cloak(z)` gives
```math
    -a+bW+cX-dWX+fY-gWY-hXY+j(WX)Y
```
This operation is equivalent to `-dagger(z)` and thus is also an involution.
"""
function cloak(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(cloak(z.l), dagger(z.r))
end

"""
    dagger{T <: Real}(z::Multivector3{T})

The dagger conjugate changes the sign of odd blades. If ``z=a+bW+cX+dWX+fY+gWY+hXY+j(WX)Y``, then `dagger(z)` gives
```math
    a-bW-cX+dWX-fY+gWY+hXY-j(WX)Y
```
This operation is an involution.
"""
function dagger(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(dagger(z.l), cloak(z.r))
end

"""
    star{T <: Real}(z::Multivector3{T})

Returns the Hodge star conjugate. If ``z=a+bW+cX+dWX+fY+gWY+hXY+j(WX)Y``, then `star(z)` gives
```math
    j+hW-gX+fWX+dY-cWY+bXY+a(WX)Y
```
This operation is an involution.
"""
function star(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(dagger(star(z.r)), star(z.l))
end

"""
    selfstar{T <: Real}(z::Multivector3{T})

The self-star-conjugate part. If ``z=a+bW+cX+dWX+fY+gWY+hXY+j(WX)Y``, then `selfstar(z)` gives
```math
    \\frac{1}{2}(a+j)(1+(WX)Y) + \\frac{1}{2}(b+h)(W+XY) + \\frac{1}{2}(c-g)(X-WY) + \\frac{1}{2}(d+f)(WX+Y)
```
This operation is idempotent.
"""
function selfstar(z::Multivector3{T}) where {T <: Real}
    (z + star(z)) / 2
end

"""
    antiselfstar{T <: Real}(z::Multivector3{T})

The anti-self-star-conjugate part. If ``z=a+bW+cX+dWX``, then `antiselfstar(z)` gives
```math
    \\frac{1}{2}(a-j)(1-(WX)Y) + \\frac{1}{2}(b-h)(W-XY) + \\frac{1}{2}(c+g)(X+WY) + \\frac{1}{2}(d-f)(WX-Y)
```
This operation is idempotent.
"""
function antiselfstar(z::Multivector3{T}) where {T <: Real}
    (z - star(z)) / 2
end

(+)(z::Multivector3{T}) where {T <: Real} = z

function (+)(x::Multivector3, y::Multivector3)
    Multivector3(x.l + y.l, x.r + y.r)
end

function (+)(x::Multivector3, y::Multivector2)
    Multivector3(x.l + y, x.r)
end

function (+)(x::Multivector2, y::Multivector3)
    Multivector3(x + y.l, y.r)
end

function (+)(x::Multivector3, y::Multivector1)
    Multivector3(x.l + y, x.r)
end

function (+)(x::Multivector1, y::Multivector3)
    Multivector3(x + y.l, y.r)
end

function (+)(z::Multivector3, a::Real)
    Multivector3(z.l + a, z.r)
end

function (+)(a::Real, z::Multivector3)
    Multivector3(z.l + a, z.r)
end

function (-)(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(-z.l, -z.r)
end

function (-)(x::Multivector3, y::Multivector3)
    Multivector3(x.l - y.l, x.r - y.r)
end

function (-)(x::Multivector3, y::Multivector2)
    Multivector3(x.l - y, x.r)
end

function (-)(x::Multivector2, y::Multivector3)
    Multivector3(x - y.l, -y.r)
end

function (-)(x::Multivector3, y::Multivector1)
    Multivector3(x.l - y, x.r)
end

function (-)(x::Multivector1, y::Multivector3)
    Multivector3(x - y.l, -y.r)
end

function (-)(z::Multivector3, a::Real)
    Multivector3(z.l - a, z.r)
end

function (-)(a::Real, z::Multivector3)
    Multivector3(a - z.l, -z.r)
end

"""
    (∧)(x::Multivector3, y::Multivector3)

Wedge product of two 2-dimensional multivectors.
"""

function (∧)(x::Multivector3, y::Multivector3)
    Multivector3(x.l ∧ y.l,  (y.r ∧ x.l) + (x.r ∧ conj(y.l)))
end

function (∧)(x::Multivector3, y::Multivector2)
    Multivector3(x.l ∧ y, x.r ∧ conj(y))
end

function (∧)(x::Multivector2, y::Multivector3)
    Multivector3(x ∧ y.l,  y.r ∧ x)
end

function (∧)(x::Multivector3, y::Multivector1)
    Multivector3(x.l ∧ y, x.r ∧ conj(y))
end

function (∧)(x::Multivector1, y::Multivector3)
    Multivector3(x ∧ y.l,  y.r ∧ x)
end

"""
    (*)(z::Multivector1, a::Real)
    (*)(a::Real, z::Multivector1)

Scaling and/or reflection of a `Multivector1` by a real number.
"""
function (*)(z::Multivector3, a::Real)
    Multivector3(z.l * a, z.r * a)
end

function (*)(a::Real, z::Multivector3)
    Multivector3(z.l * a, z.r * a)
end

function (∧)(z::AbstractMultivector3, a::Real)
    z * a
end

function (∧)(a::Real, z::AbstractMultivector3)
    a * z
end

function associator(x::AbstractMultivector, y::AbstractMultivector, z::AbstractMultivector)
    ((x ∧ y) ∧ z) - (x ∧ (y ∧ z))
end

"""
    abs(z::Multivector3)
"""
abs(z::Multivector3) = abs(z.l)

"""
    abs2(z::Multivector3)
"""
abs2(z::Multivector3) = abs2(z.l)

"""
    iszerodivisor(z::Multivector3)

Returns true if `z` is of the form ``bW+cX+dWX`` such that ``z ∧ conj(z) = 0``.
"""
iszerodivisor(z::Multivector3{T}) where {T <: Real} = iszerodivisor(z.l)

"""
    inv(z::Multivector3)

Inverse of a `Multivector3`. An error occurs if argument is a zero divisor.
"""
function inv(z::Multivector3)
    if iszerodivisor(z)
        error(ZeroDivisorInverse)
    end

    conj(z) / abs2(z)
end

function (/)(x::Multivector3, y::Multivector3)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    x * inv(y)
end

function (/)(a::Real, z::Multivector3)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    a * inv(z)
end

function (/)(z::Multivector3, a::Real)
    if iszero(a)
        error(ZeroDenominator)
    end

    Multivector3(z.l / a, z.r / a)
end

function (\)(y::Multivector3, x::Multivector3)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    inv(y) * x
end

function (\)(z::Multivector3, a::Real)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    a * inv(z)
end

function (\)(a::Real, z::Multivector3)
    if iszero(a)
        error(ZeroDenominator)
    end
    
    Multivector3(z.l / a, z.r / a)
end

"""
    crossratioL(w::AbstractMultivector3,
               x::AbstractMultivector3,
               y::AbstractMultivector3,
               z::AbstractMultivector3)

The left cross-ratio:
```julia
    ((x-y) \\ (w-y)) ∧ ((w-z) \\ (x-z))
```
The left cross-ratio is invariant under right Möbius transformations.
"""
function crossratioL(w::AbstractMultivector3,
                     x::AbstractMultivector3,
                     y::AbstractMultivector3,
                     z::AbstractMultivector3)
    (inv(x - y) ∧ (w - y)) ∧ (inv(w - z) ∧ (x - z))
end

"""
    crossratioR(w::AbstractMultivector3,
                x::AbstractMultivector3,
                y::AbstractMultivector3,
                z::AbstractMultivector3)

The right cross-ratio:
```julia
    ((w-y) / (x-y)) ∧ ((x-z) / (w-z))
```
The right cross-ratio is invariant under left Möbius transformations.
"""
function crossratioR(w::AbstractMultivector3,
                     x::AbstractMultivector3,
                     y::AbstractMultivector3,
                     z::AbstractMultivector3)
    ((w - y) ∧ inv(x - y)) ∧ ((x - z) ∧ inv(w - z))
end

"""
    möbiusL(z::AbstractMultivector3,
           a::AbstractMultivector3,
           b::AbstractMultivector3,
           c::AbstractMultivector3,
           d::AbstractMultivector3)
    möbiusL(z::AbstractMultivector3, a::Real, b::Real, c::Real, d::Real)

The left Möbius transformation:
```julia
    ((z ∧ c) + d) \\ ((z ∧ a) + b)
```
This transformation is also know as a fractional linear transformation.
"""
function möbiusL(z::AbstractMultivector3,
                a::AbstractMultivector3,
                b::AbstractMultivector3,
                c::AbstractMultivector3,
                d::AbstractMultivector3)
    inv((z ∧ c) + d) ∧ ((z ∧ a) + b)
end

function möbiusL(z::AbstractMultivector3, a::Real, b::Real, c::Real, d::Real)
    inv((c * z) + d) ∧ ((a * z) + b)
end

"""
    möbiusR(z::AbstractMultivector3,
           a::AbstractMultivector3,
           b::AbstractMultivector3,
           c::AbstractMultivector3,
           d::AbstractMultivector3)
    möbiusR(z::AbstractMultivector3, a::Real, b::Real, c::Real, d::Real)

The right Möbius transformation:
```julia
    ((a ∧ z) + b) / ((c ∧ z) + d)
```
This transformation is also know as a fractional linear transformation.
"""
function möbiusR(z::AbstractMultivector3,
                a::AbstractMultivector3,
                b::AbstractMultivector3,
                c::AbstractMultivector3,
                d::AbstractMultivector3)
    ((a ∧ z) + b) ∧ inv((c ∧ z) + d)
end

function möbiusR(z::AbstractMultivector3, a::Real, b::Real, c::Real, d::Real)
    ((a * z) + b) ∧ inv((c * z) + d)
end
