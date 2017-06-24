"""
    AbstractMultivector1{T <: Real} <: AbstractMultivector{T}

An abstract 1-dimensional multivector.
"""
abstract AbstractMultivector1{T <: Real} <: AbstractMultivector{T}

"""
    Multivector1{T <: Real} <: AbstractMultivector1{T}

An immutable pair of real numbers that represents a member of a 1-dimensional multivector algebra.

Each `Multivector1` has the form
```math
    a+bW
```
where ``a`` and ``b`` are real (and of the same type), and ``W ∧ W = 0``.
Here ``∧`` is the wedge product. Note that
```math
    (a + bW) ∧ (a + bW) = a^{2} + 2abW
```
which, in general, is not equal to zero.
"""
immutable Multivector1{T <: Real} <: AbstractMultivector1{T}
    l::T
    r::T

    Multivector1{U <: Real}(l::U, r::U) = new(l, r)
end

Multivector1{T <: Real}(a::T, b::T) = Multivector1{T}(a, b)
Multivector1(a::Real, b::Real) = Multivector1(promote(a, b)...)
Multivector1{T <: Real}(a::T) = Multivector1{T}(a, zero(T))

function show(io::IO, z::Multivector1)
    print(io, "[1: ")
    print(io, z.l)
    print(io, ", W: ")
    print(io, z.r)
    print(io, "]")
end

function zero{T <: Real}(z::Multivector1{T})
    Multivector1{T}(zero(z.l), zero(z.r))
end

function zero{T <: Real}(::Type{Multivector1{T}})
    Multivector1{T}(zero(T), zero(T))
end

function one{T <: Real}(z::Multivector1{T})
    Multivector1{T}(one(z.l), zero(z.r))
end

function one{T <: Real}(::Type{Multivector1{T}})
    Multivector1{T}(one(T), zero(T))
end

"""
    conj{T <: Real}(z::Multivector1{T})

The `Multivector1` conjugate. If ``z=a+bW``, then `conj(z)` gives
```math
    a-bW
```
This operation is an involution.
"""
function conj{T <: Real}(z::Multivector1{T})
    Multivector1{T}(z.l, -z.r)
end

"""
    dagger{T <: Real}(z::Multivector1{T})

The dagger conjugate changes the sign of odd blades. If ``z=a+bW``, then `dagger(z)` gives
```math
    a-bW
```
This operation is equivalent to `conj(z)` and thus is also an involution.
"""
function dagger(z::Multivector1)
    conj(z)
end

"""
    star{T <: Real}(z::Multivector1{T})

Returns the Hodge star conjugate. If ``z=a+bW``, then `star(z)` gives
```math
    b+aW
```
This operation is an involution.
"""
function star{T <: Real}(z::Multivector1{T})
    Multivector1{T}(z.r, z.l)
end

"""
    selfstar{T <: Real}(z::Multivector1{T})

The self-star-conjugate part. If ``z=a+bW``, then `selfstar(z)` gives
```math
    \\frac{1}{2}(a+b)(1+W)
```
The result is a `SelfStar1`. This operation is idempotent.
"""
function selfstar{T <: Real}(z::Multivector1{T})
    SelfStar1{T}((z.l + z.r)/2)
end

"""
    antiselfstar{T <: Real}(z::Multivector1{T})

The anti-self-star-conjugate part. If ``z=a+bW``, then `antiselfstar(z)` gives
```math
    \\frac{1}{2}(a-b)(1-W)
```
The result is an `AntiSelfStar1`. This operation is idempotent.
"""
function antiselfstar{T <: Real}(z::Multivector1{T})
    AntiSelfStar1{T}((z.l - z.r)/2)
end

(+){T <: Real}(z::Multivector1{T}) = z

function (+)(x::Multivector1, y::Multivector1)
    Multivector1(x.l + y.l, x.r + y.r)
end

function (+)(z::Multivector1, a::Real)
    Multivector1(z.l + a, z.r)
end

function (+)(a::Real, z::Multivector1)
    Multivector1(z.l + a, z.r)
end

function (-){T <: Real}(z::Multivector1{T})
    Multivector1{T}(-z.l, -z.r)
end

function (-)(x::Multivector1, y::Multivector1)
    Multivector1(x.l - y.l, x.r - y.r)
end

function (-)(z::Multivector1, a::Real)
    Multivector1(z.l - a, z.r)
end

function (-)(a::Real, z::Multivector1)
    Multivector1(a - z.l, -z.r)
end

"""
    (∧)(x::Multivector1, y::Multivector1)

Wedge product of two 1-dimensional multivectors.
"""

function (∧)(x::Multivector1, y::Multivector1)
    Multivector1(x.l*y.l, x.r*y.l + y.r*x.l)
end

"""
    (*)(z::Multivector1, a::Real)
    (*)(a::Real, z::Multivector1)

Scaling and/or reflection of a `Multivector1` by a real number.
"""
function (*)(z::Multivector1, a::Real)
    Multivector1(z.l*a, z.r*a)
end

function (*)(a::Real, z::Multivector1)
    Multivector1(z.l*a, z.r*a)
end

function (∧)(z::AbstractMultivector1, a::Real)
    z * a
end

function (∧)(a::Real, z::AbstractMultivector1)
    a * z
end

"""
    abs(z::Multivector1)
"""
abs(z::Multivector1) = abs(z.l)

"""
    abs2(z::Multivector1)
"""
abs2(z::Multivector1) = z.l^2

"""
    iszerodivisor(z::Multivector1)

Returns true if `z` is of the form ``bW`` such that ``z ∧ conj(z) = 0``.
Note that it follows that also ``z ∧ z = 0``.
"""
iszerodivisor{T <: Real}(z::Multivector1{T}) = z.l == zero(T)

"""
    inv(z::Multivector1)

Inverse of a `Multivector1`. An error occurs if argument is a zero divisor.
"""
function inv(z::Multivector1)
    if iszerodivisor(z)
        error(ZeroDivisorInverse)
    end

    conj(z) / abs2(z)
end

function (/)(x::Multivector1, y::Multivector1)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    Multivector1(x.l*y.l, x.r*y.l - y.r*x.l) / abs2(y)
end

function (/)(a::Real, z::Multivector1)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    a * inv(z)
end

function (/)(z::Multivector1, a::Real)
    if a == zero(a)
        error(ZeroDenominator)
    end

    Multivector1(z.l / a, z.r / a)
end

function (\)(y::Multivector1, x::Multivector1)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    Multivector1(x.l*y.l, x.r*y.l - y.r*x.l) / abs2(y)
end

function (\)(z::Multivector1, a::Real)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    a * inv(z)
end

function (\)(a::Real, z::Multivector1)
    if a == zero(a)
        error(ZeroDenominator)
    end
    
    Multivector1(z.l / a, z.r / a)
end

"""
    crossratio(w::AbstractMultivector1,
               x::AbstractMultivector1,
               y::AbstractMultivector1,
               z::AbstractMultivector1)

The cross-ratio:
```math
    \\frac{(w-y)(x-z)}{(x-y)(w-z)}
```
The cross-ratio is invariant under Möbius transformations.
"""
function crossratio(w::AbstractMultivector1,
                    x::AbstractMultivector1,
                    y::AbstractMultivector1,
                    z::AbstractMultivector1)
    inv(x-y) ∧ (w-y) ∧ inv(w-z) ∧ (x-z)
end

"""
    möbius(z::AbstractMultivector1,
           a::AbstractMultivector1,
           b::AbstractMultivector1,
           c::AbstractMultivector1,
           d::AbstractMultivector1)
    möbius(z::AbstractMultivector1, a::Real, b::Real, c::Real, d::Real)

The Möbius transformation:
```math
    \\frac{az + b}{cz + d}
```
This transformation is also know as a fractional linear transformation.
"""
function möbius(z::AbstractMultivector1,
                a::AbstractMultivector1,
                b::AbstractMultivector1,
                c::AbstractMultivector1,
                d::AbstractMultivector1)
    ((a ∧ z) + b) ∧ inv((c ∧ z) + d)
end

function möbius(z::AbstractMultivector1, a::Real, b::Real, c::Real, d::Real)
    ((a * z) + b) ∧ inv((c * z) + d)
end
