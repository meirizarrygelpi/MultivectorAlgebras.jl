"""
    Multivector2{T <: Real} <: AbstractMultivector{T}

An immutable 4-tuple of real numbers that represents a member
of a 2-dimensional multivector algebra.

Each `Multivector2` has the form
```math
    p+qB
```
where ``p`` and ``q`` are 1-dimensional multivectors, or
```math
    a+bA+cB+dAB
```
where ``a``, ``b``, ``c``, and ``d`` are real (and of the same type),
and ``A * A = 0``, ``B * B = 0``, and ``AB = A * B = -B * A``.
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
    Multivector2{T}(Multivector1{T}(a, b), Multivector1{T}(c, d))
end

function Multivector2(a::T, b::T, c::T) where {T <: Real}
    Multivector2{T}(Multivector1{T}(a, b), Multivector1{T}(c, zero(T)))
end

function Multivector2(a::T, b::T) where {T <: Real}
    Multivector2{T}(Multivector1{T}(a, b), zero(Multivector1{T}))
end

function Multivector2(a::T) where {T <: Real}
    Multivector2{T}(Multivector1{T}(a, zero(T)), zero(Multivector1{T}))
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
    print(io, ", A: ")
    print(io, z.l.r)
    print(io, ", B: ")
    print(io, z.r.l)
    print(io, ", AB: ")
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
    Multivector2(zero(T), zero(T), zero(T), zero(T))
end

function one(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(one(z.l), zero(z.r))
end

function one(::Type{Multivector2{T}}) where {T <: Real}
    Multivector2(one(T), zero(T), zero(T), zero(T))
end

"""
    conj{T <: Real}(z::Multivector2{T})

The `Multivector2` conjugate. If ``z=a+bA+cB+dAB``, then `conj(z)` gives
```math
    a-bA-cB-dAB
```
This operation is an involution.
"""
function conj(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(conj(z.l), -z.r)
end

"""
    cloak{T <: Real}(z::Multivector2{T})

The cloak conjugate changes the sign of even blades. If ``z=a+bA+cB+dAB``, then `cloak(z)` gives
```math
    -a+bA+cB-dAB
```
This operation is equivalent to `-dagger(z)` and thus is also an involution.
"""
function cloak(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(cloak(z.l), dagger(z.r))
end

"""
    dagger{T <: Real}(z::Multivector2{T})

The dagger conjugate changes the sign of odd blades. If ``z=a+bA+cB+dAB``, then `dagger(z)` gives
```math
    a-bA-cB+dAB
```
This operation is an involution.
"""
function dagger(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(dagger(z.l), cloak(z.r))
end

"""
    star{T <: Real}(z::Multivector2{T})

Returns the Hodge star conjugate. If ``z=a+bA+cB+dAB``, then `star(z)` gives
```math
    d-cA+bB+aAB
```
This operation is not an involution, but `star(star(z)) = dagger(z)`.
"""
function star(z::Multivector2{T}) where {T <: Real}
    Multivector2{T}(dagger(star(z.r)), star(z.l))
end

"""
    selfstar(z::Multivector2)

The self-star-conjugate part.
This operation is idempotent.
"""
function selfstar(z::Multivector2)
    (z + star(z) + dagger(z) + star(dagger(z))) / 4
end

"""
    antiselfstar(z::Multivector2)

The anti-self-star-conjugate part.
This operation is idempotent.
"""
function antiselfstar(z::Multivector2)
    (z - star(z) + dagger(z) - star(dagger(z))) / 4
end

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

Wedge product of 2-dimensional multivectors.
This operation is non-commutative but associative.
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
    Multivector2(a * z.l, a * z.r)
end

function (/)(x::Multivector2, y::Multivector2)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    x * inv(y)
end

function (\)(y::Multivector2, x::Multivector2)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    inv(y) * x
end

function (/)(z::Multivector2, a::Real)
    Multivector2(z.l / a, z.r / a)
end

function (\)(a::Real, z::Multivector2)
    Multivector2(a \ z.l, a \ z.r)
end

function random(::Type{Multivector2{T}}) where T <: Real
    Multivector2{T}(random(Multivector1{T}), random(Multivector1{T}))
end
