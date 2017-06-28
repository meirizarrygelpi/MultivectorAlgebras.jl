"""
    Multivector3{T <: Real} <: AbstractMultivector{T}

An immutable 8-tuple of real numbers that represents a member
of a 3-dimensional multivector algebra.

Each `Multivector3` has the form
```math
    p+qC
```
where ``p`` and ``q`` are 2-dimensional multivectors, or
```math
    a+bA+cB+dAB+fC+gAC+hBC+j(AB)C
```
where ``a``, ``b``, ``c``, ``d``, ``f``, ``g``, ``h``, and ``j`` are real
(and of the same type), and ``A * A = 0``, ``B * B = 0``, `` C * C = 0``;
``AB = A * B = -B * A``, ``AC = A * C = -C * A``, ``BC = B * C = -C * B``;
and ``(AB)C = (A * B) * C``.
Here ``*`` is the wedge product.
"""
struct Multivector3{T <: Real} <: AbstractMultivector{T}
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
    print(io, ", A: ")
    print(io, z.l.l.r)
    print(io, ", B: ")
    print(io, z.l.r.l)
    print(io, ", AB: ")
    print(io, z.l.r.r)
    print(io, ", C: ")
    print(io, z.r.l.l)
    print(io, ", AC: ")
    print(io, z.r.l.r)
    print(io, ", BC: ")
    print(io, z.r.r.l)
    print(io, ", (AB)C: ")
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

The `Multivector3` conjugate. If ``z=a+bA+cB+dAB+fC+gAC+hBC+j(AB)C``,
then `conj(z)` gives
```math
    a-bA-cB-dAB-fC-gAC-hBC-j(AB)C
```
This operation is an involution.
"""
function conj(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(conj(z.l), -z.r)
end

"""
    cloak{T <: Real}(z::Multivector3{T})

The cloak conjugate changes the sign of even blades. If ``z=a+bA+cB+dAB+fC+gAC+hBC+j(AB)C``,
then `cloak(z)` gives
```math
    -a+bA+cB-dAB+fC-gAC-hBC+j(AB)C
```
This operation is equivalent to `-dagger(z)` and thus is also an involution.
"""
function cloak(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(cloak(z.l), dagger(z.r))
end

"""
    dagger{T <: Real}(z::Multivector3{T})

The dagger conjugate changes the sign of odd blades. If ``z=a+bA+cB+dAB+fC+gAC+hBC+j(AB)C``, then `dagger(z)` gives
```math
    a-bA-cB+dAB-fC+gAC+hBC-j(AB)C
```
This operation is an involution.
"""
function dagger(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(dagger(z.l), cloak(z.r))
end

"""
    star{T <: Real}(z::Multivector3{T})

Returns the Hodge star conjugate. If ``z=a+bA+cB+dAB+fC+gAC+hBC+j(AB)C``, then `star(z)` gives
```math
    j+hA-gB+fAB+dC-cAC+bBC+a(AB)C
```
This operation is an involution.
"""
function star(z::Multivector3{T}) where {T <: Real}
    Multivector3{T}(dagger(star(z.r)), star(z.l))
end

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
    (*)(x::Multivector3, y::Multivector3)

Wedge product of two 3-dimensional multivectors.
This operation is non-commutative and non-associative.
"""
function (*)(x::Multivector3, y::Multivector3)
    Multivector3(x.l * y.l,  (y.r * x.l) + (x.r * conj(y.l)))
end

function (*)(x::Multivector3, y::Multivector2)
    Multivector3(x.l * y, x.r * conj(y))
end

function (*)(x::Multivector2, y::Multivector3)
    Multivector3(x * y.l,  y.r * x)
end

function (*)(x::Multivector3, y::Multivector1)
    Multivector3(x.l * y, x.r * conj(y))
end

function (*)(x::Multivector1, y::Multivector3)
    Multivector3(x * y.l,  y.r * x)
end

"""
    (*)(z::Multivector3, a::Real)
    (*)(a::Real, z::Multivector3)

Scaling and/or reflection of a `Multivector3` by a real number.
"""
function (*)(z::Multivector3, a::Real)
    Multivector3(z.l * a, z.r * a)
end

function (*)(a::Real, z::Multivector3)
    Multivector3(a * z.l, a * z.r)
end

function (/)(x::Multivector3, y::Multivector3)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    x * inv(y)
end

function (\)(y::Multivector3, x::Multivector3)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    inv(y) * x
end

function (/)(z::Multivector3, a::Real)
    Multivector3(z.l / a, z.r / a)
end

function (\)(a::Real, z::Multivector3)
    Multivector3(a \ z.l, a \ z.r)
end
