"""
    Multivector4{T <: Real} <: AbstractMultivector{T}

An immutable 16-tuple of real numbers that represents a member
of a 4-dimensional multivector algebra.

Each `Multivector4` has the form
```math
    p+qD
```
where ``p`` and ``q`` are 3-dimensional multivectors.
"""
struct Multivector4{T <: Real} <: AbstractMultivector{T}
    l::Multivector3{T}
    r::Multivector3{T}
    
    Multivector4{U}(l::Multivector3{U}, r::Multivector3{U}) where {U <: Real} = new(l, r)
end

function Multivector4(x::Multivector3{T}, y::Multivector3{T}) where {T <: Real}
    Multivector4{T}(x, y)
end

function Multivector4(x::Multivector3)
    Multivector4(x, zero(x))
end

function Multivector4(x::Multivector2)
    Multivector4(Multivector3(x), Multivector3(zero(x)))
end

function Multivector4(x::Multivector1)
    Multivector4(Multivector3(x), Multivector3(zero(x)))
end

function show(io::IO, z::Multivector4)
    print(io, "[1: ")
    print(io, z.l.l.l.l)
    print(io, ", A: ")
    print(io, z.l.l.l.r)
    print(io, ", B: ")
    print(io, z.l.l.r.l)
    print(io, ", AB: ")
    print(io, z.l.l.r.r)
    print(io, ", C: ")
    print(io, z.l.r.l.l)
    print(io, ", AC: ")
    print(io, z.l.r.l.r)
    print(io, ", BC: ")
    print(io, z.l.r.r.l)
    print(io, ", (AB)C: ")
    print(io, z.l.r.r.r)
    print(io, ",\n")
    print(io, " D: ")
    print(io, z.r.l.l.l)
    print(io, ", AD: ")
    print(io, z.r.l.l.r)
    print(io, ", BD: ")
    print(io, z.r.l.r.l)
    print(io, ", (AB)D: ")
    print(io, z.r.l.r.r)
    print(io, ", CD: ")
    print(io, z.r.r.l.l)
    print(io, ", (AC)D: ")
    print(io, z.r.r.l.r)
    print(io, ", (BC)D: ")
    print(io, z.r.r.r.l)
    print(io, ", ((AB)C)D: ")
    print(io, z.r.r.r.r)
    print(io, "]")
end

function real(z::Multivector4)
    real(z.l)
end

function unreal(z::Multivector4)
    vcat(unreal(z.l), asarray(z.r))
end

function zero(z::Multivector4{T}) where {T <: Real}
    Multivector4{T}(zero(z.l), zero(z.r))
end

function zero(::Type{Multivector4{T}}) where {T <: Real}
    Multivector4{T}(zero(T))
end

function one(z::Multivector4{T}) where {T <: Real}
    Multivector4{T}(one(z.l))
end

function one(::Type{Multivector4{T}}) where {T <: Real}
    Multivector4{T}(one(T))
end

"""
    conj{T <: Real}(z::Multivector4{T})

The `Multivector4` conjugate.
This operation is an involution.
"""
function conj(z::Multivector4{T}) where {T <: Real}
    Multivector4{T}(conj(z.l), -z.r)
end

"""
    cloak{T <: Real}(z::Multivector4{T})

The cloak conjugate changes the sign of even blades.
This operation is equivalent to `-dagger(z)` and thus is also an involution.
"""
function cloak(z::Multivector4{T}) where {T <: Real}
    Multivector4{T}(cloak(z.l), dagger(z.r))
end

"""
    dagger{T <: Real}(z::Multivector4{T})

The dagger conjugate changes the sign of odd blades.
This operation is an involution.
"""
function dagger(z::Multivector4{T}) where {T <: Real}
    Multivector4{T}(dagger(z.l), cloak(z.r))
end

"""
    star{T <: Real}(z::Multivector4{T})

Returns the Hodge star conjugate.
This operation is not an involution, but `star(star(z)) = dagger(z)`.
"""
function star(z::Multivector4{T}) where {T <: Real}
    Multivector4{T}(dagger(star(z.r)), star(z.l))
end

function (+)(x::Multivector4, y::Multivector4)
    Multivector4(x.l + y.l, x.r + y.r)
end

function (+)(x::Multivector4, y::Multivector3)
    Multivector4(x.l + y, x.r)
end

function (+)(x::Multivector3, y::Multivector4)
    Multivector4(x + y.l, y.r)
end

function (+)(x::Multivector4, y::Multivector2)
    Multivector4(x.l + y, x.r)
end

function (+)(x::Multivector2, y::Multivector4)
    Multivector4(x + y.l, y.r)
end

function (+)(x::Multivector4, y::Multivector1)
    Multivector4(x.l + y, x.r)
end

function (+)(x::Multivector1, y::Multivector4)
    Multivector4(x + y.l, y.r)
end

function (+)(z::Multivector4, a::Real)
    Multivector4(z.l + a, z.r)
end

function (+)(a::Real, z::Multivector4)
    Multivector4(z.l + a, z.r)
end

function (-)(z::Multivector4{T}) where {T <: Real}
    Multivector4{T}(-z.l, -z.r)
end

function (-)(x::Multivector4, y::Multivector4)
    Multivector4(x.l - y.l, x.r - y.r)
end

function (-)(x::Multivector4, y::Multivector3)
    Multivector4(x.l - y, x.r)
end

function (-)(x::Multivector3, y::Multivector4)
    Multivector4(x - y.l, -y.r)
end

function (-)(x::Multivector4, y::Multivector2)
    Multivector4(x.l - y, x.r)
end

function (-)(x::Multivector2, y::Multivector4)
    Multivector4(x - y.l, -y.r)
end

function (-)(x::Multivector4, y::Multivector1)
    Multivector4(x.l - y, x.r)
end

function (-)(x::Multivector1, y::Multivector4)
    Multivector4(x - y.l, -y.r)
end

function (-)(z::Multivector4, a::Real)
    Multivector4(z.l - a, z.r)
end

function (-)(a::Real, z::Multivector4)
    Multivector4(a - z.l, -z.r)
end

"""
    (*)(x::Multivector4, y::Multivector4)

Wedge product of two 3-dimensional multivectors.
This operation is non-commutative and non-associative.
"""
function (*)(x::Multivector4, y::Multivector4)
    Multivector4(x.l * y.l,  (y.r * x.l) + (x.r * conj(y.l)))
end

function (*)(x::Multivector4, y::Multivector3)
    Multivector4(x.l * y, x.r * conj(y))
end

function (*)(x::Multivector3, y::Multivector4)
    Multivector4(x * y.l,  y.r * x)
end

function (*)(x::Multivector4, y::Multivector2)
    Multivector4(x.l * y, x.r * conj(y))
end

function (*)(x::Multivector2, y::Multivector4)
    Multivector4(x * y.l,  y.r * x)
end

function (*)(x::Multivector4, y::Multivector1)
    Multivector4(x.l * y, x.r * conj(y))
end

function (*)(x::Multivector1, y::Multivector4)
    Multivector4(x * y.l,  y.r * x)
end

"""
    (*)(z::Multivector4, a::Real)
    (*)(a::Real, z::Multivector4)

Scaling and/or reflection of a `Multivector4` by a real number.
"""
function (*)(z::Multivector4, a::Real)
    Multivector4(z.l * a, z.r * a)
end

function (*)(a::Real, z::Multivector4)
    Multivector4(z.l * a, z.r * a)
end

function (/)(x::Multivector4, y::Multivector4)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    x * inv(y)
end

function (/)(z::Multivector4, a::Real)
    if iszero(a)
        error(ZeroDenominator)
    end

    Multivector4(z.l / a, z.r / a)
end

function (\)(y::Multivector4, x::Multivector4)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    inv(y) * x
end

function (\)(a::Real, z::Multivector4)
    if iszero(a)
        error(ZeroDenominator)
    end
    
    Multivector4(z.l / a, z.r / a)
end
