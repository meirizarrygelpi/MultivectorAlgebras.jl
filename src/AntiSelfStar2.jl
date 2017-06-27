"""
    AntiSelfStar2{T <: Real} <: AbstractMultivector2{T}

An anti-self-star-conjugate 1-dimensional multivector.
"""
struct AntiSelfStar2{T <: Real} <: AbstractMultivector2{T}
    c::T

    AntiSelfStar2{U}(c::U) where {U <: Real} = new(c)
end

AntiSelfStar2(c::T) where {T <: Real} = AntiSelfStar2{T}(c)

function Multivector2(z::AntiSelfStar2{T}) where {T <: Real}
    Multivector2{T}(z.c, zero(T), zero(T), -z.c)
end

function show(io::IO, z::AntiSelfStar2)
    print(io, "[(1-WX): ")
    print(io, z.c)
    print(io, "]")
end

function real(z::AntiSelfStar2)
    z.c
end

function unreal(z::AntiSelfStar2)
    vcat(zero(z.c), zero(z.c), -z.c)
end

function zero(z::AntiSelfStar2{T}) where {T <: Real}
    AntiSelfStar2(zero(z.c))
end

function zero(::Type{AntiSelfStar2{T}}) where {T <: Real}
    AntiSelfStar2(zero(T))
end

function one(z::AntiSelfStar2{T}) where {T <: Real}
    Multivector2(one(z.c), zero(z.c), zero(z.c), zero(z.c))
end

function one(::Type{AntiSelfStar2{T}}) where {T <: Real}
    Multivector2(one(T), zero(T), zero(T), zero(T))
end

"""
    conj{T <: Real}(z::AntiSelfStar2{T})

Turns an anti-self-star-conjugate 2-dimensional multivector into an
self-star-conjugate 2-dimensional multivector.
This operation is an involution.
"""
conj(z::AntiSelfStar2{T}) where {T <: Real} = SelfStar2{T}(z.c)

cloak(z::AntiSelfStar2{T}) where {T <: Real} = AntiSelfStar2{T}(-z.c)

dagger(z::AntiSelfStar2) = z

"""
    star(z::AntiSelfStar2)

The negation operation. This operation is an involution.
"""
star(z::AntiSelfStar2) = -z

(+)(z::AntiSelfStar2) = z

function (+)(x::AntiSelfStar2, y::AntiSelfStar2)
    AntiSelfStar2(x.c + y.c)
end

function (+)(x::Multivector2, y::AntiSelfStar2)
    Multivector2(x.l.l + y.c, zero(y.c), zero(y.c), x.r.r - y.c)
end

function (+)(x::AntiSelfStar2, y::Multivector2)
    Multivector2(x.c + y.l.l, zero(x.c), zero(x.c), y.r.r - x.c)
end

function (+)(z::AntiSelfStar2, a::Real)
    Multivector2(z.c + a, zero(a), zero(a), -z.c)
end

function (+)(a::Real, z::AntiSelfStar2)
    Multivector2(z.c + a, zero(a), zero(a), -z.c)
end

function (-)(z::AntiSelfStar2{T}) where {T <: Real}
    AntiSelfStar2{T}(-z.c)
end

function (-)(x::AntiSelfStar2, y::AntiSelfStar2)
    AntiSelfStar2(x.c - y.c)
end

function (-)(x::Multivector2, y::AntiSelfStar2)
    Multivector2(x.l.l - y.c, zero(y.c), zero(y.c), x.r.r + y.c)
end

function (-)(x::AntiSelfStar2, y::Multivector2)
    Multivector2(x.c - y.l.l, zero(y.c), zero(y.c), -(x.c + y.r.r))
end

function (-)(z::AntiSelfStar2, a::Real)
    Multivector2(z.c - a, zero(a), zero(a), -z.c)
end

function (-)(a::Real, z::AntiSelfStar2)
    Multivector2(a - z.c, zero(a), zero(a), z.c)
end

function (∧)(x::AntiSelfStar2, y::AntiSelfStar2)
    xy = x.c * y.c
    Multivector2(xy, zero(xy), zero(xy), -2xy)
end

function (∧)(x::AntiSelfStar2, y::Multivector2)
    Multivector2(x.c * y.l.l, zero(x.c), zero(x.c), x.c * (y.r.r - y.l.l))
end

function (∧)(x::Multivector2, y::AntiSelfStar2)
    Multivector2(x.l.l * y.c, zero(y.c), zero(y.c), (x.r.r - x.l.l) * y.c)
end

function (∧)(x::SelfStar2, y::AntiSelfStar2)
    x.c * y.c
end

function (∧)(x::AntiSelfStar2, y::SelfStar2)
    x.c * y.c
end

function (*)(z::AntiSelfStar2, a::Real)
    AntiSelfStar2(a * z.c)
end

function (*)(a::Real, z::AntiSelfStar2)
    AntiSelfStar2(a * z.c)
end

"""
    abs(z::AntiSelfStar2)
"""
abs(z::AntiSelfStar2) = abs(z.c)

"""
    abs2(z::AntiSelfStar2)
"""
abs2(z::AntiSelfStar2) = z.c^2

function inv(z::AntiSelfStar2)
    if iszero(z.c)
        error(ZeroInverse)
    end

    SelfStar2(1 / z.c)
end

function (/)(x::AntiSelfStar2, y::AntiSelfStar2)
    if iszero(y.c)
        error(ZeroDenominator)
    end

    x.c / y.c
end

function (/)(a::Real, z::AntiSelfStar2)
    if iszero(z.c)
        error(ZeroDenominator)
    end

    SelfStar2(a / z.c)
end

function (/)(z::AntiSelfStar2, a::Real)
    if iszero(a)
        error(ZeroDenominator)
    end

    AntiSelfStar2(z.c / a)
end

function (\)(y::AntiSelfStar2, x::AntiSelfStar2)
    if iszero(y.c)
        error(ZeroDenominator)
    end

    y.c \ x.c
end

function (\)(z::AntiSelfStar2, a::Real)
    if iszero(z.c)
        error(ZeroDenominator)
    end

    SelfStar2(z.c \ a)
end

function (\)(a::Real, z::AntiSelfStar2)
    if iszero(a)
        error(ZeroDenominator)
    end
    
    AntiSelfStar2(a \ z.c)
end
