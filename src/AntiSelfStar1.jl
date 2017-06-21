"""
    AntiSelfStar1{T <: Real} <: AbstractMultivector1{T}

An anti-self-star-conjugate 1-dimensional multivector.
"""
immutable AntiSelfStar1{T <: Real} <: AbstractMultivector1{T}
    c::T

    AntiSelfStar1{U <: Real}(c::U) = new(c)
end

AntiSelfStar1{T <: Real}(c::T) = AntiSelfStar1{T}(c)

function Multivector1{T<:Real}(z::AntiSelfStar1{T})
    Multivector1{T}(z.c, -z.c)
end

function show(io::IO, z::AntiSelfStar1)
    print(io, "[(1-W): ")
    print(io, z.c)
    print(io, "]")
end

function zero{T <: Real}(z::AntiSelfStar1{T})
    AntiSelfStar1(zero(z.c))
end

function zero{T <: Real}(::Type{AntiSelfStar1{T}})
    AntiSelfStar1(zero(T))
end

function one{T <: Real}(z::AntiSelfStar1{T})
    Multivector1(one(z.c), zero(z.c))
end

function one{T <: Real}(::Type{AntiSelfStar1{T}})
    Multivector1(one(T), zero(T))
end

"""
    conj(z::AntiSelfStar1)

Turns an anti-self-star-conjugate 1-dimensional multivector into an
self-star-conjugate 1-dimensional multivector.
This operation is an involution.
"""
conj{T <: Real}(z::AntiSelfStar1{T}) = SelfStar1{T}(z.c)

"""
    dagger(z::AntiSelfStar1)

Turns a self-star-conjugate 1-dimensional multivector into an
anti-self-star-conjugate 1-dimensional multivector.
This operation is an involution.
"""
dagger(z::AntiSelfStar1) = conj(z)

"""
    star(z::AntiSelfStar1)

The negation operation. This operation is an involution.
"""
star(z::AntiSelfStar1) = -z

(+)(z::AntiSelfStar1) = z

function (+){T <: Real}(x::AntiSelfStar1{T}, y::AntiSelfStar1{T})
    AntiSelfStar1{T}(x.c + y.c)
end

function (+){T <: Real}(x::Multivector1{T}, y::AntiSelfStar1{T})
    Multivector1{T}(x.l + y.c, x.r - y.c)
end

function (+){T <: Real}(x::AntiSelfStar1{T}, y::Multivector1{T})
    Multivector1{T}(x.c + y.l, y.r - x.c)
end

function (+){T <: Real}(z::AntiSelfStar1{T}, a::T)
    Multivector1{T}(z.c + a, -z.c)
end

function (+){T <: Real}(a::T, z::AntiSelfStar1{T})
    Multivector1{T}(z.c + a, -z.c)
end

function (+)(x::AntiSelfStar1, y::AntiSelfStar1)
    AntiSelfStar1(x.c + y.c)
end

function (+)(x::Multivector1, y::AntiSelfStar1)
    Multivector1(x.l + y.c, x.r - y.c)
end

function (+)(x::AntiSelfStar1, y::Multivector1)
    Multivector1(x.c + y.l, y.r - x.c)
end

function (+)(z::AntiSelfStar1, a::Real)
    Multivector1(z.c + a, -z.c)
end

function (+)(a::Real, z::AntiSelfStar1)
    Multivector1(z.c + a, -z.c)
end

function (-){T <: Real}(z::AntiSelfStar1{T})
    AntiSelfStar1{T}(-z.c)
end

function (-){T <: Real}(x::AntiSelfStar1{T}, y::AntiSelfStar1{T})
    AntiSelfStar1{T}(x.c - y.c)
end

function (-){T <: Real}(x::Multivector1{T}, y::AntiSelfStar1{T})
    Multivector1{T}(x.l - y.c, x.r + y.c)
end

function (-){T <: Real}(x::AntiSelfStar1{T}, y::Multivector1{T})
    Multivector1{T}(x.c - y.l, -(x.c + y.r))
end

function (-){T <: Real}(z::AntiSelfStar1{T}, a::T)
    Multivector1{T}(z.c - a, -z.c)
end

function (-){T <: Real}(a::T, z::AntiSelfStar1{T})
    Multivector1{T}(a - z.c, z.c)
end

function (-)(x::AntiSelfStar1, y::AntiSelfStar1)
    AntiSelfStar1(x.c - y.c)
end

function (-)(x::Multivector1, y::AntiSelfStar1)
    Multivector1(x.l - y.c, x.r + y.c)
end

function (-)(x::AntiSelfStar1, y::Multivector1)
    Multivector1(x.c - y.l, -(x.c + y.r))
end

function (-)(z::AntiSelfStar1, a::Real)
    Multivector1(z.c - a, -z.c)
end

function (-)(a::Real, z::AntiSelfStar1)
    Multivector1(a - z.c, z.c)
end

function (∧){T <: Real}(x::AntiSelfStar1{T}, y::AntiSelfStar1{T})
    xy = x.c * y.c
    Multivector1{T}(xy, -2xy)
end

function (∧)(x::AntiSelfStar1, y::AntiSelfStar1)
    xy = x.c * y.c
    Multivector1(xy, -2xy)
end

function (∧){T <: Real}(x::AntiSelfStar1{T}, y::Multivector1{T})
    Multivector1{T}(x.c*y.l, x.c*(y.r - y.l))
end

function (∧){T <: Real}(x::Multivector1{T}, y::AntiSelfStar1{T})
    Multivector1{T}(y.c*x.l, y.c*(x.r - x.l))
end

function (∧)(x::AntiSelfStar1, y::Multivector1)
    Multivector1(x.c*y.l, x.c*(y.r - y.l))
end

function (∧)(x::Multivector1, y::AntiSelfStar1)
    Multivector1(y.c*x.l, y.c*(x.r - x.l))
end

function (∧)(x::SelfStar1, y::AntiSelfStar1)
    x.c * y.c
end

function (∧)(x::AntiSelfStar1, y::SelfStar1)
    x.c * y.c
end

function (*){T <: Real}(z::AntiSelfStar1{T}, a::T)
    AntiSelfStar1{T}(a * z.c)
end

function (*){T <: Real}(a::T, z::AntiSelfStar1{T})
    AntiSelfStar1{T}(a * z.c)
end

function (*)(z::AntiSelfStar1, a::Real)
    AntiSelfStar1(a * z.c)
end

function (*)(a::Real, z::AntiSelfStar1)
    AntiSelfStar1(a * z.c)
end

"""
    abs(z::AntiSelfStar1)
"""
abs(z::AntiSelfStar1) = abs(z.c)

"""
    abs2(z::AntiSelfStar1)
"""
abs2(z::AntiSelfStar1) = z.c^2

function inv(z::AntiSelfStar1)
    if z.c == zero(z.c)
        error(ZeroInverse)
    end

    SelfStar1(1 / z.c)
end

function (/)(x::AntiSelfStar1, y::AntiSelfStar1)
    if y.c == zero(y.c)
        error(ZeroDenominator)
    end

    x.c / y.c
end

function (/)(a::Real, z::AntiSelfStar1)
    if z.c == zero(z.c)
        error(ZeroDenominator)
    end

    SelfStar1(a / z.c)
end

function (/)(z::AntiSelfStar1, a::Real)
    if a == zero(a)
        error(ZeroDenominator)
    end

    AntiSelfStar1(z.c / a)
end

function (\)(y::AntiSelfStar1, x::AntiSelfStar1)
    if y.c == zero(y.c)
        error(ZeroDenominator)
    end

    y.c \ x.c
end

function (\)(z::AntiSelfStar1, a::Real)
    if z.c == zero(z.c)
        error(ZeroDenominator)
    end

    SelfStar1(z.c \ a)
end

function (\)(a::Real, z::AntiSelfStar1)
    if a == zero(a)
        error(ZeroDenominator)
    end
    
    AntiSelfStar1(a \ z.c)
end