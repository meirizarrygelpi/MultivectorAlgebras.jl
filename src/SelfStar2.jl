"""
    SelfStar2{T <: Real} <: AbstractMultivector2{T}

A self-star-conjugate 1-dimensional multivector.
"""
struct SelfStar2{T <: Real} <: AbstractMultivector2{T}
    c::T

    SelfStar2{U}(c::U) where {U <: Real} = new(c)
end

SelfStar2(c::T) where {T <: Real} = SelfStar2{T}(c)

function Multivector2(z::SelfStar2{T}) where {T <: Real}
    Multivector2{T}(z.c, zero(T), zero(T), z.c)
end

function show(io::IO, z::SelfStar2)
    print(io, "[(1+WX): ")
    print(io, z.c)
    print(io, "]")
end

function zero(z::SelfStar2{T}) where {T <: Real}
    SelfStar2{T}(zero(z.c))
end

function zero(::Type{SelfStar2{T}}) where {T <: Real}
    SelfStar2{T}(zero(T))
end

function one(z::SelfStar2{T}) where {T <: Real}
    Multivector2(one(z.c), zero(z.c), zero(z.c), zero(z.c))
end

function one(::Type{SelfStar2{T}}) where {T <: Real}
    Multivector2(one(T), zero(T), zero(T), zero(T))
end

"""
    conj(z::SelfStar2)

Turns a self-star-conjugate 2-dimensional multivector into an
anti-self-star-conjugate 2-dimensional multivector. This operation is an involution.
"""
conj(z::SelfStar2{T}) where {T <: Real} = AntiSelfStar2{T}(z.c)

cloak(z::SelfStar2{T}) where {T <: Real} = SelfStar2{T}(-z.c)

dagger(z::SelfStar2) = z

"""
    star(z::SelfStar2)

The identity operation. This operation is an involution.
"""
star(z::SelfStar2) = z

# Unary +
(+)(z::SelfStar2) = z

# Addition of self-star-conjugate multivectors
function (+)(x::SelfStar2, y::SelfStar2)
    SelfStar2(x.c + y.c)
end

# Addition of a self-star-conjugate multivector and a multivector
function (+)(x::SelfStar2, y::Multivector2)
    Multivector2(x.c + y.l.l, zero(x.c), zero(x.c), x.c + y.r.r)
end

function (+)(x::Multivector2, y::SelfStar2)
    Multivector2(x.c + y.l.l, zero(x.c), zero(x.c), x.c + y.r.r)
end

# Addition of a self-star-conjugate multivector and a real number
function (+)(z::SelfStar2, a::Real)
    Multivector2(z.c + a, zero(a), zero(a), z.c)
end

function (+)(a::Real, z::SelfStar2)
    Multivector2(z.c + a, zero(a), zero(a), z.c)
end

# Unary -
function (-)(z::SelfStar2{T}) where {T <: Real}
    SelfStar2{T}(-z.c)
end

# Subtraction of self-star-conjugate multivectors
function (-)(x::SelfStar2, y::SelfStar2)
    SelfStar2(x.c - y.c)
end

# Subtraction of self-star-conjugate multivector and a multivector
function (-)(x::Multivector2, y::SelfStar2)
    Multivector2(x.l.l - y.c, x.l.r, x.r.l, x.r.r - y.c)
end

function (-)(x::SelfStar2, y::Multivector2)
    Multivector2(x.c - y.l.l, -y.l.r, -y.r.l, x.c - y.r.r)
end

# Subtraction of self-star-conjugate and a real number
function (-)(z::SelfStar2, a::Real)
    Multivector2(z.c - a, zero(a), zero(a), z.c)
end

function (-)(a::Real, z::SelfStar2)
    Multivector2(a - z.c, zero(a), zero(a), -z.c)
end

# Wedge product of two self-star-conjugate multivectors
function (∧)(x::SelfStar2, y::SelfStar2)
    xy = x.c * y.c
    Multivector2(xy, zero(xy), zero(xy), 2xy)
end

# Wedge product of a self-star-conjugate multivector and a multivector
function (∧)(x::SelfStar2, y::Multivector2)
    Multivector2(x.c * y.l.l, zero(x.c), zero(x.c), x.c * (y.l.l + y.r.r))
end

function (∧)(x::Multivector2, y::SelfStar2)
    Multivector2(x.l * y.c, zero(y.c), zero(y.c), (x.l.l + x.r.r) * y.c)
end

function (*)(z::SelfStar2, a::Real)
    SelfStar2(a * z.c)
end

function (*)(a::Real, z::SelfStar2)
    SelfStar2(a * z.c)
end

"""
    abs(z::SelfStar2)
"""
abs(z::SelfStar2) = abs(z.c)

"""
    abs2(z::SelfStar2)
"""
abs2(z::SelfStar2) = z.c^2

function inv(z::SelfStar2)
    if z.c == zero(z.c)
        error(ZeroInverse)
    end

    AntiSelfStar2(1 / z.c)
end

function (/)(x::SelfStar2, y::SelfStar2)
    if y.c == zero(y.c)
        error(ZeroDenominator)
    end

    x.c / y.c
end

function (/)(a::Real, z::SelfStar2)
    if z.c == zero(z.c)
        error(ZeroDenominator)
    end

    AntiSelfStar2(a / z.c)
end

function (/)(z::SelfStar2, a::Real)
    if a == zero(a)
        error(ZeroDenominator)
    end

    SelfStar2(z.c / a)
end

function (\)(y::SelfStar2, x::SelfStar2)
    if y.c == zero(y.c)
        error(ZeroDenominator)
    end

    y.c \ x.c
end

function (\)(z::SelfStar2, a::Real)
    if z.c == zero(z.c)
        error(ZeroDenominator)
    end

    AntiSelfStar2(z.c \ a)
end

function (\)(a::Real, z::SelfStar2)
    if a == zero(a)
        error(ZeroDenominator)
    end
    
    SelfStar2(a \ z.c)
end
