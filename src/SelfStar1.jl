"""
    SelfStar1{T <: Real} <: AbstractMultivector1{T}

A self-star-conjugate 1-dimensional multivector.
"""
immutable SelfStar1{T <: Real} <: AbstractMultivector1{T}
    c::T

    SelfStar1{U <: Real}(c::U) = new(c)
end

SelfStar1{T <: Real}(c::T) = SelfStar1{T}(c)

Multivector1{T <: Real}(z::SelfStar1{T}) = Multivector1{T}(z.c, z.c)

function show(io::IO, z::SelfStar1)
    print(io, "[(1+W): ")
    print(io, z.c)
    print(io, "]")
end

function zero{T <: Real}(z::SelfStar1{T})
    SelfStar1{T}(zero(z.c))
end

function zero{T <: Real}(::Type{SelfStar1{T}})
    SelfStar1{T}(zero(T))
end

function one{T <: Real}(z::SelfStar1{T})
    Multivector1{T}(one(z.c), zero(z.c))
end

function one{T <: Real}(::Type{SelfStar1{T}})
    Multivector1{T}(one(T), zero(T))
end

"""
    conj(z::SelfStar1)

Turns a self-star-conjugate 1-dimensional multivector into an
anti-self-star-conjugate 1-dimensional multivector. This operation is an involution.
"""
conj{T <: Real}(z::SelfStar1{T}) = AntiSelfStar1{T}(z.c)

"""
    dagger(z::SelfStar1)

Turns a self-star-conjugate 1-dimensional multivector into an
anti-self-star-conjugate 1-dimensional multivector. This operation is an involution.
"""
dagger(z::SelfStar1) = conj(z)

"""
    star(z::SelfStar1)

The identity operation. This operation is an involution.
"""
star(z::SelfStar1) = z

# Unary +
(+)(z::SelfStar1) = z

# Addition of self-star-conjugate multivectors
function (+)(x::SelfStar1, y::SelfStar1)
    SelfStar1(x.c + y.c)
end

# Addition of a self-star-conjugate multivector and a multivector
function (+)(x::SelfStar1, y::Multivector1)
    Multivector1(x.c + y.l, x.c + y.r)
end

function (+)(x::Multivector1, y::SelfStar1)
    Multivector1(x.l + y.c, x.r + y.c)
end

# Addition of a self-star-conjugate multivector and a real number
function (+)(z::SelfStar1, a::Real)
    Multivector1(z.c + a, z.c)
end

function (+)(a::Real, z::SelfStar1)
    Multivector1(z.c + a, z.c)
end

# Unary -
function (-){T <: Real}(z::SelfStar1{T})
    SelfStar1{T}(-z.c)
end

# Subtraction of self-star-conjugate multivectors
function (-)(x::SelfStar1, y::SelfStar1)
    SelfStar1(x.c - y.c)
end

# Subtraction of self-star-conjugate multivector and a multivector
function (-)(x::Multivector1, y::SelfStar1)
    Multivector1(x.l - y.c, x.r - y.c)
end

function (-)(x::SelfStar1, y::Multivector1)
    Multivector1(x.c - y.l, x.c - y.r)
end

# Subtraction of self-star-conjugate and a real number
function (-)(z::SelfStar1, a::Real)
    Multivector1(z.c - a, z.c)
end

function (-)(a::Real, z::SelfStar1)
    Multivector1(a - z.c, -z.c)
end

# Wedge product of two self-star-conjugate multivectors
function (∧)(x::SelfStar1, y::SelfStar1)
    xy = x.c * y.c
    Multivector1(xy, 2xy)
end

# Wedge product of a self-star-conjugate multivector and a multivector
function (∧)(x::SelfStar1, y::Multivector1)
    Multivector1(x.c * y.l, x.c * (y.l + y.r))
end

function (∧)(x::Multivector1, y::SelfStar1)
    Multivector1(x.l * y.c, (x.l + x.r) * y.c)
end

function (*)(z::SelfStar1, a::Real)
    SelfStar1(a * z.c)
end

function (*)(a::Real, z::SelfStar1)
    SelfStar1(a * z.c)
end

"""
    abs(z::SelfStar1)
"""
abs(z::SelfStar1) = abs(z.c)

"""
    abs2(z::SelfStar1)
"""
abs2(z::SelfStar1) = z.c^2

function inv(z::SelfStar1)
    if z.c == zero(z.c)
        error(ZeroInverse)
    end

    AntiSelfStar1(1 / z.c)
end

function (/)(x::SelfStar1, y::SelfStar1)
    if y.c == zero(y.c)
        error(ZeroDenominator)
    end

    x.c / y.c
end

function (/)(a::Real, z::SelfStar1)
    if z.c == zero(z.c)
        error(ZeroDenominator)
    end

    AntiSelfStar1(a / z.c)
end

function (/)(z::SelfStar1, a::Real)
    if a == zero(a)
        error(ZeroDenominator)
    end

    SelfStar1(z.c / a)
end

function (\)(y::SelfStar1, x::SelfStar1)
    if y.c == zero(y.c)
        error(ZeroDenominator)
    end

    y.c \ x.c
end

function (\)(z::SelfStar1, a::Real)
    if z.c == zero(z.c)
        error(ZeroDenominator)
    end

    AntiSelfStar1(z.c \ a)
end

function (\)(a::Real, z::SelfStar1)
    if a == zero(a)
        error(ZeroDenominator)
    end
    
    SelfStar1(a \ z.c)
end
