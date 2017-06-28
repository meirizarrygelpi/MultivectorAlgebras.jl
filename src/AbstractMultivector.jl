"""
    AbstractMultivector{T <: Real} <: Number

An abstract low-dimensional multivector.
"""
abstract type AbstractMultivector{T <: Real} <: Number end

function isreal(z::AbstractMultivector)
    iszero(unreal(z))
end

function asarray(z::AbstractMultivector)
    vcat(real(z), unreal(z))
end

function iszero(z::AbstractMultivector)
    iszero(asarray(z))
end

(+)(z::AbstractMultivector{T}) where {T <: Real} = z

"""
    abs(z::AbstractMultivector)

The absolute value of the real part of z.
"""
abs(z::AbstractMultivector) = abs(real(z))

"""
    abs2(z::AbstractMultivector)

The square of the real part of z.
Also known as the quadrance of z.
"""
abs2(z::AbstractMultivector) = (real(z))^2

"""
    iszerodivisor(z::AbstractMultivector)

Returns true if the real part of z is zero.
Note that it follows that also ``z * z = 0``.
That is, if z is a zero-divisor, then it is also nilpotent.
"""
iszerodivisor(z::AbstractMultivector) = iszero(real(z))

"""
    commutator(x::AbstractMultivector, y::AbstractMultivector)
Measures the failure of commutativity of multivector multiplication.
"""
function commutator(x::AbstractMultivector, y::AbstractMultivector)
    (x * y) - (y * x)
end

"""
    associator(x::AbstractMultivector, y::AbstractMultivector, z::AbstractMultivector)

Measures the failure of associativity of multivector multiplication.
"""
function associator(x::AbstractMultivector, y::AbstractMultivector, z::AbstractMultivector)
    ((x * y) * z) - (x * (y * z))
end

"""
    jacobiator(x::AbstractMultivector, y::AbstractMultivector, z::AbstractMultivector)

Measures the failure of multivector multiplication to satisfy the Jacobi identity.
"""
function jacobiator(x::AbstractMultivector, y::AbstractMultivector, z::AbstractMultivector)
    (x * (y * z)) + (y * (z * x)) + (z * (x * y))
end

"""
    inv(z::AbstractMultivector)

Inverse of a multivector. An error occurs if argument is a zero divisor.
"""
function inv(z::AbstractMultivector)
    if iszerodivisor(z)
        error(ZeroDivisorInverse)
    end

    conj(z) / abs2(z)
end

function (/)(a::Real, z::AbstractMultivector)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    a * inv(z)
end

function (\)(z::AbstractMultivector, a::Real)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    inv(z) * a
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
    ((x - y) \ (w - y)) * ((w - z) \ (x - z))
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
    ((w - y) / (x - y)) * ((x - z) / (w - z))
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
    ((z * c) + d) \ ((z * a) + b)
end

function möbiusL(z::AbstractMultivector, a::Real, b::Real, c::Real, d::Real)
    ((c * z) + d) \ ((a * z) + b)
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
    ((a * z) + b) / ((c * z) + d)
end

function möbiusR(z::AbstractMultivector, a::Real, b::Real, c::Real, d::Real)
    ((a * z) + b) / ((c * z) + d)
end
