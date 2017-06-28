"""
    Multivector1{T <: Real} <: AbstractMultivector{T}

An immutable 2-tuple of real numbers that represents a member of a 1-dimensional multivector algebra.

Each `Multivector1` has the form
```math
    a+bA
```
where ``a`` and ``b`` are real (and of the same type), and ``A * A = 0``.
Here ``*`` is the wedge product. Note that
```math
    (a + bA) * (a + bA) = a^{2} + 2abA
```
which, in general, is not equal to zero.
"""
struct Multivector1{T <: Real} <: AbstractMultivector{T}
    l::T
    r::T

    Multivector1{U}(l::U, r::U) where {U <: Real} = new(l, r)
end

Multivector1(a::T, b::T) where {T <: Real} = Multivector1{T}(a, b)
Multivector1(a::Real, b::Real) = Multivector1(promote(a, b)...)
Multivector1(a::T) where {T <: Real} = Multivector1{T}(a, zero(T))

function show(io::IO, z::Multivector1)
    print(io, "[1: ")
    print(io, z.l)
    print(io, ", A: ")
    print(io, z.r)
    print(io, "]")
end

function real(z::Multivector1)
    z.l
end

function unreal(z::Multivector1)
    z.r
end

function isreal(z::AbstractMultivector)
    iszero(unreal(z))
end

function asarray(z::AbstractMultivector)
    vcat(real(z), unreal(z))
end

function zero(z::Multivector1{T}) where {T <: Real}
    Multivector1{T}(zero(z.l), zero(z.r))
end

function zero(::Type{Multivector1{T}}) where {T <: Real}
    Multivector1{T}(zero(T), zero(T))
end

function one(z::Multivector1{T}) where {T <: Real}
    Multivector1{T}(one(z.l), zero(z.r))
end

function one(::Type{Multivector1{T}}) where {T <: Real}
    Multivector1{T}(one(T), zero(T))
end

"""
    conj{T <: Real}(z::Multivector1{T})

The `Multivector1` conjugate. If ``z=a+bA``, then `conj(z)` gives
```math
    a-bA
```
This operation is an involution.
"""
function conj(z::Multivector1{T}) where {T <: Real}
    Multivector1{T}(z.l, -z.r)
end

"""
    cloak{T <: Real}(z::Multivector1{T})

The cloak conjugate changes the sign of even blades. If ``z=a+bA``, then `cloak(z)` gives
```math
    -a+bA
```
This operation is equivalent to `-conj(z)` and thus is also an involution.
"""
function cloak(z::Multivector1{T}) where {T <: Real}
    Multivector1{T}(-z.l, z.r)
end

"""
    dagger{T <: Real}(z::Multivector1{T})

The dagger conjugate changes the sign of odd blades. If ``z=a+bA``, then `dagger(z)` gives
```math
    a-bA
```
This operation is equivalent to `conj(z)` and thus is also an involution.
"""
function dagger(z::Multivector1)
    conj(z)
end

"""
    star{T <: Real}(z::Multivector1{T})

Returns the Hodge star conjugate. If ``z=a+bA``, then `star(z)` gives
```math
    b+aA
```
This operation is an involution.
"""
function star(z::Multivector1{T}) where {T <: Real}
    Multivector1{T}(z.r, z.l)
end

"""
    selfstar(z::AbstractMultivector)

The self-star-conjugate part.
In odd number of dimensions, this operation is idempotent.
"""
function selfstar(z::AbstractMultivector)
    (z + star(z)) / 2
end

"""
    antiselfstar(z::AbstractMultivector)

The anti-self-star-conjugate part.
In odd number of dimensions, this operation is idempotent.
"""
function antiselfstar(z::AbstractMultivector)
    (z - star(z)) / 2
end

(+)(z::AbstractMultivector{T}) where {T <: Real} = z

function (+)(x::Multivector1, y::Multivector1)
    Multivector1(x.l + y.l, x.r + y.r)
end

function (+)(z::Multivector1, a::Real)
    Multivector1(z.l + a, z.r)
end

function (+)(a::Real, z::Multivector1)
    Multivector1(z.l + a, z.r)
end

function (-)(z::Multivector1{T}) where {T <: Real}
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
    (*)(x::Multivector1, y::Multivector1)

Wedge product of two 1-dimensional multivectors.
This operation is commutative and associative.
"""
function (*)(x::Multivector1, y::Multivector1)
    Multivector1(x.l * y.l, (x.r * y.l) + (y.r * x.l))
end

"""
    (*)(z::Multivector1, a::Real)
    (*)(a::Real, z::Multivector1)

Scaling and/or reflection of a `Multivector1` by a real number.
"""
function (*)(z::Multivector1, a::Real)
    Multivector1(z.l * a, z.r * a)
end

function (*)(a::Real, z::Multivector1)
    Multivector1(z.l * a, z.r * a)
end

"""
    abs(z::Multivector1)
"""
abs(z::Multivector1) = abs(z.l)

"""
    abs2(z::Multivector1)
"""
abs2(z::Multivector1) = (z.l)^2

"""
    iszerodivisor(z::Multivector1)

Returns true if `z` is of the form ``bA`` such that ``z * conj(z) = 0``.
Note that it follows that also ``z * z = 0``.
"""
iszerodivisor(z::Multivector1{T}) where {T <: Real} = iszero(z.l)

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

    Multivector1(x.l * y.l, (x.r * y.l) - (y.r * x.l)) / abs2(y)
end

function (/)(a::Real, z::Multivector1)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    a * inv(z)
end

function (/)(z::Multivector1, a::Real)
    if iszero(a)
        error(ZeroDenominator)
    end

    Multivector1(z.l / a, z.r / a)
end

function (\)(y::Multivector1, x::Multivector1)
    if iszerodivisor(y)
        error(ZeroDivisorDenominator)
    end

    Multivector1(x.l * y.l, (x.r * y.l) - (y.r * x.l)) / abs2(y)
end

function (\)(z::Multivector1, a::Real)
    if iszerodivisor(z)
        error(ZeroDivisorDenominator)
    end

    a * inv(z)
end

function (\)(a::Real, z::Multivector1)
    if iszero(a)
        error(ZeroDenominator)
    end
    
    Multivector1(z.l / a, z.r / a)
end
