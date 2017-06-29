using MultivectorAlgebras
using Base.Test: @test, @test_throws

# Multivector1:same typed arguments
@test begin
    l = Multivector1(1, 2)
    r = Multivector1{Int64}(1, 2)
    l == r
end

# Multivector1: mixed typed arguments
@test begin
    l = Multivector1(1, 2.3)
    r = Multivector1{Float64}(1.0, 2.3)
    l == r
end

# Multivector1: single argument
@test begin
    l = Multivector1(2.3)
    r = Multivector1{Float64}(2.3, 0.0)
    l == r
end

# real
@test begin
    l = real(Multivector1(1, 2))
    r = 1
    l == r
end

# unreal
@test begin
    l = unreal(Multivector1(3, 4))
    r = 4
    l == r
end

# zero: zero-like
@test begin
    l = zero(Multivector1(1, 2))
    r = Multivector1(0, 0)
    l == r
end

# zero: type argument
@test begin
    l = zero(Multivector1{Float64})
    r = Multivector1(0.0, 0.0)
    l == r
end

# one: one-like
@test begin
    l = one(Multivector1(1, 2))
    r = Multivector1(1, 0)
    l == r
end

# one: type argument
@test begin
    l = one(Multivector1{Float64})
    r = Multivector1(1.0, 0.0)
    l == r
end

# conj
@test begin
    l = conj(Multivector1(5, 6))
    r = Multivector1(5, -6)
    l == r
end

# Distributivity of conj
@test begin
    a = Multivector1(1, 2)
    b = Multivector1(3, 4)
    l = conj(a * b)
    r = conj(a) * conj(b)
    l == r
end

# cloak
@test begin
    l = cloak(Multivector1(5, 6))
    r = Multivector1(-5, 6)
    l == r
end

# Anti-distributivity of cloak
@test begin
    a = Multivector1(1, 2)
    b = Multivector1(3, 4)
    l = -cloak(a * b)
    r = cloak(a) * cloak(b)
    l == r
end

# cloak: equivalent to -dagger
@test begin
    a = Multivector1(7, 8)
    l = cloak(a)
    r = -dagger(a)
    l == r
end

# dagger
@test begin
    l = dagger(Multivector1(1, 2))
    r = Multivector1(1, -2)
    l == r
end

# Distributivity of dagger
@test begin
    a = Multivector1(1, 2)
    b = Multivector1(3, 4)
    l = dagger(a * b)
    r = dagger(a) * dagger(b)
    l == r
end

# dagger: equivalent to conj
@test begin
    a = Multivector1(9, 8)
    l = dagger(a)
    r = conj(a)
    l == r
end

# star
@test begin
    l = star(Multivector1(2, 3))
    r = Multivector1(3, 2)
    l == r
end

# Involutivity of conj
@test begin
    l = conj(conj(Multivector1(1, 2)))
    r = Multivector1(1, 2)
    l == r
end

# Involutivity of cloak
@test begin
    l = cloak(cloak(Multivector1(1, 2)))
    r = Multivector1(1, 2)
    l == r
end

# Involutivity of dagger
@test begin
    l = dagger(dagger(Multivector1(1, 2)))
    r = Multivector1(1, 2)
    l == r
end

# Involutivity of star
@test begin
    l = star(star(Multivector1(1, 2)))
    r = Multivector1(1, 2)
    l == r
end

# +
@test begin
    l = Multivector1(1, 2) + Multivector1(3, 4)
    r = Multivector1(4, 6)
    l == r
end

# +: addition of real
@test begin
    l = 2 + Multivector1(8, 9)
    r = Multivector1(10, 9)
    l == r
end

# +: addition of real
@test begin
    l = Multivector1(8, 9) + 1
    r = Multivector1(9, 9)
    l == r
end

# zero is addition identity
@test begin
    l = Multivector1(5, 6) + zero(Multivector1(5, 6))
    r = Multivector1(5, 6)
    l == r
end

# Commutativity of addition
@test begin
    l = Multivector1(1, 2) + Multivector1(3, 4)
    r = Multivector1(3, 4) + Multivector1(1, 2)
    l == r
end

# Associativity of addition
@test begin
    l = (Multivector1(1, 2) + Multivector1(3, 4)) + Multivector1(5, 6)
    r = Multivector1(1, 2) + (Multivector1(3, 4) + Multivector1(5, 6))
    l == r
end

# -
@test begin
    l = Multivector1(1, 2) - Multivector1(3, 4)
    r = Multivector1(-2, -2)
    l == r
end

# -: subtraction from real
@test begin
    l = 2 - Multivector1(8, 9)
    r = Multivector1(-6, -9)
    l == r
end

# -: subtraction of real
@test begin
    l = Multivector1(8, 9) - 1
    r = Multivector1(7, 9)
    l == r
end

# -: minus
@test begin
    l = -Multivector1(7, 8)
    r = Multivector1(-7, -8)
    l == r
end

# Anti-commutativity of subtraction
@test begin
    l = Multivector1(1, 2) - Multivector1(3, 4)
    r = -(Multivector1(3, 4) - Multivector1(1, 2))
    l == r
end

# Associativity of subtraction
@test begin
    l = (Multivector1(1, 2) - Multivector1(3, 4)) - Multivector1(5, 6)
    r = Multivector1(1, 2) - (Multivector1(3, 4) + Multivector1(5, 6))
    l == r
end

# *
@test begin
    l = Multivector1(2, 3) * Multivector1(4, 5)
    r = Multivector1(8, 22)
    l == r
end

# *: multiplication by real
@test begin
    l = 5 * Multivector1(2, 3)
    r = Multivector1(10, 15)
    l == r
end

# *: multiplication by real
@test begin
    l = Multivector1(3, 4) * 3
    r = Multivector1(9, 12)
    l == r
end

# one is multiplication identity
@test begin
    l = Multivector1(2, 3) * one(Multivector1(2, 3))
    r = Multivector1(2, 3)
    l == r
end

# Commutativity of multiplication
@test begin
    l = Multivector1(1, 2) * Multivector1(3, 4)
    r = Multivector1(3, 4) * Multivector1(1, 2)
    l == r
end

# Associativity of multiplication
@test begin
    l = (Multivector1(1, 2) * Multivector1(3, 4)) * Multivector1(5, 6)
    r = Multivector1(1, 2) * (Multivector1(3, 4) * Multivector1(5, 6))
    l == r
end

# Nilpotency
@test begin
    l = Multivector1(0, 8) * Multivector1(0, 8)
    r = Multivector1(0, 0)
    l == r
end

# Positivity
@test begin
    l = abs2(Multivector1(6, 7)) > 0
    r = true
    l == r
end

# Composition
@test begin
    a = Multivector1(1, 2)
    b = Multivector1(3, 4)
    l = abs2(a) * abs2(b)
    r = abs2(a * b)
    l == r
end

# \
@test begin
    a = Multivector1(1.0, 2.0)
    b = Multivector1(3.0, 4.0)
    c = a * b
    l = a \ c
    r = b
    l == r
end

# \: zero divisor error
@test_throws ErrorException begin
    Multivector1(0, 1) \ Multivector1(2, 3)
end

# \: zero divisor error
@test_throws ErrorException begin
    Multivector1(0, 1) \ 2
end

# /
@test begin
    a = Multivector1(1.0, 2.0)
    b = Multivector1(3.0, 4.0)
    c = a * b
    l = c / b
    r = a
    l == r
end

# /: zero divisor error
@test_throws ErrorException begin
    Multivector1(1, 2) / Multivector1(0, 3)
end

# /: zero divisor error
@test_throws ErrorException begin
    2 / Multivector1(0, 3)
end

# /: division by real
@test begin
    l = Multivector1(1.2, 3.4) / 5
    r = Multivector1(1.2 / 5, 3.4 / 5)
    l == r
end

# /: zero denominator
@test_throws ErrorException begin
    Multivector1(1, 2) / 0
end

# \: division by real
@test begin
    l = 5 \ Multivector1(1.2, 3.4)
    r = Multivector1(5 \ 1.2, 5 \ 3.4)
    l == r
end

# \: zero denominator
@test_throws ErrorException begin
    0 \ Multivector1(1, 2)
end

# inv
@test begin
    a = Multivector1(1 // 2, 3 // 4)
    b = Multivector1(5 // 6, 7 // 8)
    l = a / b
    r = inv(b / a)
    l == r
end

# inv: zero divisor error
@test_throws ErrorException begin
    inv(Multivector1(0, 1))
end

# isreal: false case
@test begin
    l = isreal(Multivector1(1, 2))
    r = false
    l == r
end

# isreal: true case
@test begin
    l = isreal(Multivector1(3, 0))
    r = true
    l == r
end

# asarray
@test begin
    l = asarray(Multivector1(4, 5))
    r = [4, 5]
    l == r
end

# iszero: false case
@test begin
    l = iszero(Multivector1(1, 2))
    r = false
    l == r
end

# iszero: true case
@test begin
    l = iszero(Multivector1(0, 0))
    r = true
    l == r
end

# abs: negative case
@test begin
    l = abs(Multivector1(-2, 3))
    r = 2
    l == r
end

# abs: positive case
@test begin
    l = abs(Multivector1(2, 3))
    r = 2
    l == r
end

# abs2
@test begin
    l = abs2(Multivector1(4, 5))
    r = 16
    l == r
end

# iszerodivisor: false case
@test begin
    l = iszerodivisor(Multivector1(2, 3))
    r = false
    l == r
end

# iszerodivisor: true case
@test begin
    l = iszerodivisor(Multivector1(0, 3))
    r = true
    l == r
end
