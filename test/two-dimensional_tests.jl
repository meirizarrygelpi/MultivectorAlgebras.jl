using MultivectorAlgebras
using Base.Test: @test, @test_throws

# Multivector2: same typed arguments
@test begin
    l = Multivector2(1.2, 3.4, 5.6, 7.8)
    r = Multivector2{Float64}(Multivector1{Float64}(1.2, 3.4), Multivector1{Float64}(5.6, 7.8))
    l == r
end

# Multivector2: mixed typed arguments
@test begin
    l = Multivector2(1, 2.3, 4, 5.6)
    r = Multivector2{Float64}(Multivector1{Float64}(1.0, 2.3), Multivector1{Float64}(4.0, 5.6))
    l == r
end

# Multivector2: single argument
@test begin
    l = Multivector2(2.3)
    r = Multivector2{Float64}(Multivector1{Float64}(2.3, 0.0), Multivector1{Float64}(0.0, 0.0))
    l == r
end

# Multivector2: two one-dimensionals
@test begin
    l = Multivector2(Multivector1(1, 2), Multivector1(3, 4))
    r = Multivector2(1, 2, 3, 4)
    l == r
end

# Multivector2: one one-dimensionals
@test begin
    l = Multivector2(Multivector1(1, 2))
    r = Multivector2(1, 2, 0, 0)
    l == r
end

# Multivector2: three arguments
@test begin
    l = Multivector2(1, 2, 3)
    r = Multivector2(1, 2, 3, 0)
    l == r
end

# real
@test begin
    l = real(Multivector2(1, 2, 3, 4))
    r = 1
    l == r
end

# unreal: Multivector2
@test begin
    l = unreal(Multivector2(3, 4, 5, 6))
    r = [4, 5, 6]
    l == r
end

# zero: zero-like
@test begin
    l = zero(Multivector2(1, 2, 3, 4))
    r = Multivector2(0, 0, 0, 0)
    l == r
end

# zero: type argument
@test begin
    l = zero(Multivector2{Float64})
    r = Multivector2(0.0, 0.0, 0.0, 0.0)
    l == r
end

# one: one-like
@test begin
    l = one(Multivector2(1, 2, 3, 4))
    r = Multivector2(1, 0, 0, 0)
    l == r
end

# one: type argument
@test begin
    l = one(Multivector2{Float64})
    r = Multivector2(1.0, 0.0, 0.0, 0.0)
    l == r
end

# conj
@test begin
    l = conj(Multivector2(5, 6, 7, 8))
    r = Multivector2(5, -6, -7, -8)
    l == r
end

# Anti-distributivity of conj
@test begin
    a = Multivector2(1, 2, 3, 4)
    b = Multivector2(5, 6, 7, 8)
    l = conj(a * b)
    r = conj(b) * conj(a)
    l == r
end

# cloak
@test begin
    l = cloak(Multivector2(5, 6, 7, 8))
    r = Multivector2(-5, 6, 7, -8)
    l == r
end

# Anti-distributivity of cloak
@test begin
    a = Multivector2(1, 2, 5, 6)
    b = Multivector2(3, 4, 7, 8)
    l = -cloak(a * b)
    r = cloak(a) * cloak(b)
    l == r
end

# cloak: equivalent to -dagger
@test begin
    a = Multivector2(7, 8, 5, 6)
    l = cloak(a)
    r = -dagger(a)
    l == r
end

# dagger
@test begin
    l = dagger(Multivector2(1, 2, 3, 4))
    r = Multivector2(1, -2, -3, 4)
    l == r
end

# Anti-distributivity of dagger
@test begin
    a = Multivector2(1, 2, 5, 6)
    b = Multivector2(3, 4, 7, 8)
    l = dagger(a * b)
    r = dagger(a) * dagger(b)
    l == r
end

# star
@test begin
    l = star(Multivector2(1, 2, 3, 4))
    r = Multivector2(4, -3, 2, 1)
    l == r
end

# Involutivity of conj
@test begin
    l = conj(conj(Multivector2(1, 2, 3, 4)))
    r = Multivector2(1, 2, 3, 4)
    l == r
end

# Involutivity of cloak
@test begin
    l = cloak(cloak(Multivector2(1, 2, 3, 4)))
    r = Multivector2(1, 2, 3, 4)
    l == r
end

# Involutivity of dagger
@test begin
    l = dagger(dagger(Multivector2(1, 2, 3, 4)))
    r = Multivector2(1, 2, 3, 4)
    l == r
end

# star twice is dagger
@test begin
    l = star(star(Multivector2(1, 2, 3, 4)))
    r = dagger(Multivector2(1, 2, 3, 4))
    l == r
end

# +
@test begin
    l = Multivector2(1, 2, 5, 6) + Multivector2(3, 4, 7, 8)
    r = Multivector2(4, 6, 12, 14)
    l == r
end

# +: addition of real
@test begin
    l = 2 + Multivector2(8, 9, 1, 2)
    r = Multivector2(10, 9, 1, 2)
    l == r
end

# +: addition of real
@test begin
    l = Multivector2(8, 9, 2, 1) + 1
    r = Multivector2(9, 9, 2, 1)
    l == r
end

# zero is addition identity
@test begin
    l = Multivector2(5, 6, 1, 2) + zero(Multivector2(5, 6, 1, 2))
    r = Multivector2(5, 6, 1, 2)
    l == r
end

# Commutativity of addition
@test begin
    l = Multivector2(1, 2, 5, 6) + Multivector2(3, 4, 7, 8)
    r = Multivector2(3, 4, 7, 8) + Multivector2(1, 2, 5, 6)
    l == r
end

# Associativity of addition
@test begin
    l = (Multivector2(1, 2, 7, 8) + Multivector2(3, 4, 7, 8)) + Multivector2(5, 6, 7, 8)
    r = Multivector2(1, 2, 7, 8) + (Multivector2(3, 4, 7, 8) + Multivector2(5, 6, 7, 8))
    l == r
end

# -
@test begin
    l = Multivector2(1, 2, 7, 8) - Multivector2(3, 4, 5, 6)
    r = Multivector2(-2, -2, 2, 2)
    l == r
end

# -: subtraction from real
@test begin
    l = 2 - Multivector2(8, 9, 7, 8)
    r = Multivector2(-6, -9, -7, -8)
    l == r
end

# -: subtraction of real
@test begin
    l = Multivector2(8, 9, 7, 8) - 1
    r = Multivector2(7, 9, 7, 8)
    l == r
end

# -: minus
@test begin
    l = -Multivector2(7, 8, 7, 8)
    r = Multivector2(-7, -8, -7, -8)
    l == r
end

# Anti-commutativity of subtraction
@test begin
    l = Multivector2(1, 2, 7, 8) - Multivector2(3, 4, 7, 8)
    r = -(Multivector2(3, 4, 7, 8) - Multivector2(1, 2, 7, 8))
    l == r
end

# Associativity of subtraction
@test begin
    l = (Multivector2(1, 2, 3, 7) - Multivector2(3, 4, 3, 7)) - Multivector2(5, 6, 3, 7)
    r = Multivector2(1, 2, 3, 7) - (Multivector2(3, 4, 3, 7) + Multivector2(5, 6, 3, 7))
    l == r
end

# *
@test begin
    l = Multivector2(2, 3, 6, 7) * Multivector2(4, 5, 8, 9)
    r = Multivector2(8, 22, 40, 40)
    l == r
end

# *: multiplication by real
@test begin
    l = 5 * Multivector2(2, 3, 4, 5)
    r = Multivector2(10, 15, 20, 25)
    l == r
end

# *: multiplication by real
@test begin
    l = Multivector2(3, 4, 5, 6) * 3
    r = Multivector2(9, 12, 15, 18)
    l == r
end

# one is multiplication identity
@test begin
    l = Multivector2(2, 3, 4, 5) * one(Multivector2(2, 3, 4, 5))
    r = Multivector2(2, 3, 4, 5)
    l == r
end

# Non-commutativity of multiplication
@test begin
    l = commutator(Multivector2(1, 2, 3, 4), Multivector2(5, 6, 7, 8))
    r = zero(Multivector2{Int64})
    l != r
end

# Associativity of multiplication
@test begin
    l = associator(Multivector2(1, 2, 3, 4), Multivector2(5, 6, 7, 8), Multivector2(9, 10, 11, 12))
    r = zero(Multivector2{Int64})
    l == r
end

# Nilpotency
@test begin
    l = Multivector2(0, 1, 2, 3) * Multivector2(0, 1, 2, 3)
    r = Multivector2(0, 0, 0, 0)
    l == r
end

# Positivity
@test begin
    l = abs2(Multivector2(6, 7, 8, 9)) > 0
    r = true
    l == r
end

# Composition
@test begin
    a = Multivector2(1, 2, 5, 6)
    b = Multivector2(3, 4, 7, 8)
    l = abs2(a) * abs2(b)
    r = abs2(a * b)
    l == r
end

# \
@test begin
    a = Multivector2(1 // 2, 3 // 4, 5 // 6, 7 // 8)
    b = Multivector2(9 // 10, 11 // 12, 13 // 14, 15 // 16)
    c = a * b
    l = a \ c
    r = b
    l == r
end

# \: zero divisor error
@test_throws ErrorException begin
    Multivector2(0, 1, 3, 4) \ Multivector2(2, 3, 9, 8)
end

# \: zero divisor error
@test_throws ErrorException begin
    Multivector2(0, 1, 2, 3) \ 2
end

# /
@test begin
    a = Multivector2(1 // 2, 3 // 4, 5 // 6, 7 // 8)
    b = Multivector2(9 // 10, 11 // 12, 13 // 14, 15 // 16)
    c = a * b
    l = c / b
    r = a
    l == r
end

# /: zero divisor error
@test_throws ErrorException begin
    Multivector2(1, 2) / Multivector2(0, 3)
end

# /: zero divisor error
@test_throws ErrorException begin
    2 / Multivector2(0, 3)
end

# /: division by real
@test begin
    l = Multivector2(1.2, 3.4) / 5
    r = Multivector2(1.2 / 5, 3.4 / 5)
    l == r
end

# /: division of real
@test begin
    l = 2 / Multivector2(3 // 4, 5 // 6)
    r = 2 * inv(Multivector2(3 // 4, 5 // 6))
    l == r
end

# /: zero denominator
@test_throws ErrorException begin
    Multivector2(1, 2, 3, 4) / 0
end

# \: division by real
@test begin
    l = 5 \ Multivector2(1.2, 3.4, 5.6, 7.8)
    r = Multivector2(5 \ 1.2, 5 \ 3.4, 5 \ 5.6, 5 \ 7.8)
    l == r
end

# \: division of real
@test begin
    l = Multivector2(3 // 4, 5 // 6, 7 // 8, 9 // 10) \ 2
    r = inv(Multivector2(3 // 4, 5 // 6, 7 // 8, 9 // 10)) * 2
    l == r
end

# \: zero denominator
@test_throws ErrorException begin
    0 \ Multivector2(1, 2, 3, 4)
end

# inv
@test begin
    a = Multivector2(1 // 2, 3 // 4, 7 // 8, 9 // 10)
    b = Multivector2(5 // 6, 7 // 8, 7 // 8, 9 // 10)
    l = a / b
    r = inv(b / a)
    l == r
end

# inv: zero divisor error
@test_throws ErrorException begin
    inv(Multivector2(0, 1, 2, 3))
end

# isreal: false case
@test begin
    l = isreal(Multivector2(1, 2, 3, 4))
    r = false
    l == r
end

# isreal: true case
@test begin
    l = isreal(Multivector2(3, 0, 0, 0))
    r = true
    l == r
end

# asarray
@test begin
    l = asarray(Multivector2(4, 5, 6, 7))
    r = [4, 5, 6, 7]
    l == r
end

# iszero: false case
@test begin
    l = iszero(Multivector2(1, 2, 3, 4))
    r = false
    l == r
end

# iszero: true case
@test begin
    l = iszero(Multivector2(0, 0, 0, 0))
    r = true
    l == r
end

# abs: negative case
@test begin
    l = abs(Multivector2(-2, 3, 5, 4))
    r = 2
    l == r
end

# abs: positive case
@test begin
    l = abs(Multivector2(2, 3, 4, 5))
    r = 2
    l == r
end

# abs2
@test begin
    l = abs2(Multivector2(4, 5, 8, 9))
    r = 16
    l == r
end

# iszerodivisor: false case
@test begin
    l = iszerodivisor(Multivector2(2, 3, 4, 5))
    r = false
    l == r
end

# iszerodivisor: true case
@test begin
    l = iszerodivisor(Multivector2(0, 3, 5, 6))
    r = true
    l == r
end

# selfstar
@test begin
    a = 1.2
    b = 3.4
    l = selfstar(Multivector2(a, 5.0, 6.0, b))
    r = Multivector2((a + b) / 2, 0.0, 0.0, (a + b) / 2)
    l == r
end

# Involutivity of selfstar
@test begin
    z = Multivector2(1 // 2, 3 // 4, 5 // 6, 7 // 8)
    l = selfstar(selfstar(z))
    r = selfstar(z)
    l == r
end

# antiselfstar
@test begin
    a = 1 // 2
    b = 3 // 4
    l = antiselfstar(Multivector2(a, 3, 4, b))
    r = Multivector2((a - b) / 2, 0, 0, (b - a) / 2)
    l == r
end

# Involutivity of antiselfstar
@test begin
    z = Multivector2(1 // 2, 3 // 4, 5 // 6, 7 // 8)
    l = antiselfstar(antiselfstar(z))
    r = antiselfstar(z)
    l == r
end

# show
@test begin
    io = IOBuffer()
    show(io, Multivector2(1, 2, 3, 4))
    l = String(take!(io))
    r = "[1: 1, A: 2, B: 3, AB: 4]"
    l == r
end

# crossratioL and möbiusR
@test begin
    a = 1
    b = 2
    c = 3
    d = 4
    x1 = Multivector2(1 // 2, 3 // 4, 5 // 6, 7 // 8)
    x2 = Multivector2(9 // 2, 10 // 4, 11 // 6, 12 // 8)
    x3 = Multivector2(13 // 2, 14 // 4, 15 // 6, 16 // 8)
    x4 = Multivector2(17 // 2, 18 // 4, 19 // 6, 20 // 8)
    y1 = möbiusR(x1, a, b, c, d)
    y2 = möbiusR(x2, a, b, c, d)
    y3 = möbiusR(x3, a, b, c, d)
    y4 = möbiusR(x4, a, b, c, d)
    l = real(crossratioL(x1, x2, x3, x4))
    r = real(crossratioL(y1, y2, y3, y4))
    l == r
end

# crossratioR and möbiusL
@test begin
    a = 1
    b = 2
    c = 3
    d = 4
    x1 = Multivector2(1 // 2, 3 // 4, 5 // 6, 7 // 8)
    x2 = Multivector2(9 // 2, 10 // 4, 11 // 6, 12 // 8)
    x3 = Multivector2(13 // 2, 14 // 4, 15 // 6, 16 // 8)
    x4 = Multivector2(17 // 2, 18 // 4, 19 // 6, 20 // 8)
    y1 = möbiusL(x1, a, b, c, d)
    y2 = möbiusL(x2, a, b, c, d)
    y3 = möbiusL(x3, a, b, c, d)
    y4 = möbiusL(x4, a, b, c, d)
    l = real(crossratioR(x1, x2, x3, x4))
    r = real(crossratioR(y1, y2, y3, y4))
    l == r
end
