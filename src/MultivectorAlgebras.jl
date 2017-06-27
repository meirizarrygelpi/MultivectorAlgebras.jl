__precompile__()

module MultivectorAlgebras

importall Base

"Error when finding the inverse of zero."
const ZeroInverse = "inverse of zero"

"Error when finding the inverse of a zero divisor."
const ZeroDivisorInverse = "inverse of zero divisor"

"Error when the denominator in a quotient is zero."
const ZeroDenominator = "denominator is zero"

"Error when the denominator in a quotient is a zero divisor."
const ZeroDivisorDenominator = "denominator is zero divisor"

"""
    AbstractMultivector{T <: Real} <: Number

An abstract low-dimensional multivector.
"""
abstract type AbstractMultivector{T <: Real} <: Number end

include("Multivector1.jl")
include("Multivector2.jl")
include("Multivector3.jl")
# include("Multivector4.jl")

export AbstractMultivector
export Multivector1, unreal, asarray, cloak, dagger, star, selfstar, antiselfstar, iszerodivisor
export Multivector2, commutator, crossratioL, crossratioR, möbiusL, möbiusR
export Multivector3, associator

end
