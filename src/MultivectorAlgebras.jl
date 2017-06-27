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
include("SelfStar1.jl")
include("AntiSelfStar1.jl")
include("Multivector2.jl")
include("SelfStar2.jl")
include("AntiSelfStar2.jl")
# include("Multivector3.jl")
# include("SelfStar3.jl")
# include("AntiSelfStar3.jl")
# include("Multivector4.jl")
# include("SelfStar4.jl")
# include("AntiSelfStar4.jl")

export AbstractMultivector, AbstractMultivector1
export Multivector1, cloak, dagger, star, selfstar, antiselfstar, ∧, iszerodivisor, crossratio, möbius
export SelfStar1
export AntiSelfStar1
export Multivector2, commutator, crossratioL, crossratioR, möbiusL, möbiusR
export SelfStar2
export AntiSelfStar2

end
