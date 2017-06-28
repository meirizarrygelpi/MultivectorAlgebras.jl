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

include("AbstractMultivector.jl")
include("Multivector1.jl")
include("Multivector2.jl")
include("Multivector3.jl")
include("Multivector4.jl")

export AbstractMultivector, asarray, selfstar, antiselfstar, iszerodivisor, commutator, associator, jacobiator, crossratioL, crossratioR, möbiusL, möbiusR
export Multivector1, unreal, cloak, dagger, star
export Multivector2
export Multivector3
export Multivector4

end
