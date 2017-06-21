# MultivectorAlgebras

[![Build Status](https://travis-ci.org/meirizarrygelpi/MultivectorAlgebras.jl.svg?branch=master)](https://travis-ci.org/meirizarrygelpi/MultivectorAlgebras.jl)

[![Coverage Status](https://coveralls.io/repos/meirizarrygelpi/MultivectorAlgebras.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/meirizarrygelpi/MultivectorAlgebras.jl?branch=master)

[![codecov.io](http://codecov.io/github/meirizarrygelpi/MultivectorAlgebras.jl/coverage.svg?branch=master)](http://codecov.io/github/meirizarrygelpi/MultivectorAlgebras.jl?branch=master)

This module provides arithmetic for low-dimensional real multivector algebras (dimensions 1, 2, 3, and 4).

Three concrete types are implemented:

1. `Multivector1`: a 1-dimensional multivector, with 2 real components
1. `SelfStar1`: a 1-dimensional self-star-conjugate multivector, with 1 real component
1. `AntiSelfStar1`: a 1-dimensional anti-self-star-conjugate multivector, with 1 real component

Future concrete types are:

1. `Multivector2`: a 2-dimensional multivector, with 4 real components
1. `SelfStar2`: a 2-dimensional self-star-conjugate multivector, with 1 real component
1. `AntiSelfStar2`: a 2-dimensional anti-self-star-conjugate multivector, with 1 real component
1. `Multivector3`: a 3-dimensional multivector, with 8 real components
1. `SelfStar3`: a 3-dimensional self-star-conjugate multivector, with 4 real components
1. `AntiSelfStar3`: a 3-dimensional anti-self-star-conjugate multivector, with 4 real components
1. `Multivector4`: a 4-dimensional multivector, with 16 real components
1. `SelfStar4`: a 4-dimensional self-star-conjugate multivector, with 4 real components
1. `AntiSelfStar4`: a 4-dimensional anti-self-star-conjugate multivector, with 4 real components

Work-in-progress!