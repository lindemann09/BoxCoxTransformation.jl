# BoxCoxTransformations

[![Build Status](https://github.com/lindemann09/BoxCoxTransformations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lindemann09/BoxCoxTransformations.jl/actions/workflows/CI.yml?query=branch%3Amain)


Wraps currently merely [BoxCoxTrans.jl](https://github.com/tk3369/BoxCoxTrans.jl)
and provides confidence intervals.

```
using BoxCoxTransformations
using Distributions

x = rand(Gamma(2,2), 10000) .+ 1;

bc = BoxCoxTransformation(x)

# lambda parameter
bc.lambda

# confidence intervall for lambda
confint(bc)

# transformed data
transform(bc)
```