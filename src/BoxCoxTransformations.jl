module BoxCoxTransformations

using Reexport
using Optim: optimize, minimizer, UnivariateOptimizationResults
using Distributions: Chisq, quantile
import StatsAPI: confint

using BoxCoxTrans


include("struct.jl")

export BoxCoxTransformation,
    confint,
    transform

end
