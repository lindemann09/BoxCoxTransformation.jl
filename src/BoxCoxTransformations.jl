module BoxCoxTransformations

using Distributions: Chisq, quantile
using Optim: optimize, minimizer, UnivariateOptimizationResults
using BoxCoxTrans
import StatsAPI: confint

include("struct.jl")

export BoxCoxTransformation,
    confint,
    transform

end
