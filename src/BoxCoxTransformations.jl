module BoxCoxTransformations

using Distributions: Chisq, quantile
using Optim: optimize, minimizer, UnivariateOptimizationResults
using BoxCoxTrans
using StatsAPI
import StatsAPI: confint, loglikelihood, fit

export BoxCoxTransformation,
    confint,
    transform,
    loglikelihood

include("boxcox.jl")


end
