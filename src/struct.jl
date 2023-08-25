struct BoxCoxTransformation{T <: AbstractFloat, N <: Real}
    value::T
    details::UnivariateOptimizationResults
    method::Symbol
    dat::Vector{N}
end

Base.propertynames(::BoxCoxTransformation) =
            (:lambda, :dat, :method, :details, :loglik)

function Base.getproperty(bc::BoxCoxTransformation, s::Symbol)
    if s === :lambda || s === :λ
        return bc.value
    elseif s === :𝐱
        return bc.dat
    elseif s === :loglik
        return log_likelihood(bc)
    else
        return getfield(bc, s)
    end
end

"""
    log_likelihood(bc::BoxCoxTransformation)
    log_likelihood(bc::BoxCoxTransformation, λ)

Returns log likelihood for optimal λ as specified in bc or for a particular λ,
if specified.
"""
log_likelihood(bc::BoxCoxTransformation) = BoxCoxTrans.log_likelihood(bc.dat, bc.λ; method = bc.method)
log_likelihood(bc::BoxCoxTransformation, λ) = BoxCoxTrans.log_likelihood(bc.dat, λ; method = bc.method)


"""
    transform(bc::BoxCoxTransformation; α = 0, scaled = false)
    transform(bc::BoxCoxTransformation, λ::Real; α = 0, scaled = false)

Returns transformed data for optimal λ as specified in bc or for a particular λ,
if specified.

Keyword arguments:
- α: added to all values in 𝐱 before transformation. Default = 0.
- scaled: scale transformation results.  Default = false.

"""
transform(bc::BoxCoxTransformation; α = 0, scaled = false) =
    BoxCoxTrans.transform(bc.𝐱, bc.λ; α, scaled)
transform(bc::BoxCoxTransformation, λ::Real; α = 0, scaled = false) =
    BoxCoxTrans.transform(bc.𝐱, λ; α, scaled)

"""
    BoxCoxTransformation(𝐱::Vector{<:Real}; interval = (-2.0, 2.0),
        method = :geomean, kwargs...)

Returns BoxCoxTransformation.

Lambda parameter is calculated from 𝐱 using a log-likelihood estimator.

Keyword arguments:
- method: for the cacluation of the log-likelihood
    - :geomean => -N / 2.0 * log(2 * π * σ² / gm ^ (2 * (λ - 1)) + 1)
    - :normal => -N / 2.0 * log(σ²) + (λ - 1) * sum(log.(𝐱))
- any other keyword arguments accepted by Optim.optimize function e.g. abs_tol
"""
function BoxCoxTransformation(𝐱::Vector{<:Real};
    interval = (-2.0, 2.0), method = :geomean, kwargs...)
    bc = BoxCoxTrans.lambda(𝐱; interval, method, kwargs...)
    return BoxCoxTransformation(bc.value, bc.details, method, 𝐱)
end

"""
    confint(bc::BoxCoxTransformation; alpha=0.05)

Returns confidence interval for the estimated maximum power parameter λ of a
given BoxCoxTransformation.

Note: The confidence interval is the range of λs for which the
log-likelihood(λ) >= log_likelihood(λ_max) - 0.5 * quantile(Chisq(1), 1-alpha)

Reference: Linear Models With R by Julian Faraway Section 7.1.
"""
function confint(bc::BoxCoxTransformation; alpha = 0.05)
    interval = (bc.details.initial_lower, bc.details.initial_upper)
    τ = log_likelihood(bc, bc.λ) - 0.5 * quantile(Chisq(1), 1 - alpha) # target log-likelihood
    dist_fnc(x) = abs(log_likelihood(bc, x) - τ)
    res_upper = optimize(dist_fnc, bc.λ, interval[2])
    res_lower = optimize(dist_fnc, interval[1], bc.λ)
    return (minimizer(res_lower), minimizer(res_upper))
end


function Base.show(io::IO, ::MIME"text/plain", x::BoxCoxTransformation)
	println(io, "BoxCoxTransformation: λ = " * string(x.value))
end;