module WilliamsonTransform

import Distributions
import TaylorSeries
import Base.minimum
import Roots

export 𝒲, 𝒲₋₁

"""
    𝒲(X,d)(x)

Computes the Williamson d-transform of the random variable X, taken at point x, as described in 

WILLIAMSON, R. E. (1956). Multiply monotone functions and their Laplace transforms. Duke Math. J. 23 189–207. MR0077581

and used in 

McNeil, Alexander J., and Johanna Nešlehová. "Multivariate Archimedean copulas, d-monotone functions and ℓ 1-norm symmetric distributions." (2009): 3059-3097.

For a univariate non-negative random variable ``X`` for distribution function ``F`` and ``d\\ge 2`` and integer, the williamson-d-transform of ``X`` is the real function supported on ``[0,\\infty[`` given by:

```math
𝒲_{X,d}(x) = \\int_{x}^{\\infty} \\left(1 - \\frac{x}{t}\\right)^{d-1} dF(t) = \\mathbb E\\left( (1 - \\frac{x}{X})^{d-1}_+\\right) \\mathbb 1_{x > 0} + \\left(1 - F(0)\\right)\\mathbb 1_{x <0}
```
"""
struct 𝒲{TX}
    X::TX
    d::Int64
    # E::TE
    function 𝒲(X::TX,d) where TX<:Distributions.UnivariateDistribution
        @assert minimum(X) ≥ 0 && maximum(X) ≤ Inf 
        @assert d ≥ 2 && isinteger(d) 
        return new{typeof(X)}(X,d)
    end
end

function (ϕ::𝒲)(x)
    if x <= 0
        return 1 - Distributions.cdf(ϕ.X,0)
    else
        return Distributions.expectation(y -> (1 - x/y)^(ϕ.d-1) * (y > x), ϕ.X)
        # We need to compute the expectation of (1 - x/X)^{d-1}
        # return ϕ.E(y -> (y > x) * (1 - x/y)^(ϕ.d-1))
    end
end

"""
    𝒲₋₁(ϕ,d)

Computes the inverse Williamson d-transform of the d-monotone archimedean generator ϕ. This inverse is a CDF, and we return it on the form of a random variable `<:Distributions.ContinuousUnivariateDistribution` from `Distributions.jl`. The result can be sampled through `Distributions.rand()`. See 

WILLIAMSON, R. E. (1956). Multiply monotone functions and their Laplace transforms. Duke Math. J. 23 189–207. MR0077581

and moreover

McNeil, Alexander J., and Johanna Nešlehová. "Multivariate Archimedean copulas, d-monotone functions and ℓ 1-norm symmetric distributions." (2009): 3059-3097.

for details. 

The cumulative distribution function of this random variable is given by:

```math
𝒲₋₁(X,d)(x) = 1 - \\frac{(-x)^{d-1} \\phi_+^{(d-1)}(x)}{k!} - \\sum_{k=0}^{d-2} \\frac{(-x)^k \\phi^{(k)}(x)}{k!}
```
"""
function taylor(f, x, d, T)
    return f(x + TaylorSeries.Taylor1(T,d)).coeffs
end
struct 𝒲₋₁{Tϕ} <: Distributions.ContinuousUnivariateDistribution
    ϕ::Tϕ
    d::Int64
    function 𝒲₋₁(ϕ,d)
        @assert ϕ(0) == 1
        @assert ϕ(Inf) == 0
        # And assertion about d-monotony... how can this be check ? this is hard. 
        return new{typeof(ϕ)}(ϕ,d)
    end
end
function Distributions.cdf(d::𝒲₋₁, x::Real)
    rez = zero(x)
    c_ϕ = taylor(d.ϕ, x, d.d, typeof(x))
    c_ϕ[end] = max(c_ϕ[end], 0)
    for k in 0:(d.d-1)
        rez += (-1)^k * x^k * c_ϕ[k+1]
    end
    return 1-rez
end
function Distributions.logpdf(d::𝒲₋₁, x::Real)
    ϕ_d = taylor(d.ϕ, x, d.d+1, typeof(x))[end]
    r = (d.d-1)*log(x) - sum(log.(1:(d.d-1)))
    return log(ϕ_d) + r
end
function Distributions.rand(rng::Distributions.AbstractRNG, d::𝒲₋₁)
    u = rand(rng)
    Roots.find_zero(x -> (Distributions.cdf(d,x) - u), (0, Inf))
end
end
