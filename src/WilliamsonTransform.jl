module WilliamsonTransform

import Distributions
import Expectations
import TaylorSeries
import Base.minimum
import SpecialFunctions
import Roots
import StatsBase


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
struct 𝒲{TX,TE}
    X::TX
    d::Int64
    E::TE
    function 𝒲(X::TX,d) where TX<:Distributions.UnivariateDistribution
        # S = support(X)
        # @assert S.lb ≥ 0 && S.ub ≤ Inf # check that X is indeed non-negative. 
        # @assert d ≥ 2 && isinteger(d) # check that d is an integer greater than 2.
        E = Expectations.expectation(X) 
        return new{typeof(X),typeof(E)}(X,d,E)
    end
end

function (ϕ::𝒲)(x)
    if x <= 0
        return 1 - Distributions.cdf(ϕ.X,0)
    else
        # We need to compute the expectation of (1 - x/X)^{d-1}
        return ϕ.E(y -> (y > x) * (1 - x/y)^(ϕ.d-1))
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
struct 𝒲₋₁{Tϕ,TF} <: Distributions.ContinuousUnivariateDistribution
    ϕ::Tϕ
    d::Int64
    F::TF
    function 𝒲₋₁(ϕ,d)
        # This gives the CDF of the random variable that is the inverse williamson d-transform of ϕ. 
        # The issue is that it might take a while to compute. 
        function F(x)
            rez = one(x)
            t_taylor = TaylorSeries.Taylor1(eltype(x),d)
            ϕ_taylor = ϕ(x + t_taylor).coeffs
            ϕ_taylor[end] = max(ϕ_taylor[end], 0)
            for k in 1:(d-1)
                rez -= (-1)^k * x^k * ϕ_taylor[k+1]
            end
            return rez
        end
        return new{typeof(ϕ),typeof(F)}(ϕ,d,F)
    end
end
function Distributions.rand(rng::Distributions.AbstractRNG, d::𝒲₋₁)
    u = rand(rng)
    Roots.find_zero(x -> (d.F(x) - u), (0, Inf))
end


end
