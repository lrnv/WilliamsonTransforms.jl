module WilliamsonTransforms

@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end WilliamsonTransforms

import Distributions
import TaylorSeries
import Base.minimum
import Roots

export 𝒲, 𝒲₋₁

"""
    𝒲(X,d)(x)

Computes the Williamson d-transform of the random variable X, taken at point x.

For a univariate non-negative random variable ``X``, with cumulative distribution function ``F`` and an integer ``d\\ge 2``, the Williamson-d-transform of ``X`` is the real function supported on ``[0,\\infty[`` given by:

```math
\\phi(t) = 𝒲_{d}(X)(t) = \\int_{t}^{\\infty} \\left(1 - \\frac{t}{x}\\right)^{d-1} dF(x) = \\mathbb E\\left( (1 - \\frac{t}{X})^{d-1}_+\\right) \\mathbb 1_{t > 0} + \\left(1 - F(0)\\right)\\mathbb 1_{t <0}
```

This function has several properties: 
    - We have that ``\\phi(0) = 1`` and ``\\phi(Inf) = 0``
    - ``\\phi`` is ``d-2`` times derivable, and the signs of its derivatives alternates : ``\\forall k \\in 0,...,d-2, (-1)^k \\phi^{(k)} \\ge 0``.
    - ``\\phi^{(d-2)}`` is convex.

These properties makes this function what is called an *archimedean generator*, able to generate *archimedean copulas* in dimensions up to ``d``. 

References: 
- Williamson, R. E. (1956). Multiply monotone functions and their Laplace transforms. Duke Math. J. 23 189–207. MR0077581
- McNeil, Alexander J., and Johanna Nešlehová. "Multivariate Archimedean copulas, d-monotone functions and ℓ 1-norm symmetric distributions." (2009): 3059-3097.
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

function taylor(f, x, d, T)
    return f(x + TaylorSeries.Taylor1(T,d)).coeffs
end

"""
    𝒲₋₁(ϕ,d)


Computes the inverse Williamson d-transform of the d-monotone archimedean generator ϕ. 

A ``d``-monotone archimedean generator is a function ``\\phi`` on ``\\mathbb R_+`` that has these three properties:
- ``\\phi(0) = 1`` and ``\\phi(Inf) = 0``
- ``\\phi`` is ``d-2`` times derivable, and the signs of its derivatives alternates : ``\\forall k \\in 0,...,d-2, (-1)^k \\phi^{(k)} \\ge 0``.
- ``\\phi^{(d-2)}`` is convex.

For such a function ``\\phi``, the inverse Williamson-d-transform of ``\\phi`` is the cumulative distribution function ``F`` of a non-negative random variable ``X``, defined by : 

```math
F(x) = 𝒲_{d}^{-1}(\\phi)(x) = 1 - \\frac{(-x)^{d-1} \\phi_+^{(d-1)}(x)}{k!} - \\sum_{k=0}^{d-2} \\frac{(-x)^k \\phi^{(k)}(x)}{k!}
```

We return this cumulative distribution function in the form of the corresponding random variable `<:Distributions.ContinuousUnivariateDistribution` from `Distributions.jl`. You may then compute : 
    - The cdf via `Distributions.cdf`
    - The pdf via `Distributions.pdf` and the logpdf via `Distributions.logpdf`
    - Samples from the distribution via `rand(X,n)`

References: 
    - Williamson, R. E. (1956). Multiply monotone functions and their Laplace transforms. Duke Math. J. 23 189–207. MR0077581
    - McNeil, Alexander J., and Johanna Nešlehová. "Multivariate Archimedean copulas, d-monotone functions and ℓ 1-norm symmetric distributions." (2009): 3059-3097.
"""
struct 𝒲₋₁{Tϕ} <: Distributions.ContinuousUnivariateDistribution
    ϕ::Tϕ
    d::Int64
    function 𝒲₋₁(ϕ,d)
        @assert ϕ(0.0) == 1.0
        @assert ϕ(float(Inf)) == 0.0
        # And assertion about d-monotony... how can this be check ? this is hard. 
        return new{typeof(ϕ)}(ϕ,d)
    end
end
function Distributions.cdf(d::𝒲₋₁, x::Real)
    rez = zero(x)
    c_ϕ = taylor(d.ϕ, x, d.d, typeof(x))
    c_ϕ[end] = max(c_ϕ[end], 0)
    for k in 0:(d.d-1)
        if c_ϕ[k+1] != 0 # We need c_ϕ = 0 to dominate x = Inf
            rez += (-1)^k * x^k * c_ϕ[k+1]
        end
    end
    # simple hack to ensure convergence :
    return isnan(rez) ? one(x) : 1 - rez
    # return 1-rez
end
function Distributions.logpdf(d::𝒲₋₁, x::Real)
    ϕ_d = taylor(d.ϕ, x, d.d+1, typeof(x))[end]
    r = (d.d-1)*log(x) - sum(log.(1:(d.d-1)))
    return log(ϕ_d) + r
end
function Distributions.rand(rng::Distributions.AbstractRNG, d::𝒲₋₁)
    u = rand(rng)
    Roots.find_zero(x -> (Distributions.cdf(d,x) - u), (0.0, Inf))
end
end
