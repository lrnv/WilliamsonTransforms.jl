@testitem "Exemple 2.1 - Minimal Frechet (W) Generator, dimension 2" begin
    using Distributions
    d=2
    X = Dirac(1)
    ϕ(x, d) = max((1-x)^(d-1),zero(x))
    Xhat = 𝒲₋₁(x -> ϕ(x,d),d)
    ϕhat = 𝒲(X,d)

    rand(Xhat,10)
    
    @test maximum(abs.([cdf(X,x) - cdf(Xhat,x) for x in 0:0.01:10*d])) <= sqrt(eps(Float64))
    @test maximum(abs.([ϕ(x, d) - ϕhat(x) for x in 0:0.01:10])) <= sqrt(eps(Float64))
end

@testitem "Exemple 3.2 - Independence generator" begin
    using Distributions
    for d in 3:20
        X = Erlang(d)
        ϕ(x) = exp(-x)
        Xhat = 𝒲₋₁(ϕ,d)
        ϕhat = 𝒲(X,d)

        rand(Xhat,10)
    
        @test maximum(abs.([cdf(X,x) - cdf(Xhat,x) for x in 0:0.01:3*d])) <= sqrt(eps(Float64))
        @test maximum(abs.([ϕ(x) - ϕhat(x) for x in 0:0.01:10])) <= sqrt(eps(Float64))
    end
end

@testitem "Exemple 3.3: Clayton Generator" begin
    using SpecialFunctions, Distributions
    
    # exemple 3.3. : back to clayton. 
    ϕ(x, θ) = max((1 + θ * x),zero(x))^(-1/θ)
    function F(x, θ, d)
        if x < 0
            return zero(x)
        end
        α = -1/θ
        if θ < 0
            if x >= α
                return one(x)
            end
            rez = zero(x)
            x_α = x/α
            for k in 0:(d-1)
                rez += gamma(α+1)/gamma(α-k+1)/gamma(k+1) * (x_α)^k * (1 - x_α)^(α-k)
            end
            return 1-rez
        elseif θ == 0
            return exp(-x)
        else
            rez = zero(x)
            for k in 0:(d-1)
                pr = one(θ)
                for j in 0:(k-1)
                    pr *= (1+j*θ)
                end
                rez += pr / gamma(k+1) * x^k * (1 + θ * x)^(-(1/θ+k))
            end
            return 1-rez
        end
    end

    for (d, θ) in (
        (3, 1/7),
        (2, -0.2), 
        (10, -1/10),
        (2, -1.0)
    )
        Xhat = 𝒲₋₁(x -> ϕ(x,θ),d)
        rand(Xhat,10)
        @test maximum(abs.([F(x,θ,d) - cdf(Xhat,x) for x in 0:0.01:10])) <= sqrt(eps(Float64))
    end

end


@testitem "AMH Generator - test theta=-1" begin
    ϕ(t) = 2 / (1+exp(t))
    d=2
    X = 𝒲₋₁(ϕ,d)
    rand(X,100)
    @test true
end

@testitem "Quantile test - Independence Generator" begin
    using Distributions
    for d in 3:20
        X = Erlang(d)
        ϕ(x) = exp(-x)
        Xhat = 𝒲₋₁(ϕ, d)

        # Perform 10 tests of the quantile function
        for _ in 1:10
            p = rand()
            x = quantile(X, p)
            xhat = quantile(Xhat, p)

            # Verify that the difference between the quantiles is small
            @test abs(x - xhat) <= sqrt(eps(Float64))
        end
    end
end

@testitem "Expectation test - IndependantCopula" begin
    using Distributions
    for d in 3:20
        X = Erlang(d)
        ϕ(x) = exp(-x)
        Xhat = 𝒲₋₁(ϕ, d)
        truth = Distributions.expectation(ϕ, X)
        estimated = Distributions.expectation(ϕ,Xhat)
        @test truth ≈ estimated
    end
end


@testitem "testing one-dimensional williamson transformation" begin
    using Distributions
    ϕ(x) = exp(-x)
    Xhat = 𝒲₋₁(ϕ, 1)
    @assert all(Distributions.cdf(Xhat,x) == 1 - ϕ(x) for x in -log.(rand(1000)))
end


@testitem "testing fractional-dimensional williamson transformation" begin
    using Distributions
    ϕ(x) = exp(-x)
    @test_throws AssertionError 𝒲₋₁(ϕ, 0.7)
end
