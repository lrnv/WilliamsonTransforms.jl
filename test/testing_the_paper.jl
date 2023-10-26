@testitem "Exemple 2.1" begin
    using Distributions
    d=10
    ϕ = 𝒲(Dirac(1),d)
    generator_ex2_2(x) = max((1-x)^(d-1),0)
    @test all(ϕ(x) ≈ generator_ex2_2(x) for x in 0:0.01:10)
end

@testitem "Exemple 3.2" begin
    using Distributions
    d=10
    ϕ = 𝒲(Erlang(10),d)
    gen_indep(x) = exp(-x)
    @test maximum(abs.([ϕ(x) - gen_indep(x) for x in 0:0.01:10])) <= sqrt(eps(Float64))
    # @test all(ϕ(x) ≈ gen_indep(x) for x in 0:0.01:10)
end


@testitem "Exemple 3.3: inverse williamson clayton" begin
    using SpecialFunctions
    
    # exemple 3.3. : back to clayton. 
    gen_clayton(x,θ) = (1 + θ * x)^(-1/θ)
    function true_radial_cdf_for_clayton(x,θ,d)
        if x < 0
            return zero(x)
        end
        if θ < 0
            α = -1/θ
            if x >= α
                return one(x)
            end
            rez = zero(x)
            θx = θ*x
            cst = log(-θx/(1+θx))
            @show x, cst
            for k in 0:(d-1)
                rez += exp(loggamma(α+k+1) - loggamma(k+1) + k*cst)
            end
            rez *= (1+θx)^α/gamma(α+1)
            return 1-rez
        elseif θ == 0
            return exp(-x)
        else
            rez = zero(x)
            for k in 0:(d-1)
                rez +=  prod(1+j*θ for j in 0:(k-1))/factorial(k) * x^k * (1+θ*x)^(-(1/θ+k))
            end
            return 1-rez
        end
    end
    θ = -0.3
    X = 𝒲₋₁(x -> gen_clayton(x,θ),2)

    @test maximum(abs.([true_radial_cdf_for_clayton(x,θ,2) - X.F(x) for x in 0:0.01:10])) <= sqrt(eps(Float64))
end






# # An easy case:
# F1(x) = (1 - exp(-x)) * (x > 0)
# X = FromCDF(F1)
# x = rand(X,10000)
# Plots.plot(t -> StatsBase.ecdf(x)(t), 0, 10)
# Plots.plot!(F1)

# # A more involved one: 
# F2(x) = 1*(x >= 2) # Dirac(1)
# X = FromCDF(F2)
# x = rand(X,10000)
# Plots.plot(t -> StatsBase.ecdf(x)(t), 0, 4)
# Plots.plot!(F2)

# # A more involved one: 
# F3(x) = Distributions.cdf(Distributions.Binomial(10,0.7),x)
# X = FromCDF(F3)
# x = rand(X,10000)
# Plots.plot(t -> StatsBase.ecdf(x)(t), 0, 12)
# Plots.plot!(F3)


# # Final try: 
# F4(x) = (F1(x)+F2(x)+F3(x) + x>=>)/3
# X = FromCDF(F4)
# x = rand(X,10000)
# Plots.plot(t -> StatsBase.ecdf(x)(t), 0, 12)
# Plots.plot!(F4)







# # ϕ = 𝒲(Gamma(1,2),10)
# ϕ = 𝒲(Dirac(1),10)
# ϕ(0.0)
# ϕ(Inf)
# using Plots
# plot(x -> ϕ(x),xlims=(0,10))


# # Now we could implement the inverse case
# # which allows to construct the random variable corresponding to a generator. 
# # some precomputations might be necessary. 
# # and a struct. 

# # a certain generator from Ex 2.1: 
# gen_ex22(x,d) = max((1-x)^(d-1),0)
# # the clayton gen from ex 2.3: 
# gen_clayton(x,d,θ) = (1 + θ * x)^(-1/θ)
# # Note: θ = 0 generates the independence copula. 
# # θ < 0 is interesting, as it is equal to gen_ex22 at the lower bound θ = -1/(d-1)

# # example 3.1 : 
# # X = Dirac(1) correspond to gen_ex22
# ϕ = 𝒲(Dirac(1),10)
# plot(x -> ϕ(x),xlims=(0,10))
# plot!(x->gen_ex22(x,10))
# # indeed same plot ! 


# # exemple 3.2: 
# # independent copula correspond to erlang distribution with parameter d
# ϕ = 𝒲(Erlang(10),10)
# gen_indep(x) = exp(-x)
# plot(x -> ϕ(x),xlims=(0,10))
# plot!(gen_indep)


# we should check if a williamson d-transform of this distribution recovers the generator correctly. 


