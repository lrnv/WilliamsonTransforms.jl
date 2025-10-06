@testitem "taylor() composes ForwardDiff-on-TaylorDiff" begin
    using ForwardDiff
    using WilliamsonTransforms: taylor

    # A few representative scalar functions
    fs = (
        (a, x) -> a + x,
        (a, x) -> a * x,
        (a, x) -> a * x^2,
        (a, x) -> expm1(a * x^2),
        (a, x) -> a / (1 + x^2),
    )

    # Helper: compute df/da of p-th derivative w.r.t. x at (a0, x0)
    function test_case(f; a0=0.3, x0=0.2, p=4)
        # Inner AD: TaylorDiff w.r.t x
        df(a, x, ::Val{P}) where {P} = taylor(y -> f(a, y), x, Val(P))[end] * factorial(P)

        # Outer AD: ForwardDiff w.r.t a
        g(a) = df(a, x0, Val(p))

        # Expect this not to error and produce a finite derivative
        d = ForwardDiff.derivative(g, a0)
        @test isfinite(d)

        # Cross-check via finite differences in a (no AD nesting)
        h = 1e-6
        d_fd = (g(a0 + h) - g(a0 - h)) / (2h)
        @test isapprox(d, d_fd; atol=1e-6, rtol=1e-6)
    end

    for f in fs
        test_case(f)
    end
end

@testitem "ForwardDiff-on-TaylorDiff with parameter vector" begin
    using ForwardDiff
    using WilliamsonTransforms: taylor

    # Functions of parameters θ = (a,b) and variable x
    fs = (
        (θ, x) -> θ[1]*x + θ[2]*x^3,
        (θ, x) -> θ[1]*sin(x) + θ[2]*cos(x),
        (θ, x) -> θ[1]*exp(x) + θ[2]/(1 + x^2),
        (θ, x) -> log1p(θ[1]*x) + θ[2]*x^2,
    )

    function test_case(f; θ0=(0.2, -0.4), x0=0.15, p=3)
        df(θ, x, ::Val{P}) where {P} = taylor(y -> f(θ, y), x, Val(P))[end] * factorial(P)
        g(θ) = df(θ, x0, Val(p))

        # Compare ForwardDiff gradient to finite differences
        grad = ForwardDiff.gradient(g, collect(θ0))
        @test all(isfinite, grad)

        h = 1e-6
        θv = collect(θ0)
        fd = similar(grad)
        for i in eachindex(θv)
            θp = copy(θv); θp[i] += h
            θm = copy(θv); θm[i] -= h
            fd[i] = (g(θp) - g(θm)) / (2h)
        end
        @test isapprox(grad, fd; atol=1e-6, rtol=1e-6)
    end

    for f in fs
        test_case(f)
    end
end

@testitem "Differentiate w.r.t. x after Taylor (ForwardDiff outer)" begin
    using ForwardDiff
    using WilliamsonTransforms: taylor

    # For smooth f, d/ dx of ∂^p/∂x^p f equals ∂^{p+1}/∂x^{p+1} f.
    # We compute g(x) = ∂^p/∂x^p f(a0, x) via Taylor, then ForwardDiff.derivative(g, x0)
    # and compare to ∂^{p+1}/∂x^{p+1} f(a0, x0) via Taylor.
    fs = (
        (a, x) -> a + x + x^3,
        (a, x) -> a * x^2 + sin(x),
        (a, x) -> a * exp(x) + cos(x),
        (a, x) -> a / (1 + x) + x^4,
    )

    function test_case(f; a0=0.7, x0=0.1, p=3)
        fpx(a, x, ::Val{P}) where {P} = taylor(y -> f(a, y), x, Val(P))[end] * factorial(P)
        g(x) = fpx(a0, x, Val(p))
        d1 = ForwardDiff.derivative(g, x0)

        d2 = fpx(a0, x0, Val(p+1))
        @test isapprox(d1, d2; atol=1e-8, rtol=1e-8)
    end

    for f in fs
        test_case(f)
    end
end

@testitem "Mixed: derivative w.r.t. a and x (ForwardDiff outer)" begin
    using ForwardDiff
    using WilliamsonTransforms: taylor

    # Check commuting derivatives on smooth functions: ∂/∂a ∂^p/∂x^p f = ∂^p/∂x^p ∂/∂a f
    fs = (
        (a, x) -> expm1(a * x) + x^2,
        (a, x) -> a * sin(x) + a^2 * x,
        (a, x) -> log1p(a * x) + a / (1 + x^2),
    )

    function dedual(x)
        x isa ForwardDiff.Dual ? ForwardDiff.value(x) : x
    end

    function test_case(f; a0=0.3, x0=0.2, p=2)
        fpx(a, x, ::Val{P}) where {P} = taylor(y -> f(a, y), x, Val(P))[end] * factorial(P)
        g(a) = fpx(a, x0, Val(p))
        da_then_px = ForwardDiff.derivative(g, a0)

        # px then da: first differentiate f w.r.t a (outer), then take p-th derivative w.r.t x via finite differences
        hfun(x) = ForwardDiff.derivative(a -> f(a, x), a0)
        hstep = 1e-6
        if p == 1
            px_then_da = (hfun(x0 + hstep) - hfun(x0 - hstep)) / (2hstep)
        elseif p == 2
            px_then_da = (hfun(x0 + hstep) - 2hfun(x0) + hfun(x0 - hstep)) / (hstep^2)
        else
            # fall back to comparing finite-difference pth derivative in x
            # naive stencil for demonstration
            px_then_da = (hfun(x0 + hstep) - hfun(x0 - hstep)) / (2hstep)
        end

    @test isapprox(dedual(da_then_px), dedual(px_then_da); atol=1e-3, rtol=1e-3)
    end

    for f in fs
        test_case(f)
    end
end
