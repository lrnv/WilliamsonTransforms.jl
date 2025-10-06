module WilliamsonTransformsForwardDiffExt

using WilliamsonTransforms
import ForwardDiff
import TaylorDiff

const Dual = ForwardDiff.Dual
const TaylorScalar = TaylorDiff.TaylorScalar

# Global methods to make ForwardDiff (outer) × TaylorDiff work.
# This is clearly type piracy ;)
Base.:+(a::Dual, b::TaylorScalar) = TaylorScalar(TaylorDiff.value(b) + a, map(p -> p * one(a), TaylorDiff.partials(b)))
Base.:+(b::TaylorScalar, a::Dual) = TaylorScalar(TaylorDiff.value(b) + a, map(p -> p * one(a), TaylorDiff.partials(b)))

Base.:-(a::Dual, b::TaylorScalar) = TaylorScalar(a - TaylorDiff.value(b), map(p -> p * one(a), TaylorDiff.partials(b)))
Base.:-(b::TaylorScalar, a::Dual) = TaylorScalar(TaylorDiff.value(b) - a, map(p -> p * one(a), TaylorDiff.partials(b)))

Base.:*(a::Dual, b::TaylorScalar) = TaylorScalar(TaylorDiff.value(b) * a, map(p -> p * a, TaylorDiff.partials(b)))
Base.:*(b::TaylorScalar, a::Dual) = TaylorScalar(TaylorDiff.value(b) * a, map(p -> p * a, TaylorDiff.partials(b)))

Base.:/(b::TaylorScalar, a::Dual) = (inva = inv(a); TaylorScalar(TaylorDiff.value(b) / a, map(p -> p * inva, TaylorDiff.partials(b))))
Base.:/(a::Dual, b::TaylorScalar) = a * inv(b)

end # module
