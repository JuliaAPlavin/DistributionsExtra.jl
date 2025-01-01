# these would be nice to upstream but Distributions devs don't like it, see https://github.com/JuliaStats/Distributions.jl/pull/1809
Distributions.Uniform(i::Interval) = Uniform(leftendpoint(i), rightendpoint(i))
Distributions.LogUniform(i::Interval) = LogUniform(leftendpoint(i), rightendpoint(i))
Distributions.truncated(d0, i::Interval) = truncated(d0, leftendpoint(i), rightendpoint(i))
Distributions.censored(d0, i::Interval) = censored(d0, leftendpoint(i), rightendpoint(i))
Base.convert(::Type{T}, ri::RealInterval) where {T<:Interval} = T(minimum(ri), maximum(ri))
(::Type{T})(d::UnivariateDistribution) where {T<:Interval} = convert(T, support(d))
