struct SpherePiecewiseLatUniformArea{PU<:PiecewiseUniform,AS<:AbstractVector} <: ContinuousMultivariateDistribution
    pu::PU
    areas::AS
end

SpherePiecewiseLatUniformArea(; kwargs...) = SpherePiecewiseLatUniformArea(values(kwargs))

SpherePiecewiseLatUniformArea((;edges, probs)::NamedTuple{(:edges,:probs)}) =
    SpherePiecewiseLatUniformArea(PiecewiseUniform(;edges, probs), map(sphere_area_of_lats, edges_to_ints(edges)))
SpherePiecewiseLatUniformArea((;edges, relprobs)::NamedTuple{(:edges,:relprobs)}) =
    SpherePiecewiseLatUniformArea(;edges, probs=relprobs ./ sum(relprobs))
SpherePiecewiseLatUniformArea((;edges, reldensities)::NamedTuple{(:edges,:reldensities)}) =
    SpherePiecewiseLatUniformArea(;edges, relprobs=reldensities .* map(sphere_area_of_lats, edges_to_ints(edges)))

Base.length(d::SpherePiecewiseLatUniformArea) = 2


function Distributions.pdf(d::SpherePiecewiseLatUniformArea, x::AbstractVector{<:Real})
    lat = length(x) == 2 ? x[2] : error("")
    i = _findi(d.pu, lat)
    isnothing(i) && return zero(promote_type(eltype(d.pu.edges), typeof(lat)))
    pi = d.pu.dnp.p[i]
    return pi * 1/d.areas[i]
end

# function Base.rand(d::DecPiecewiseConstantDistribution)
# 	weights = map(d.ints_muls) do (;int, mul)
# 		sphere_area_of_lats(int) * mul
# 	end
# 	@assert sum(weights) ≈ 1
# 	i = rand(DiscreteNonParametric(eachindex(d.ints_muls), weights))
# 	int = d.ints_muls[i].int
	
#     dec = asin(1 - 2 * rand())
#     while dec ∉ int
#         dec = asin(1 - 2 * rand())
#     end
# 	return dec
# end

function sphere_area_of_lats(i::Interval)
    a, b = endpoints(i)
    2π*(sin(b) - sin(a))
end

edges_to_ints(es) = Interval.((@delete last(es)), (@delete first(es)))
