struct PiecewiseUniform{EV<:AbstractVector,DNP} <: ContinuousUnivariateDistribution
    edges::EV
    dnp::DNP
end

PiecewiseUniform(; kwargs...) = PiecewiseUniform(values(kwargs))

PiecewiseUniform((;edges, probs)::NamedTuple{(:edges,:probs)}) =
    PiecewiseUniform(edges, DiscreteNonParametric(Base.OneTo(length(edges) - 1), probs))
PiecewiseUniform((;edges, relprobs)::NamedTuple{(:edges,:relprobs)}) =
    PiecewiseUniform(;edges, probs=relprobs ./ sum(relprobs))
PiecewiseUniform((;edges, reldensities)::NamedTuple{(:edges,:reldensities)}) =
    PiecewiseUniform(;edges, relprobs=reldensities .* diff(edges))

function Base.rand(rng::AbstractRNG, d::PiecewiseUniform)
    i = rand(rng, d.dnp)
    return rand(rng, Uniform(d.edges[i], d.edges[i+1]))
end

function Distributions.pdf(d::PiecewiseUniform, x::Real)
    i = _findi(d, x)
    isnothing(i) && return zero(promote_type(eltype(d.edges), typeof(x)))
    pi = d.dnp.p[i]
    return pi * pdf(Uniform(d.edges[i], d.edges[i+1]), x)
end

Distributions.logpdf(d::PiecewiseUniform, x::Real) = log(pdf(d, x))

function Distributions.cdf(d::PiecewiseUniform, x::Real)
    first(d.edges) ≤ x || return zero(promote_type(eltype(d.edges), typeof(x)))
    x ≤ last(d.edges)  || return one(promote_type(eltype(d.edges), typeof(x)))
    i = _findi(d, x)
    cdf(d.dnp, i - 1) + pdf(d.dnp, i) * cdf(Uniform(d.edges[i], d.edges[i+1]), x)
end

@inline function _findi(d::PiecewiseUniform, x::Real)
    first(d.edges) ≤ x ≤ last(d.edges) || return nothing
    i = findlast(≤(x), d.edges)  # searchsortedlast(d.edges, x)
    i = i == lastindex(d.edges) ? i - 1 : i  # does -1 when x == last(d.edges)
    return i
end
