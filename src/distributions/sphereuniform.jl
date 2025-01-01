struct SphereUniformArea <: ContinuousMultivariateDistribution end

Base.length(::SphereUniformArea) = 3

Distributions.insupport(s::SphereUniformArea, x::AbstractVector{<:Real}) =
    if length(x) == 3
        Distributions.norm(x) ≈ 1
    elseif length(x) == 2
        true
    else
        error("Invalid dimensionality: $(length(x))")
    end

function Distributions._rand!(rng::AbstractRNG, ::SphereUniformArea, x::AbstractVector{<:Real})
    Distributions.randn!(rng, x)
    normalize!(x)
end

Distributions.pdf(s::SphereUniformArea, x::AbstractVector{<:Real}) = insupport(s, x) ? 1/4π : 0.0


struct SphereUniformLonLat <: ContinuousMultivariateDistribution end

Base.length(::SphereUniformLonLat) = 2
# Distributions.insupport(s::SphereUniformLonLat, x::AbstractVector{<:Real}) = length(x) == length(s) && Distributions.isunitvec(x)

function Distributions._rand!(rng::AbstractRNG, ::SphereUniformLonLat, x::AbstractVector{<:Real})
    x[1] = rand(rng, Uniform(-π, π))
    x[2] = asin(1 - 2 * rand(rng))
    x
end

Distributions.pdf(s::SphereUniformLonLat, x::AbstractVector{<:Real}) = cos(x[2])/4π
