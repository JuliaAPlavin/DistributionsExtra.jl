uniform(x::AbstractInterval) = Uniform(x)
uniform(x::AbstractUnitRange) = DiscreteUniform(first(x), last(x))
uniform(x::AbstractVector) = DiscreteNonParametric(x, fill(1/length(x), length(x)))
