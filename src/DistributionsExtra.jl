module DistributionsExtra

using Reexport
@reexport using Distributions
@reexport using IntervalSets
@reexport using IntervalUnions
using IntervalUnions: intervals
using DataPipes
@reexport using AccessorsExtra
using AccessorsExtra: getproperties
using InverseFunctions
using QuadGK
using Distributions: AbstractRNG
using LinearAlgebra: normalize!

export
    ℙ, quadgk,
    PiecewiseUniform, SphereUniformArea, SphereUniformLonLat

include("intervals_integrations.jl")
include("distributions/piecewise.jl")
include("distributions/sphereuniform.jl")


"""    ℙ(pred, dist)

Compute the probability of a predicate being true under a given distribution.
In mathematical terms, it computes `ℙ(pred(x))` given that `x ~ dist`.

`ℙ` is "exact" in the sense that it doesn't perform numerical integration.
Instead, the predicate is transformed into a union of intervals, and their probabilities are computed with `cdf`/`ccdf`/`pdf` functions.

# Examples
```
julia> ℙ(>(0), Normal(0, 1))
0.5

# predicate has to be an introspectable function composed from elements like Base.Fix2
# for example, >(0) works but anonymous x -> x > 0 doesn't:
julia> ℙ(x -> x > 0, Normal(0, 1))
ERROR: MethodError ...

julia> ℙ(>(2) ⩔ <(-2), Normal(0, 1))
0.04550026309183032

julia> ℙ((@o abs(_) > 2), Normal(0, 1))
0.04550026309183032

julia> ℙ(∈(1..5), Normal(0, 1))
0.15865497088552993

julia> ℙ(∈(0:2), Poisson(1))
0.9196986029286058
```
"""
function ℙ end

ℙ(f, d; method) = ℙ(f, d, method)
ℙ(f, d, ::typeof(cdf)) = ℙ(f, d)

ℙ(f::Base.Fix2{typeof(< )}, d::UnivariateDistribution) = cdf(d, f.x - √eps(float(typeof(f.x))))
ℙ(f::Base.Fix2{typeof(<=)}, d::UnivariateDistribution) = cdf(d, f.x)
ℙ(f::Base.Fix2{typeof(==)}, d::DiscreteUnivariateDistribution) = pdf(d, f.x)
ℙ(f::Base.Fix2{typeof(!=)}, d::DiscreteUnivariateDistribution) = 1 - pdf(d, f.x)
ℙ(f::Base.Fix2{typeof(>=)}, d::UnivariateDistribution) = ccdf(d, f.x - √eps(float(typeof(f.x))))
ℙ(f::Base.Fix2{typeof(> )}, d::UnivariateDistribution) = ccdf(d, f.x)

ℙ(f::Base.Fix2{typeof(∈)}, d::UnivariateDistribution) = d isa DiscreteUnivariateDistribution ? sum(x -> pdf(d, x), f.x) : error("Only discrete distributions supported")
function ℙ(f::Base.Fix2{typeof(∈), <:Interval}, d::UnivariateDistribution)
    Pr = isrightclosed(f.x) ? ℙ(<=(rightendpoint(f.x)), d) : ℙ(<(rightendpoint(f.x)), d)
    Pl = isleftclosed(f.x)  ? ℙ(<(leftendpoint(f.x)), d)   : ℙ(<=(leftendpoint(f.x)), d)
    return Pr - Pl
end

ℙ(f, d::UnivariateDistribution) = @p pred_to_intervals(f) |> intervals |> sum(int -> ℙ(∈(int), d))


ℙ(f, d, method::Base.Fix2{typeof(rand)}) = mean(f, method(d))

function ℙ(f, d::ContinuousUnivariateDistribution, method::AccessorsExtra.FixArgs{typeof(quadgk)})
    @assert method.args === (AccessorsExtra.Placeholder(),)
    val, err = quadgk(x -> pdf(d, x) * f(x), extrema(d)...; method.kwargs...)
    @assert err < 1e-2 * val
    val
end

include("preimage.jl")

pred_to_intervals(f::Base.Fix2{typeof(> )}) = Interval{:open,   :closed}( f.x, Inf)
pred_to_intervals(f::Base.Fix2{typeof(>=)}) = Interval{:closed, :closed}( f.x, Inf)
pred_to_intervals(f::Base.Fix2{typeof(==)}) = Interval{:closed, :closed}( f.x, f.x)
pred_to_intervals(f::Base.Fix2{typeof(<=)}) = Interval{:closed, :closed}(-Inf, f.x)
pred_to_intervals(f::Base.Fix2{typeof(< )}) = Interval{:closed,   :open}(-Inf, f.x)
pred_to_intervals(f::Base.Fix2{typeof(∈), <:Union{Interval, IntervalUnion}}) = f.x

pred_to_intervals(f::Base.Fix2{typeof(∉), <:Union{Interval, IntervalUnion}}) = pred_to_intervals(!∈(f.x))
pred_to_intervals(f::Base.Fix2{typeof(!=)}) = pred_to_intervals(!(==)(f.x))

pred_to_intervals(f::AccessorsExtra.FixArgs{typeof(isapprox), <:Tuple{AccessorsExtra.Placeholder, <:Any}, <:NamedTuple{(:atol,)}}) = f.args[2] ± f.kwargs.atol
pred_to_intervals(f::AccessorsExtra.FixArgs{typeof(isapprox), <:Tuple{<:Any, AccessorsExtra.Placeholder}, <:NamedTuple{(:atol,)}}) = f.args[1] ± f.kwargs.atol

pred_to_intervals(f::⩔) =
    @p getproperties(f) |> map(IntervalUnion(pred_to_intervals(_))) |> reduce(∪)
pred_to_intervals(f::⩓) =
    @p getproperties(f) |> map(pred_to_intervals) |> reduce(∩)

pred_to_intervals(f::ComposedFunction{typeof(!)}) = setdiff(-Inf..Inf, IntervalUnion(pred_to_intervals(f.inner)))

function pred_to_intervals(f::ComposedFunction)
    ints = pred_to_intervals(f.outer)
    preimage(f.inner, ints)
end

end
