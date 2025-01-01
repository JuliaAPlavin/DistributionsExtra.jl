module DistributionsExtra

using Reexport
@reexport using Distributions
@reexport using IntervalSets
@reexport using IntervalUnions
using IntervalUnions: intervals
using DataPipes
using AccessorsExtra
using AccessorsExtra: getproperties
using InverseFunctions
using QuadGK

export ℙ, ℙᵋ


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


function ℙᵋ(f, d::ContinuousUnivariateDistribution; kwargs...)
    val, err = quadgk(x -> pdf(d, x) * f(x), extrema(d)...; kwargs...)
    @assert err < 1e-3 * val
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
