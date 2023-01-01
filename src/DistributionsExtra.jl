module DistributionsExtra

using Reexport
@reexport using Distributions
@reexport using IntervalSets
using ChainedFixes
using DataPipes
using AccessorsExtra
using InverseFunctions

export ℙ, ⩔, ⩓


ℙ(f::Base.Fix2{typeof(< )}, d) = cdf(d, f.x - eps(float(typeof(f.x))))
ℙ(f::Base.Fix2{typeof(<=)}, d) = cdf(d, f.x)
ℙ(f::Base.Fix2{typeof(==)}, d::DiscreteDistribution) = pdf(d, f.x)
ℙ(f::Base.Fix2{typeof(!=)}, d::DiscreteDistribution) = 1 - pdf(d, f.x)
ℙ(f::Base.Fix2{typeof(>=)}, d) = ccdf(d, f.x - eps(float(typeof(f.x))))
ℙ(f::Base.Fix2{typeof(> )}, d) = ccdf(d, f.x)

ℙ(f::Base.Fix2{typeof(∈)}, d) = d isa DiscreteDistribution ? sum(x -> pdf(d, x), f.x) : error("Only discrete distributions supported")
function ℙ(f::Base.Fix2{typeof(∈), <:Interval}, d)
    Pr = isrightclosed(f.x) ? ℙ(<=(rightendpoint(f.x)), d) : ℙ(<(rightendpoint(f.x)), d)
    Pl = isleftclosed(f.x)  ? ℙ(<(leftendpoint(f.x)), d)   : ℙ(<=(leftendpoint(f.x)), d)
    return Pr - Pl
end

ℙ(f, d) = @p pred_to_intervals(f) |> intervals |> sum(int -> ℙ(∈(int), d))


include("intervalsunion.jl")
include("preimage.jl")


pred_to_intervals(f::Base.Fix2{typeof(> )}) = Interval{:open,   :closed}( f.x, Inf)
pred_to_intervals(f::Base.Fix2{typeof(>=)}) = Interval{:closed, :closed}( f.x, Inf)
pred_to_intervals(f::Base.Fix2{typeof(==)}) = Interval{:closed, :closed}( f.x, f.x)
pred_to_intervals(f::Base.Fix2{typeof(<=)}) = Interval{:closed, :closed}(-Inf, f.x)
pred_to_intervals(f::Base.Fix2{typeof(< )}) = Interval{:closed,   :open}(-Inf, f.x)
pred_to_intervals(f::Base.Fix2{typeof(∈), <:Union{Interval, IntervalsUnion}}) = f.x

pred_to_intervals(f::Base.Fix2{typeof(∉), <:Union{Interval, IntervalsUnion}}) = pred_to_intervals(!∈(f.x))
pred_to_intervals(f::Base.Fix2{typeof(!=)}) = pred_to_intervals(!(==)(f.x))

pred_to_intervals(f::ChainedFixes.Or) =
    @p ChainedFixes.getargs(f) |> map(IntervalsUnion(pred_to_intervals(_))) |> reduce(∪)
pred_to_intervals(f::ChainedFixes.And) =
	@p ChainedFixes.getargs(f) |> map(pred_to_intervals) |> reduce(∩)

pred_to_intervals(f::ComposedFunction{typeof(!)}) = setdiff(-Inf..Inf, IntervalsUnion(pred_to_intervals(f.inner)))
if !(!identity isa ComposedFunction)  # Julia pre-1.9
    pred_to_intervals(f::typeof(!identity).name.wrapper) = pred_to_intervals((!) ∘ f.f)
end

function pred_to_intervals(f::ComposedFunction)
	ints = pred_to_intervals(f.outer)
	preimage(f.inner, ints)
end

end
