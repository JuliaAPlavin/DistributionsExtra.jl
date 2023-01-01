module DistributionsExtra

using Reexport
@reexport using Distributions
@reexport using IntervalSets
using ChainedFixes
using DataPipes

export ℙ


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

ℙ(f, d) = sum(int -> ℙ(∈(int), d), (pred_to_intervals(f)::IntervalsUnion).ints)


pred_to_intervals(f::Base.Fix2{typeof(> )}) = IntervalsUnion((Interval{:open,     :open}( f.x, Inf),))
pred_to_intervals(f::Base.Fix2{typeof(>=)}) = IntervalsUnion((Interval{:closed,   :open}( f.x, Inf),))
pred_to_intervals(f::Base.Fix2{typeof(==)}) = IntervalsUnion((Interval{:closed, :closed}( f.x, f.x),))
pred_to_intervals(f::Base.Fix2{typeof(<=)}) = IntervalsUnion((Interval{:open,   :closed}(-Inf, f.x),))
pred_to_intervals(f::Base.Fix2{typeof(< )}) = IntervalsUnion((Interval{:open,     :open}(-Inf, f.x),))
pred_to_intervals(f::Base.Fix2{typeof(∈), <:Interval}) = IntervalsUnion((f.x,))
pred_to_intervals(f::Base.Fix2{typeof(∉)}) = setdiff(-Inf..Inf, IntervalsUnion((f.x,)))

pred_to_intervals(f::ChainedFixes.Or) =
    @p ChainedFixes.getargs(f) |> map(pred_to_intervals) |> reduce(∪)
pred_to_intervals(f::ChainedFixes.And) =
	@p ChainedFixes.getargs(f) |> map(pred_to_intervals) |> reduce(∩)

pred_to_intervals(f::ComposedFunction{typeof(!)}) = setdiff(-Inf..Inf, pred_to_intervals(f.inner))
if !(!identity isa ComposedFunction)  # Julia pre-1.9
    pred_to_intervals(f::typeof(!identity).name.wrapper) = pred_to_intervals((!) ∘ f.f)
end


struct IntervalsUnion{T}
    ints::T
end

_dropempty(iu::IntervalsUnion) = IntervalsUnion(filter(!isempty, iu.ints))

_opposite_closedness(x::Symbol) = x == :closed ? :open : x == :open ? :closed : error()
_setdiff(a::Interval{La, Ra}, b::Interval{Lb, Rb}) where {La, Ra, Lb, Rb} = IntervalsUnion((
    Interval{La, _opposite_closedness(Lb)}(leftendpoint(a), leftendpoint(b)),
    Interval{_opposite_closedness(Rb), Ra}(rightendpoint(b), rightendpoint(a)),
))

Base.setdiff(a::Interval, b::IntervalsUnion) = @p b.ints |> map(_setdiff(a, _)) |> reduce(∩)
Base.:∪(a::IntervalsUnion, b::IntervalsUnion) = IntervalsUnion((a.ints..., b.ints...))
Base.:∩(a::Interval, b::IntervalsUnion) = @p b.ints |> map(a ∩ _) |> IntervalsUnion |> _dropempty
Base.:∩(a::IntervalsUnion, b::IntervalsUnion) = @p a.ints |> map(_ ∩ b) |> reduce(∪)

end
