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

ℙ(f::Base.Fix2{typeof(∈)}, d::DiscreteDistribution) = sum(x -> pdf(d, x), f.x)
function ℙ(f::Base.Fix2{typeof(∈), <:Interval}, d)
    Pr = isrightclosed(f.x) ? ℙ(<=(rightendpoint(f.x)), d) : ℙ(<(rightendpoint(f.x)), d)
    Pl = isleftclosed(f.x)  ? ℙ(<(leftendpoint(f.x)), d)   : ℙ(<=(leftendpoint(f.x)), d)
    return Pr - Pl
end

ℙ(f, d) = sum(int -> ℙ(∈(int), d), pred_to_intervals(f))


pred_to_intervals(f::Base.Fix2{typeof(> )}) = (Interval{:open,     :open}( f.x, Inf),)
pred_to_intervals(f::Base.Fix2{typeof(>=)}) = (Interval{:closed,   :open}( f.x, Inf),)
pred_to_intervals(f::Base.Fix2{typeof(==)}) = (Interval{:closed, :closed}( f.x, f.x),)
pred_to_intervals(f::Base.Fix2{typeof(<=)}) = (Interval{:open,   :closed}(-Inf, f.x),)
pred_to_intervals(f::Base.Fix2{typeof(< )}) = (Interval{:open,     :open}(-Inf, f.x),)
pred_to_intervals(f::Base.Fix2{typeof(∈), <:Interval}) = (f.x,)
pred_to_intervals(f::Base.Fix2{typeof(∉), <:Interval{L, R}}) where {L, R} = (
    Interval{:open, _opposite_closedness(L)}(-Inf, leftendpoint(f.x)),
    Interval{_opposite_closedness(R), :open}(rightendpoint(f.x), Inf),
)

_opposite_closedness(x::Symbol) = x == :closed ? :open : x == :open ? :closed : error()

pred_to_intervals(f::ChainedFixes.Or) =
    @p ChainedFixes.getargs(f) |>
        map(a -> pred_to_intervals(a)) |>
        reduce((a, b) -> (a..., b...))

pred_to_intervals(f::ChainedFixes.And) =
	@p ChainedFixes.getargs(f) |>
		map(__ -> reduce(∩, pred_to_intervals(__))) |>
		reduce(∩) |>
		(__,)

# should use setdiff, waiting for https://github.com/JuliaMath/IntervalSets.jl/issues/106
# pred_to_intervals(f::ComposedFunction{typeof(!)}) =
#     @p pred_to_intervals(f.inner) |>

end
