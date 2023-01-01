struct IntervalsUnion{T}
    ints::T
end

Base.:(==)(a::IntervalsUnion, b::IntervalsUnion) = a.ints == b.ints

IntervalsUnion(x::Union{Interval, IntervalsUnion}) = convert(IntervalsUnion, x)
Base.convert(::Type{IntervalsUnion}, x::IntervalsUnion) = x
Base.convert(::Type{IntervalsUnion}, x::Interval) = IntervalsUnion((x,))

_dropempty(iu::IntervalsUnion) = IntervalsUnion(filter(!isempty, iu.ints))

_opposite_closedness(x::Symbol) = x == :closed ? :open : x == :open ? :closed : error()
_setdiff(a::Interval{La, Ra}, b::Interval{Lb, Rb}) where {La, Ra, Lb, Rb} = IntervalsUnion((
    Interval{La, _opposite_closedness(Lb)}(leftendpoint(a), leftendpoint(b)),
    Interval{_opposite_closedness(Rb), Ra}(rightendpoint(b), rightendpoint(a)),
))

Base.setdiff(a::Interval, b::IntervalsUnion) = @p intervals(b) |> map(_setdiff(a, _)) |> reduce(∩) |> _dropempty
Base.setdiff(a::IntervalsUnion, b::Union{Interval, IntervalsUnion}) = @p intervals(a) |> map(setdiff(_, IntervalsUnion(b))) |> reduce(∪)
Base.:∪(a::Interval, b::IntervalsUnion) = IntervalsUnion((intervals(a)..., intervals(b)...))
Base.:∪(a::IntervalsUnion, b::Union{Interval, IntervalsUnion}) = IntervalsUnion((intervals(a)..., intervals(b)...))
Base.:∩(a::Interval, b::IntervalsUnion) = @modify(i -> a ∩ i, intervals(b) |> Elements()) |> _dropempty
Base.:∩(a::IntervalsUnion, b::Interval) = @modify(i -> b ∩ i, intervals(a) |> Elements()) |> _dropempty
Base.:∩(a::IntervalsUnion, b::IntervalsUnion) = @p intervals(a) |> map(_ ∩ b) |> reduce(∪)

intervals(x::Interval) = (x,)
intervals(x::IntervalsUnion) = x.ints

Accessors.set(x::IntervalsUnion, ::typeof(intervals), v) = IntervalsUnion(v)
