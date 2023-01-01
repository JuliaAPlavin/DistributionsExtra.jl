using DistributionsExtra
using DistributionsExtra: pred_to_intervals, IntervalsUnion
using ChainedFixes
using Accessors
using Test

@testset begin
    @test pred_to_intervals(∈(0±5)) == Interval{:closed, :closed}(-5, 5)
    @test pred_to_intervals(∉(0±5)) == IntervalsUnion((Interval{:closed, :open}(-Inf, -5), Interval{:open, :closed}(5, Inf)))
    @test pred_to_intervals(>(10) ⩔ <(2) ⩔ ∈(4..5)) == IntervalsUnion((Interval{:open, :closed}(10, Inf), Interval{:closed, :open}(-Inf, 2), Interval{:closed, :closed}(4, 5)))
    @test pred_to_intervals(>(10) ⩓ <(21) ⩓ >(50)) == 50..21
    @test pred_to_intervals(@optic(_^3 > 27)) == Interval{:open, :closed}(3, Inf)
    @test pred_to_intervals(@optic(_^3 >= 27)) == Interval{:closed, :closed}(3, Inf)
    @test pred_to_intervals(@optic(-_^3 <= -27)) == Interval{:closed, :closed}(3, Inf)
    @test pred_to_intervals(@optic(-_^3 ∈ -1..8)) == Interval{:closed, :closed}(-2, 1)
end


import CompatHelperLocal as CHL
CHL.@check()
