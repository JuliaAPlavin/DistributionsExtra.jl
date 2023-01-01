using DistributionsExtra
using DistributionsExtra: pred_to_intervals, IntervalsUnion
using ChainedFixes
using Accessors
using Test

@testset begin
    @test pred_to_intervals(∈(0±5)) == -5..5
    @test pred_to_intervals(∉(0±5)) == IntervalsUnion((Interval{:closed, :open}(-Inf, -5), Interval{:open, :closed}(5, Inf)))
    @test pred_to_intervals(>(10) ⩔ <(2) ⩔ ∈(4..5) ⩔ ==(6)) == IntervalsUnion((Interval{:open, :closed}(10, Inf), Interval{:closed, :open}(-Inf, 2), 4..5, 6..6))
    @test pred_to_intervals(>(10) ⩓ ∉(7..13) ⩓ <(21) ⩓ ∉(7..13) ⩓ >(5)) == IntervalsUnion(Interval{:open, :open}(13., 21.))
    @test pred_to_intervals(@optic(_^3 > 27)) == Interval{:open, :closed}(3, Inf)
    @test pred_to_intervals(@optic(_^3 >= 27)) == 3..Inf
    @test pred_to_intervals(@optic(-_^3 <= -27)) == 3..Inf
    @test pred_to_intervals(@optic(-_^3 ∈ -1..8)) == -2..1
    @test pred_to_intervals(@optic(_^3 ∈ IntervalsUnion((-1..8, 2..27)))) == -1..3
    @test pred_to_intervals(!@optic(_^3 > 27)) == IntervalsUnion(-Inf..3)
end


import CompatHelperLocal as CHL
CHL.@check()
