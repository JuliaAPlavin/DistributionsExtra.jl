using TestItems
using TestItemRunner
@run_package_tests


@testitem "to intervals" begin
    using DistributionsExtra: pred_to_intervals, IntervalUnion

    @test pred_to_intervals(∈(0±5.)) == -5..5
    @test pred_to_intervals(∉(0±5.)) == IntervalUnion((Interval{:open, :open}(-Inf, -5), Interval{:open, :open}(5, Inf)))
    @test pred_to_intervals(>(10) ⩔ <(2) ⩔ ∈(4..5) ⩔ ==(6)) == IntervalUnion((Interval{:open, :closed}(10, Inf), Interval{:closed, :open}(-Inf, 2), 4..5, 6..6))
    @test pred_to_intervals(>(10) ⩓ ∉(7..13) ⩓ <(21) ⩓ ∉(7..13) ⩓ >(5)) == IntervalUnion(Interval{:open, :open}(13., 21.))
    @test pred_to_intervals((>(10.) ⩔ <(2) ⩔ ∈(4..5.0) ⩔ ==(6.0)) ⩓ !=(0)) == IntervalUnion((Interval{:open, :open}(10, Inf), Interval{:open, :open}(-Inf, 0), Interval{:open, :open}(0, 2), 4..5, 6..6))
    @test pred_to_intervals(@o(_^3 > 27)) == Interval{:open, :closed}(3, Inf)
    @test pred_to_intervals(@o(_^3 >= 27)) == 3..Inf
    @test pred_to_intervals(@o(-_^3 <= -27)) == 3..Inf
    @test pred_to_intervals(@o(-_^3 ∈ -1..8)) == -2..1
    @test pred_to_intervals(@o(_^3 + 1 >= 28)) == 3..Inf
    @test pred_to_intervals(@o(_^3 ∈ IntervalUnion((-1..8, 2..27)))) == -1..3
    @test pred_to_intervals(!@o(_^3 > 27)) == IntervalUnion(Interval{:open,:closed}(-Inf, 3))
    @test pred_to_intervals(@o(abs(_) > 5)) == IntervalUnion((Interval{:open, :closed}(5, Inf), Interval{:closed, :open}(-Inf, -5)))
    @test pred_to_intervals(@o(abs(_) < 5)) == Interval{:open, :open}(-5, 5)
    @test_broken pred_to_intervals(@o(_^2 < 25)) == Interval{:open, :open}(-5, 5)
    @test pred_to_intervals(@o(log10(_) ∈ -1..2)) == 0.1..100
    @test_broken pred_to_intervals(@o(10^_ ∈ -1..100)) == -Inf..2
end

@testitem "P" begin
    d = Poisson(1)
    @test ℙ(>=(0), d) == 1
    @test ℙ(<(0), d) == 0
    @test ℙ(<=(0), d) == ℙ(==(0), d) == 1/ℯ
    @test ℙ(==(2.5), d) == 0
    @test ℙ(>(0), d) == ℙ(>(0.5), d) == ℙ(>=(1), d) ≈ 1 - 1/ℯ
    @test ℙ(<(2), d) == ℙ(<(1.5), d) == ℙ(<=(1), d) == ℙ(∈(0:1), d) == ℙ(∈(0..1), d) ≈ 0.7357588823428847
    @test ℙ(∈(1), d) == ℙ(==(1), d) ≈ 1/ℯ
    @test ℙ(∈([0, 2]), d) == ℙ(==(0) ⩔ ==(2), d) ≈ 0.5518191617571635
    @test ℙ(∈(0..2.1), d) == ℙ(∈(0:2), d) == ℙ(>=(0) ⩓ <(2.5), d) ≈ 0.9196986029286058
    @test ℙ(>=(10), Poisson(1.45)) ≈ 3.053405450128097e-6

    d = censored(Normal(0, 1), 0..Inf)
    @test ℙ(==(0), d) == 0.5
    @test ℙ(!=(0), d) == 0.5
    @test ℙ(<(0), d) == 0
    @test ℙ(≤(0), d) == 0.5
    @test ℙ(>(0), d) == 0.5
    @test ℙ(≥(0), d) == 1
    @test ℙ(≥(5), d) ≈ 2.866515940330412e-7  rtol=∛eps()
    @test ℙ((@o ≈(_, 3, atol=1e-3)), d) ≈ pdf(d, 3)*2*1e-3  rtol=1e-4

    @testset for P in [ℙ, ℙᵋ]
        d = Normal(0, 1)
        @test P(>=(5), d) ≈ 2.866515940330412e-7  rtol=∛eps()
        @test P(∉(0±5), d) == P(@o(abs(_) > 5), d)
        @test P(∉(0±5), d) ≈ 5.73303121648882e-7  rtol=∛eps()
        @test P(@o(abs(_) > 2), d) + P(@o(abs(_) < 2), d) ≈ 1
    end
end

@testitem "_" begin
    import CompatHelperLocal as CHL
    CHL.@check()

    using Aqua
    Aqua.test_all(DistributionsExtra, ambiguities=false)
end
