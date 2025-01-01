function preimage(f, image::IntervalUnion)
    mapreduce(int -> preimage(f, int), ∪, image.ints)
end

function preimage(f, image::Interval)
    f⁻¹ = inverse(f)
    eps = f⁻¹.(endpoints(image))
    if eps[2] >= eps[1]
        @set endpoints(image) = eps
    else
        @p begin
            image
            @set endpoints(__) = reverse(eps)
            @modify(reverse, closedendpoints(__))
        end
    end
end

preimage(::typeof(abs), image::Interval) =
    if 0 ∈ image
        OC = isrightclosed(image) ? :closed : :open
        Interval{OC, OC}(-rightendpoint(image), rightendpoint(image))
    else
        neg = @p image |>
            @modify(reverse ∘ (.-), endpoints(__)) |>
            @modify(reverse, closedendpoints(__))
        IntervalUnion((image, neg))
    end
