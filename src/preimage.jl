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

function preimage(f, image::DistributionsExtra.IntervalsUnion)
	mapreduce(int -> preimage(f, int), ∪, image.ints)
end
