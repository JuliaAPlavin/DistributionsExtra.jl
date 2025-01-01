# DistributionsExtra.jl

`Distributions.jl`-related functions that are too niche, experimental, or could not be included into `Distributions.jl` proper for other reasons.

Current content:
- `ℙ(pred, dist)`. Compute the probability of a predicate being true under a given distribution. In mathematical terms, it computes `ℙ(pred(x))` given that `x ~ dist`.
    ```julia
    julia> ℙ(>(0), Normal(0, 1))
    0.5

    julia> ℙ((@o abs(_) > 2), Normal(0, 1))
    0.04550026309183032
    ```
    See more details in the `ℙ` docstring.

- Support for intervals in distributions constructors, whenever it makes sense: `Uniform(1..5)`, `truncated(dist, 0..Inf)`, ... .
