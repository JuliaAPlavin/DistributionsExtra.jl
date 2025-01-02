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

- Support for intervals from `IntervalSets.jl` whenever it makes sense *(unfortunately, [the corresponding PR](https://github.com/JuliaStats/Distributions.jl/pull/1809) to `Distributions.jl` wasn't accepted)*:
    - `Uniform(1..5)`, `truncated(dist, 0..Inf)`, ...
    - convert to/from `Distributions.RealInterval` custom type
    - `uniform(X)` function that returns the uniform distribution on `X`; works for intervals, ranges, and vectors

- Finally, some specialized distribution types. They are intended to only be here temporarily, and upstreamed to `Distributions.jl` or `Manifolds.jl` at some point.
    - `PiecewiseUniform` ([Distributions.jl PR](https://github.com/JuliaStats/Distributions.jl/pull/1367) stuck for years)
    - Uniform on a sphere: `SphereUniformArea` for 3d xyz parametrization, and `SphereUniformLonLat` for 2d lon-lat parametrization (a useful distribution but the specific semantics is tricky, [[1]](https://github.com/JuliaManifolds/Manifolds.jl/issues/708))
    - Piecewise-uniform on a sphere: `SpherePiecewiseLatUniformArea` for 2d lon-lat parametrization
