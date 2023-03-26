# Ocean Visions Summit 2023

Here you will find the code used to generate the figures for my poster at Ocean Visions Summit 2023 as well as a copy of the poster.

## How to run:
The required dependencies should be detailed in `Project.toml`, so you can instantiate the environment and then run the simulations:

```julia
julia> using Pkg
julia> Pkg.instantiate()
julia> include("eady/simulaiton.jl")
```

You can then run the plotting files:
```julia
julia> include("eady/animation.jl")
```

# Results
## Baroclinical instability with kelp-seeded buoys

## Column