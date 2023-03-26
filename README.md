# Ocean Visions Summit 2023

Here you will find the code used to generate the figures for my poster at Ocean Visions Summit 2023.

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
https://user-images.githubusercontent.com/26657828/227775767-61fbf724-6b71-4181-8304-0891f8cdd5c5.mp4

## Column
![ovs_column](https://user-images.githubusercontent.com/26657828/227775254-a510014d-91ad-4b33-9cc3-2278a9082d20.png)
