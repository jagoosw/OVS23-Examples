using OceanBioME, Oceananigans, Printf, JLD2
using OceanBioME.SLatissimaModel: SLatissima
using Oceananigans.Units

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / 365days))*(1 / (1 + 0.2 * exp(-((mod(t, 365days) - 200days) / 50days) ^ 2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, 365days) * (1 / (1 +exp(-(t - 100days) / (5days)))) * (1 / (1 + exp((t - 330days) / (25days))))

@inline MLD(t) = - (10 + 340 * (1 - fmld1(365days-eps(365days)) * exp(-mod(t, 365days) / 25days) - fmld1(mod(t, 365days))))

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t))/10)) / 2 + 1e-4

@inline t_function(x, y, z, t) = 2.4 * cos(t * 2π / year + 50day) + 10

# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
Lx, Ly = 20, 20
grid = RectilinearGrid(size=(1, 1, 50), extent=(Lx, Ly, 200)) 

# ## Kelp Particle setup
@info "Setting up kelp particles"
n = 5 # number of kelp fronds
z₀ = [-21:5:-1;]*1.0 # depth of kelp fronds

particles = SLatissima(; x = ones(n) * Lx / 2, y = ones(n) * Ly / 2, z = z₀, 
                         A = ones(n) * 10.0, N = ones(n) * 0.014, C =  ones(n) * 0.4, 
                         latitude = 57.5,
                         scalefactor = 500.0, 
                         pescribed_temperature = t_function)

# Specify the boundary conditions for DIC and O₂ based on the air-sea CO₂ and O₂ flux
CO₂_flux = GasExchange(; gas = :CO₂, temperature = t_function, salinity = (args...) -> 35)

N_fixation(x, y, t) = - 10 / year
N_upwelling = Relaxation(rate = 1/10days, mask = (x, y, z) -> ifelse(z < -150, 1.0, 0.0), target = 1.0)

model = NonhydrostaticModel(; grid,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ), 
                              biogeochemistry = LOBSTER(; grid,
                                                          surface_phytosynthetically_active_radiation = PAR⁰,
                                                          carbonates = true,
                                                          variable_redfield = true,
                                                          particles,
                                                          advection_schemes = (sPOM = WENO(grid), bPOM = WENO(grid))),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),
                                                     NO₃ = FieldBoundaryConditions(top = FluxBoundaryCondition(N_fixation))),
                              forcing = (NO₃ = N_upwelling, ),
                              advection = nothing)

# Running first time round for 10 years with `variable_redfield` and `particles` commented out and the following defaults from t =0
#set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0)

# Then uncommenting the below and adding back in variable redfield and particles from below time to `t = 3years`
# Also need to uncommnet particles writer and change negative scaling from `OM` to `ON` and change output name
model.clock.time = 1year - 30days

file = jldopen("column.jld2")
final_it = iterations = keys(file["timeseries/t"])[end-30]

model.tracers.P[1, 1, 1:50] = file["timeseries/P/$final_it"][1, 1, 1:50];
model.tracers.Z[1, 1, 1:50] = file["timeseries/Z/$final_it"][1, 1, 1:50];
model.tracers.NO₃[1, 1, 1:50] = file["timeseries/NO₃/$final_it"][1, 1, 1:50];
model.tracers.NH₄[1, 1, 1:50] = file["timeseries/NH₄/$final_it"][1, 1, 1:50];
model.tracers.sPON[1, 1, 1:50] = file["timeseries/sPOM/$final_it"][1, 1, 1:50];
model.tracers.sPOC[1, 1, 1:50] = file["timeseries/sPOM/$final_it"][1, 1, 1:50] .* model.biogeochemistry.organic_redfield;
model.tracers.bPON[1, 1, 1:50] = file["timeseries/bPOM/$final_it"][1, 1, 1:50];
model.tracers.bPOC[1, 1, 1:50] = file["timeseries/bPOM/$final_it"][1, 1, 1:50] .* model.biogeochemistry.organic_redfield;
model.tracers.DON[1, 1, 1:50] = file["timeseries/DOM/$final_it"][1, 1, 1:50];
model.tracers.DOC[1, 1, 1:50] = file["timeseries/DOM/$final_it"][1, 1, 1:50] .* model.biogeochemistry.organic_redfield;
model.tracers.DIC[1, 1, 1:50] = file["timeseries/DIC/$final_it"][1, 1, 1:50];
model.tracers.Alk[1, 1, 1:50] = file["timeseries/Alk/$final_it"][1, 1, 1:50];

close(file)

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

simulation = Simulation(model, Δt = 10minutes, stop_time=3years) 

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "column_particles"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)
simulation.output_writers[:particles] = JLD2OutputWriter(model, (; particles), filename = "$(filename)_particles.jld2", schedule = TimeInterval(1day), overwrite_existing=true)

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())

wizard = TimeStepWizard(cfl = 0.15, diffusive_cfl = 0.15, max_change = 2.0, min_change = 0.5, cell_diffusion_timescale = column_diffusion_timescale, cell_advection_timescale = column_advection_timescale)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# ## Run!
# Finally we run the simulation
run!(simulation)