# Submesoscale eddies in the Eady model
# The parameters in the model roughly follow Taylor, 2016, GRL
using Oceananigans, Printf, JLD2, OceanBioME
using Oceananigans.Units

# Set the domain size and grid spacing
const Lx=1000;
const Ly=1000;
const Lz=140;
const Nx=128;
const Ny=128;
const Nz=16;

# Set the duration of the simulation 
duration = 10days;

# Construct a grid with uniform grid spacing
grid = RectilinearGrid(GPU(); size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz))

# Set the Coriolis parameter
const f = 1e-4

# Specify parameters that are used to construct the background state
background_state_parameters = (; M2 = 3e-8, # s⁻¹, geostrophic shear
                                 f,      # s⁻¹, Coriolis parameter
                                 Lz,
                                 N = 1e-4)            # s⁻¹, buoyancy frequency

# Here, B is the background buoyancy field and V is the corresponding thermal wind
@inline V(x, y, z, t, p) = + p.M2/p.f * (z - p.Lz/2)
@inline B(x, y, z, t, p) = p.M2 * x + p.N^2 * (z - p.Lz/2)

V_field = BackgroundField(V, parameters = background_state_parameters)
B_field = BackgroundField(B, parameters = background_state_parameters)

# Specify the horizontal and vertical viscosity/diffusivity
κ₂z = 1e-4 # [m² s⁻¹] Vertical vertical viscosity and diffusivity
κ₂h = 1e-2 # [m² s⁻¹] Horizontal viscosity and diffusivity

vertical_diffusivity = VerticalScalarDiffusivity(ν=κ₂z, κ=κ₂z)
horizontal_diffusivity = HorizontalScalarDiffusivity(ν=κ₂h, κ=κ₂h)

# Setup BGC

# nitrate - restoring to a (made up) climatology otherwise we depleat Nutrients
NO₃_forcing = Relaxation(rate = 1/10days, target = 4.0)

# alkalinity - restoring to a (made up) climatology otherwise we depleat it
Alk_forcing = Relaxation(rate = 1/day, target = 2409.0)

biogeochemistry = LOBSTER(; grid, 
                            carbonates = true, 
                            open_bottom = true,
                            surface_phytosynthetically_active_radiation = (x, y, t) -> 100)

# Fix temperature
DIC_bcs = FieldBoundaryConditions(top = GasExchange(; gas = :CO₂, temperature = (args...) -> 12, salinity = (args...) -> 35))

# Model instantiation
model = NonhydrostaticModel(;
                   grid,
                   biogeochemistry,
                   forcing = (NO₃ = NO₃_forcing, Alk = Alk_forcing),
                   boundary_conditions = (DIC = DIC_bcs, ),
              advection = CenteredSecondOrder(),
            timestepper = :RungeKutta3,
               coriolis = FPlane(; f),
                tracers = :b,
               buoyancy = BuoyancyTracer(),
      background_fields = (b = B_field, v = V_field),
                closure = (vertical_diffusivity, horizontal_diffusivity))

# ## Initial conditions
# Start with a bit of random noise added to the background thermal wind
Ξ(z) = randn() * z/grid.Lz * (z/grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u=uᵢ, v=vᵢ)


@load "validation/LOBSTER/steady_state_PAR_100_NO3_4_Alk_2409_deep.jld2" steady_state

for k in 1:grid.Nz
    model.tracers.P[:, :, k] .= steady_state.P[k]
    model.tracers.Z[:, :, k] .= steady_state.Z[k]
    model.tracers.NO₃[:, :, k] .= steady_state.NO₃[k]
    model.tracers.NH₄[:, :, k] .= steady_state.NH₄[k]
    model.tracers.DIC[:, :, k] .= steady_state.DIC[k]
    model.tracers.Alk[:, :, k] .= steady_state.Alk[k]
    model.tracers.sPOM[:, :, k] .= steady_state.sPOM[k]
    model.tracers.bPOM[:, :, k] .= steady_state.bPOM[k]
end

# specify a maximum timestep size
max_Δt = 15minutes

simulation = Simulation(model, Δt = max_Δt, stop_time = duration)

# Adapt the time step while keeping the CFL number fixed
wizard = TimeStepWizard(cfl=0.85, max_change=1.1, max_Δt=max_Δt)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# Create a progress message 
start_time = time_ns()
progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(1e-9 * (time_ns() - start_time)),
                        prettytime(sim.Δt),
                        AdvectiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# Here, add some diagnostics to calculate and output

u, v, w = model.velocities # unpack velocity `Field`s

# Periodically write the velocity, vorticity, and divergence out to a file
simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.tracers, (; u, v, w));
                                                      schedule = TimeInterval(4hours),
                                                      filename = "eady_turbulence_bgc",
                                                      overwrite_existing = true)
nothing # hide

# negativity protections 

using KernelAbstractions
using KernelAbstractions.Extras: @unroll
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device

@kernel function _remove_NaN_tendencies!(fields)
    i, j, k = @index(Global, NTuple)
    for field in fields
        if @inbounds isnan(field[i, j, k])
            field[i, j, k] = 0.0
        end
    end
end

@inline function remove_NaN_tendencies!(model)
    workgroup, worksize = work_layout(model.grid, :xyz)
    remove_NaN_tendencies_kernel! = _remove_NaN_tendencies!(device(model.grid.architecture), workgroup, worksize)
    event = remove_NaN_tendencies_kernel!(values(model.timestepper.Gⁿ))
    wait(event)
end

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())
simulation.callbacks[:nan_tendencies] = Callback(remove_NaN_tendencies!; callsite = TendencyCallsite())

@inline function zero_negative_tracers!(model)
    @unroll for tracer in values(model.tracers)
        parent(tracer) .= max.(0.0, parent(tracer))
    end
end

simulation.callbacks[:abort_zeros] = Callback(zero_negative_tracers!; callsite = UpdateStateCallsite())

run!(simulation)
#=
#####
##### Didn't actually run consecutivly
#####

n = 25
x₀ = [repeat([-200.0, -100.0, 0.0, 100.0, 200.0], 1, 5)...] .+ Lx / 2
y₀ = [repeat([-200.0, -100.0, 0.0, 100.0, 200.0], 1, 5)'...] .+ Ly / 2
z₀ = zeros(Float64, n)

function reset_location!(particles, model, Δt)
  particles.properties.z .= 0.0
end

particles = SLatissima.setup(; n, x₀, y₀, z₀, 
                               A₀ = 5.0, N₀ = 0.01, C₀ = 0.18, 
                               latitude = 57.5,
                               scalefactor = 10.0 ^ 8, 
                               optional_tracers = (:NH₄, :DIC, :bPON, :bPOC, :DON, :DOC),
                               custom_dynamics = reset_location!,
                               T = (args...) -> 12.0, S = (args...) -> 35.0, urel = 0.2)

biogeochemistry2 = LOBSTER(; grid, 
                             carbonates = true, 
                             open_bottom = true,
                             variable_redfield = true,
                             surface_phytosynthetically_active_radiation = (x, y, t) -> 100)

model2 = NonhydrostaticModel(; grid,
                               biogeochemistry = biogeochemistry2,
                               forcing = (NO₃ = NO₃_forcing, Alk = Alk_forcing),
                               boundary_conditions = (DIC = DIC_bcs, ),
                               advection = CenteredSecondOrder(),
                               timestepper = :RungeKutta3,
                               coriolis = coriolis,
                               tracers = :b,
                               buoyancy = BuoyancyTracer(),
                               background_fields = (b = B_field, v = V_field),
                               closure = (vertical_diffusivity, horizontal_diffusivity),
                               particles)

model2.clock.time = 50days

@load "eady_model.jld2" velocities tracers

model2.velocities.u[1:Nx, 1:Ny, 1:Nz] = velocities.u[1:Nx, 1:Ny, 1:Nz];
model2.velocities.v[1:Nx, 1:Ny, 1:Nz] = velocities.v[1:Nx, 1:Ny, 1:Nz];
model2.velocities.w[1:Nx, 1:Ny, 1:Nz] = velocities.w[1:Nx, 1:Ny, 1:Nz];

model2.tracers.b[1:Nx, 1:Ny, 1:Nz] = tracers.b[1:Nx, 1:Ny, 1:Nz];
model2.tracers.P[1:Nx, 1:Ny, 1:Nz] = tracers.P[1:Nx, 1:Ny, 1:Nz];
model2.tracers.Z[1:Nx, 1:Ny, 1:Nz] = tracers.Z[1:Nx, 1:Ny, 1:Nz];
model2.tracers.NO₃[1:Nx, 1:Ny, 1:Nz] = tracers.NO₃[1:Nx, 1:Ny, 1:Nz];
model2.tracers.NH₄[1:Nx, 1:Ny, 1:Nz] = tracers.NH₄[1:Nx, 1:Ny, 1:Nz];
model2.tracers.sPON[1:Nx, 1:Ny, 1:Nz] = tracers.sPOM[1:Nx, 1:Ny, 1:Nz];
model2.tracers.sPOC[1:Nx, 1:Ny, 1:Nz] = tracers.sPOM[1:Nx, 1:Ny, 1:Nz] .* biogeochemistry2.organic_redfield;
model2.tracers.bPON[1:Nx, 1:Ny, 1:Nz] = tracers.bPOM[1:Nx, 1:Ny, 1:Nz];
model2.tracers.bPOC[1:Nx, 1:Ny, 1:Nz] = tracers.bPOM[1:Nx, 1:Ny, 1:Nz] .* biogeochemistry2.organic_redfield;
model2.tracers.DON[1:Nx, 1:Ny, 1:Nz] = tracers.DOM[1:Nx, 1:Ny, 1:Nz];
model2.tracers.DOC[1:Nx, 1:Ny, 1:Nz] = tracers.DOM[1:Nx, 1:Ny, 1:Nz] .* biogeochemistry2.organic_redfield;
model2.tracers.DIC[1:Nx, 1:Ny, 1:Nz] = tracers.DIC[1:Nx, 1:Ny, 1:Nz];
model2.tracers.Alk[1:Nx, 1:Ny, 1:Nz] = tracers.Alk[1:Nx, 1:Ny, 1:Nz];

simulation2 = Simulation(model2, Δt = max_Δt, stop_time = 70days)

# Adapt the time step while keeping the CFL number fixed
wizard = TimeStepWizard(cfl=0.8, max_change=1.5)
simulation2.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# Create a progress message 
start_time = time_ns()

simulation2.callbacks[:progress] = Callback(progress, IterationInterval(10))

# Here, add some diagnostics to calculate and output

u, v, w = model2.velocities # unpack velocity `Field`s

# calculate the vertical vorticity [s⁻¹]
ζ = Field(∂x(v) - ∂y(u))

# horizontal divergence [s⁻¹]
δ = Field(∂x(u) + ∂y(v))

# Periodically write the velocity, vorticity, and divergence out to a file
simulation2.output_writers[:fields] = JLD2OutputWriter(model2, merge(model2.tracers, (; u, v, w, ζ, δ));
                                                      schedule = TimeInterval(1hours),
                                                      filename = "eady_turbulence_bgc_with_particles_dense",
                                                      overwrite_existing = true)


scale_negative_tracers = ScaleNegativeTracers(; model = model2, tracers = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON))
simulation2.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())
simulation2.callbacks[:nan_tendencies] = Callback(remove_NaN_tendencies!; callsite = TendencyCallsite())

simulation2.callbacks[:abort_zeros] = Callback(zero_negative_tracers!; callsite = UpdateStateCallsite())

simulation2.output_writers[:particles] = JLD2OutputWriter(model2, (particles = model2.particles, );
                                                                   schedule = TimeInterval(1hours),
                                                                   filename = "eady_turbulence_bgc_with_particles_dense_particles",
                                                                   overwrite_existing = true)

#run!(simulation)
=#