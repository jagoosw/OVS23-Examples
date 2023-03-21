using JLD2, GLMakie, Oceananigans
using Oceananigans.Units

#####
##### Load tracers
#####

file = jldopen("eady_turbulence_bgc_with_particles_dense.jld2")

iterations = keys(file["timeseries/t"])

grid = file["serialized/grid"]

N, P, OC, DIC = ntuple(n -> ones(grid.Nx, grid.Ny, grid.Nz, length(iterations)) .* NaN, 4)
times = ones(length(iterations)) .* NaN

for (idx, it) in enumerate(iterations)
    N[:, :, :, idx] = file["timeseries/NO₃/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz] .+ file["timeseries/NH₄/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz]
    P[:, :, :, idx] = file["timeseries/P/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz] .* 6.56
    OC[:, :, :, idx] = file["timeseries/sPOC/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz] .+ file["timeseries/bPOC/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz] .+ file["timeseries/DOC/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz]
    DIC[:, :, :, idx] = file["timeseries/DIC/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz]

    times[idx] = file["timeseries/t/$it"]
end

close(file)

#####
##### Load particles
#####

file = jldopen("eady_turbulence_bgc_with_particles_dense_particles.jld2")

iterations = keys(file["timeseries/t"])

x, y, z, A, N_kelp, C = ntuple(n -> ones(25, length(iterations)) .* NaN, 6)

for (idx, it) in enumerate(iterations)
    particles = file["timeseries/particles/$it"]

    x[:, idx] = particles.x
    y[:, idx] = particles.y
    z[:, idx] = particles.z
    A[:, idx] = particles.A
    N_kelp[:, idx] = particles.N
    C[:, idx] = particles.C
end

close(file)

#####
##### Plot
#####

fig = Figure(resolution = (2665, 1344))

n = Observable(1)

N_plt = @lift N[1:16, :, :, $n]
P_plt = @lift P[17:32, :, :, $n]
OC_plt = @lift OC[33:48, :, :, $n]
DIC_plt = @lift DIC[49:end, :, :, $n]

x_plt = @lift x[:, $n]
y_plt = @lift y[:, $n]
z_plt = @lift z[:, $n]
A_plt = @lift A[:, $n]

xs = xnodes(Center, grid)[1:grid.Nx]
ys = ynodes(Center, grid)[1:grid.Ny]
zs = znodes(Center, grid)[1:grid.Nz]

lims = [(minimum(T), maximum(T)) for T in (N, P, OC, DIC)]

Aₘᵢ, Aₘₐ = minimum(A), maximum(A)

vm1 = contour(fig[1, 1], xs[1:17], ys, zs, N_plt, levels = 24, colormap = Reverse(:bamako))
vm2 = contour!(fig[1, 1], xs[17:33], ys, zs, P_plt, levels = 24, colormap = Reverse(:batlow))
vm3 = contour!(fig[1, 1], xs[33:49], ys, zs, OC_plt, levels = 24, colormap = :lajolla)
vm4 = contour!(fig[1, 1], xs[49:end], ys, zs, DIC_plt, levels = 24, colormap = Reverse(:devon))

sc = scatter!(x_plt, y_plt, z_plt, color = A_plt, colormap=:grayC)

txt = text!(
    [Point3f(xs[16 + (i - 1) * 16], ys[end], -100) for i in 1:4],
    text = ["Nutrients", "Phytoplankton", "Organic carbon", "Inorganic carbon"],
    rotation = [π for i in 1:4],
    align = (:left, :top),
    fontsize = 25,
    markerspace = :data
)

#=
Colorbar(fig[2, 1], vm1.plot, vertical = false, label = "Nutrient (NO₃ + NH₄) concentration (mmol N / m³)")
Colorbar(fig[3, 1], vm2, vertical = false, label = "Phytoplankton concentration (mmol C / m³)")
Colorbar(fig[4, 1], vm3, vertical = false, label = "Organic carbon concentration (mmol C / m³)")
Colorbar(fig[5, 1], vm4, vertical = false, label = "Inorganic carbon concentration (mmol C / m³)")
=#
n[] = 45

rotate_cam!(vm1.axis.scene, (0.2, 0, 0))
zoom!(vm1.axis.scene, 0.7 ^ 2)

sc1 = vm1.axis.scene.plots[1]
sc1.padding = 0.01
sc1.names.axisnames = ("x (m)", "y (m)", "z (m)")
sc1.names.fontsize = (20, 20, 20)
sc1.ticks.fontsize = (10, 10, 10)
save("ovs_3d.png", fig)

#=
record(fig, "ovs_3d.mp4", 1:length(iterations), framerate = 5) do i
    n[] = i
    msg = string("Plotting frame ", i, " of ", length(iterations))
    print(msg * " \r")
    #vm.axis.title = "t=$(prettytime(iterations[i]))"
end
=#


#####
##### Area plot inset (manually need to make background transparent)
#####

fig = Figure(resolution = (550, 420))
ax = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Frond area (dm²)", limits = (0, 19, 4, 6))
l = [lines!(ax, times[10:end]./day .- times[10]./day, A[idx, 10:end], color = A[idx, 10:end] .* 0.5 .* (C[idx, 10:end] .+ 0.2), colorrange = (0.5, 1.75), colormap = :lajolla) for idx in 1:25]
Colorbar(fig[0, 1], l[1], label = "Carbon storage (gC / frond)", vertical = false)

save("ovs_area.png", fig)

#####
##### colorbar inset
#####
fig = Figure(resolution = (550, 420))
Colorbar(fig[1, 1], vm1.plot, vertical = false, label = "Nutrient (NO₃ + NH₄) concentration (mmol N / m³)")
Colorbar(fig[2, 1], vm2, vertical = false, label = "Phytoplankton concentration (mmol C / m³)")
Colorbar(fig[3, 1], vm3, vertical = false, label = "Organic carbon (DOC + sPOC + bPOC) (mmol C / m³)")
Colorbar(fig[4, 1], vm4, vertical = false, label = "Inorganic carbon (mmol C / m³)")

save("ovs_colorbars.png", fig)

