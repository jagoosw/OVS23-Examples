using JLD2, CairoMakie, Statistics
using Oceananigans.Units

@inline t_function(x, y, z, t) = 2.4 * cos(t * 2π / year + 50day) + 10

P, NO₃, DIC, Alk, sPOC, bPOC = ntuple(n -> zeros(50, 5 * 365 - 183), 6)

times = zeros(1825 - 183)

file = jldopen("column.jld2")
iterations = keys(file["timeseries/t"])[end-1095:end-30]

for (idx, it) in enumerate(iterations)
    P[:, idx] = file["timeseries/P/$it"][1, 1, 1:50]
    NO₃[:, idx] = file["timeseries/NO₃/$it"][1, 1, 1:50]
    DIC[:, idx] = file["timeseries/DIC/$it"][1, 1, 1:50]
    Alk[:, idx] = file["timeseries/Alk/$it"][1, 1, 1:50]
    sPOC[:, idx] = file["timeseries/sPOM/$it"][1, 1, 1:50] .* 6.56
    bPOC[:, idx] = file["timeseries/bPOM/$it"][1, 1, 1:50] .* 6.56
    times[idx] = file["timeseries/t/$it"]
end

close(file)

initial_offset = length(iterations)

file = jldopen("column_particles.jld2")
iterations = keys(file["timeseries/t"])[2:2*365 + 30 - 183]

for (idx, it) in enumerate(iterations)
    P[:, idx + initial_offset] = file["timeseries/P/$it"][1, 1, 1:50]
    NO₃[:, idx + initial_offset] = file["timeseries/NO₃/$it"][1, 1, 1:50]
    DIC[:, idx + initial_offset] = file["timeseries/DIC/$it"][1, 1, 1:50]
    Alk[:, idx + initial_offset] = file["timeseries/Alk/$it"][1, 1, 1:50]
    sPOC[:, idx + initial_offset] = file["timeseries/sPOC/$it"][1, 1, 1:50]
    bPOC[:, idx + initial_offset] = file["timeseries/bPOC/$it"][1, 1, 1:50]
    times[idx + initial_offset] = file["timeseries/t/$it"] + times[initial_offset] - file["timeseries/t/0"]
end

times .-= times[1]

close(file)

file = jldopen("column_particles_particles.jld2")
iterations = keys(file["timeseries/t"])[1:2*365 + 30 - 183]

A, N, C = ntuple(n -> zeros(5, 760 - 183), 3)
kelp_times = zeros(760 - 183)

for (idx, it) in enumerate(iterations)
    particles = file["timeseries/particles/$it"]
    A[:, idx] = particles.A
    N[:, idx] = particles.N
    C[:, idx] = particles.C
    kelp_times[idx] = file["timeseries/t/$it"]
end

close(file)

kelp_times .-= kelp_times[1]

air_sea_CO₂_flux = zeros(1825 - 183)
carbon_export = zeros(1825 - 183)
for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = CO₂_flux.condition.parameters(0.0, 0.0, t, DIC[50, i], Alk[50, i], t_function(1, 1, 0, t), 35)
    carbon_export[i] = (sPOC[end - 20, i] * 3.47e-5 + bPOC[end - 20, i] * 200/day) 
end

zs = [-198:4:-2;]

fig = Figure(resolution = (1600, 960))

axP = Axis(fig[1:3, 1:2], xlabel = "Time (years)", ylabel = "Depth (m)", title = "Phytoplankton Concentraiton (mmol N / m³)")

hmP = heatmap!(axP, times ./ year, zs[10:50], P[10:50, :]', colormap = Reverse(:batlow), interpolate=true)
Colorbar(fig[1:3, 3], hmP)

lines!(axP, [3 - 30/365, 3 - 30/365], [zs[9], 0], color=:black)

axN = Axis(fig[4:6, 1:2], xlabel = "Time (years)", ylabel = "Depth (m)", title = "Nutrient Concentraiton (mmol N / m³)")

hmN = heatmap!(axN, times ./ year, zs[10:50], NO₃[10:50, :]', colormap = Reverse(:batlow), interpolate=true)
Colorbar(fig[4:6, 3], hmN)

lines!(axN, [3 - 30/365, 3 - 30/365], [zs[9], 0], color=:black)

axS = Axis(fig[7:9, 1:4], xlabel = "Time (years)", ylabel = "Carbon Flux (kg CO₂ / m² / year)", limits = (0, times[end] /years, 1.1 * (12 + 16 * 2) * year /(1000 * 1000) * min(minimum(air_sea_CO₂_flux), minimum(-carbon_export)), 1.2 * (12 + 16 * 2) * year /(1000 * 1000) * max(maximum(air_sea_CO₂_flux), maximum(-carbon_export))))

lnAS = lines!(axS, times ./ year, air_sea_CO₂_flux .* (12 + 16 * 2) .* year /(1000 * 1000), label = "Air-sea exchange")
lnSI = lines!(axS, times ./ year, - carbon_export .* (12 + 16 * 2) .* year /(1000 * 1000), label = "Sinking export")

Legend(fig[7:9, 1:4], [lnAS, lnSI], ["Air-sea CO₂ exchange", "Sinking export"], halign = :left, valign = :bottom)

lines!(axS, [3 - 30/365, 3 - 30/365], [1, -2], color=:black)

axA = Axis(fig[1:2, 4], xlabel = "Time (years)", ylabel = "Frond area (dm² / frond)", title = "Kelp growth (equivilant to growing 500 frond / m² in the top 50m of water)")

lines!(axA, kelp_times ./ year, mean(A, dims=1)[1, :])

axC = Axis(fig[3:4, 4], xlabel = "Time (years)", ylabel = "Carbon stored (kg CO₂ / m²)")

lines!(axC, kelp_times ./ year, sum((C .+ 0.2) .* A .* 0.5 .* 100 * (12 + 16 * 2)/ (12 * 1000), dims = 1)[1, :])

axN = Axis(fig[5:6, 4], xlabel = "Time (years)", ylabel = "Nitrogen stored (mmol N / m²)")

lines!(axN, kelp_times ./ year, sum((N .+ 0.0146) .* A .* 0.5 .* 100 / 14, dims = 1)[1, :])
