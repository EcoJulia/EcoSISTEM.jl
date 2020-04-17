using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase

# Larger grid
birth = [0.00001/day; fill(1e-10/day, 3)]
death = [0.005/day; fill(1e-10/day, 3)]
beta = 1e-2/day
sigma = 0.1/14days
param = EpiGrowth{typeof(unit(beta))}(birth, death, beta, sigma)

scotpop = Array{Float64, 2}(readfile("/Users/claireh/Documents/COVID/ScotlandPopulation2011.tif", 0.0, 7e5, 5e5, 1.25e6))
active = Array{Bool, 2}(.!isnan.(scotpop))
popdiv = [sum(scotpop .== i) for i in unique(scotpop)]
un_pops = unique(scotpop)
for i in 1:length(un_pops)
    scotpop[scotpop .== un_pops[i]] ./= popdiv[i]
end
scotpop = round.(scotpop)

grid = (700, 750)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, active, NoControl())

abun = fill(0, 4)

dispersal_dists = [2.0km; fill(0.01km, 3)]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

traits = GaussTrait(fill(298.0K, 4), fill(0.1K, 4))
epilist = SIR(traits, abun, movement, param)

rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)
scotpop[isnan.(scotpop)] .= 0
epi.abundances.grid[2, :, :] .+= scotpop
samp = sample(1:525_000, weights(1.0 .* epi.abundances.matrix[2, :]), 100)
epi.abundances.matrix[1, samp] .= 100
epi.abundances.matrix[3, samp] .= 10

abuns = zeros(Int64, 4, 525_000, 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

using Plots
plotlyjs()
abuns = Float64.(abuns)
infecteds = abuns[3, :, :]
infecteds[.!active[1:end], :] .= NaN
display(heatmap(reshape(infecteds[:, 1], 700, 750)', layout = (@layout [a b; c d]), subplot = 1, title = "Day 1", clim = (0, 2000), background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, size = (1500, 1000), color = :deep))
display(heatmap!(reshape(infecteds[:, 7], 700, 750)', subplot = 2, title = "Day 7", clim = (0, 2000), background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, color = :deep))
display(heatmap!(reshape(infecteds[:, 14], 700, 750)', subplot = 3, title = "Day 14", clim = (0, 2000), background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, color = :deep))
display(heatmap!(reshape(infecteds[:, 30], 700, 750)', subplot = 4, title = "Day 30", clim = (0, 2000), background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, color = :deep))
Plots.pdf("/Users/claireh/Documents/COVID/ScotlandSim2.pdf")


display(plot(mapslices(sum, abuns[2, :, :], dims = 1)[1, :], color = :Blue, label = "Susceptible"))
display(plot!(mapslices(sum, abuns[3, :, :], dims = 1)[1, :], color = :Red, label = "Infected"))
display(plot!(mapslices(sum, abuns[4, :, :], dims = 1)[1, :], color = :Black, label = "Recovered"))
Plots.pdf("/Users/claireh/Documents/COVID/ScotlandSimSIR.pdf")


@animate 
