using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase

do_plot = false

# Set simulation parameters
birth = [0.0/day; fill(1e-5/day, 3); 0.0/day]
death = [0.0/day; fill(1e-5/day, 3); 0.0/day]
virus_growth = 0.1/day
virus_decay = 1.0/day
beta = 1e-2/day
sigma = 1/14days
param = SIRGrowth{typeof(unit(beta))}(birth, death, virus_growth, virus_decay, beta, sigma)
param = transition(param)

# Read in population sizes for Scotland
scotpop = Array{Float64, 2}(readfile("test/examples/ScotlandDensity2011.tif", 0.0, 7e5, 5e5, 1.25e6))

# Set up simple gridded environment with initial population scotpop
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, area, NoControl(), scotpop)

# Set population to initially have no individuals
abun = fill(0, 5)

# Dispersal kernels for virus and disease classes
dispersal_dists = [1e-2km; fill(2.0km, 3); 1e-2km]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, 5), fill(0.1K, 5))
epilist = SIR(traits, abun, movement, param)
rel = Gauss{eltype(epienv.habitat)}()

# Create epi system with all information
# This will fill in initial susceptible population (scotpop) from epienv
epi = EpiSystem(epilist, epienv, rel)

# Add in initial infections randomly (samples weighted by population size)
samp = sample(1:525_000, weights(1.0 .* epi.abundances.matrix[2, :]), 100)
epi.abundances.matrix[1, samp] .= 100 # Virus pop
epi.abundances.matrix[3, samp] .= 10 # Infected pop

# Run simulation
abuns = zeros(Int64, 5, 525_000, 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

# leave pkgs out of the `do_plot` part as otherwise we will come
#   across errors like  `UndefVarError: @layout not defined` when
#   `do_plot=false`
using Plots
plotlyjs()
if do_plot
    # Plot heatmap of spatial disease progression
    abuns = Float64.(abuns)
    infecteds = abuns[3, :, :]
    infecteds[.!epi.epienv.active[1:end], :] .= NaN
    display(heatmap(reshape(infecteds[:, 1], 700, 750)', layout = (@layout [a b; c d]), subplot = 1, title = "Day 1", clim = (0, 2600), background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, size = (1500, 1000), color = :fire_r))
    display(heatmap!(reshape(infecteds[:, 7], 700, 750)', subplot = 2, title = "Day 7", clim = (0, 2600), background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, color = :fire_r))
    display(heatmap!(reshape(infecteds[:, 14], 700, 750)', subplot = 3, title = "Day 14", clim = (0, 2600), background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, color = :fire_r))
    display(heatmap!(reshape(infecteds[:, 30], 700, 750)', subplot = 4, title = "Day 30", clim = (0, 2600), background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, color = :fire_r))

    # View summed SIR dynamics for whole area
    display(plot(mapslices(sum, abuns[2, :, :], dims = 1)[1, :], color = :Blue, label = "Susceptible"))
    display(plot!(mapslices(sum, abuns[3, :, :], dims = 1)[1, :], color = :Red, label = "Infected"))
    display(plot!(mapslices(sum, abuns[4, :, :], dims = 1)[1, :], color = :Black, label = "Recovered"))

    # Create animation of both
    abuns = Float64.(abuns)
    infecteds = abuns[3, :, :]
    infecteds[.!epi.epienv.active[1:end], :] .= NaN
    anim = @animate for i ∈ 1:365
        heatmap(reshape(infecteds[:, i], 700, 750)', layout = (@layout [a; b]), subplot = 1, title = "Day $i", clim = (0, 2500), grid = false, aspect_ratio = 1, size = (1200, 1500), color = :magma_r)
        plot!(mapslices(sum, abuns[2, :, 1:i], dims = 1)[1, :], color = :Blue, label = "Susceptible", subplot = 2, xlim = (0, 365), ylim = (0, 5.5e6), legend = :bottom, bottom_margin = 20* Plots.mm)
        plot!(mapslices(sum, abuns[3, :, 1:i], dims = 1)[1, :], color = :Red, label = "Infected", subplot = 2, ylim = (0, 5.5e6), xlim = (0, 365))
        plot!(mapslices(sum, abuns[4, :, 1:i], dims = 1)[1, :], color = :Black, label = "Recovered", subplot = 2, ylim = (0, 5.5e6), xlim = (0, 365))
    end every 7
    gif(anim, "test/examples/ScotlandSIRSim.gif", fps = 15)
end
