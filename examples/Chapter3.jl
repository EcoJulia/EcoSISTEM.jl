using Diversity
using MyUnitful
using ClimatePref
using Simulation
#using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using Plots
pyplot()
using AxisArrays
function create_eco(numSpecies::Int64, grid::Tuple{Int64, Int64}, energy::Float64, req::Float64, individuals::Int64)
    # Set up initial parameters for ecosystem

    # Set up how much energy each species consumes
    energy_vec = SolarRequirement(fill(req .* kJ, numSpecies))
    #energy_vec = SimpleRequirement([2.0])
    # Set probabilities
    birth = 0.6/year
    death = 0.6/year
    l = 1.0
    s = 0.0
    boost = 1.0

    # Collect model parameters together (in this order!!)
    param = EqualPop(birth, death, l, s, boost)

    area = 4.0km^2
    #dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
    #dir = "/home/claireh/Documents/RawWorldclim"
    #srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
    #srad.array = srad.array[-10° .. 60°, 35° .. 80°]
    #meansrad = mean(srad.array[.!isnan.(srad.array)])
    #totalK = uconvert(kJ, meansrad * month * (area/(grid[1]*grid[2])))/5
    totalK = energy * kJ

    # Create ecosystem
    kernel = GaussianKernel(0.1km, numSpecies, 10e-10)
    movement = BirthOnlyMovement(kernel, Torus())

    opts = fill(274.0K, numSpecies)
    vars = fill(0.5K, numSpecies)
    traits = GaussTrait(opts, vars)
    native = fill(true, numSpecies)
    abun = rand(Multinomial(individuals, numSpecies))
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    abenv = simplehabitatAE(274.0K, grid, totalK, area)
    rel = Gauss{typeof(1.0K)}()
    eco = Ecosystem(sppl,abenv,rel)
end
function runsim(eco::Ecosystem, times::Unitful.Time, reps::Int64)
    burnin = 1month; interval = 1month; timestep = 1month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, reps)
    totalK = sum(eco.abenv.budget.matrix)
    grid = size(eco.abenv.habitat.matrix)
    for j in 1:reps
        repopulate!(eco)
        reenergise!(eco, totalK, grid)
        thisstore = view(abun, :, :, :, j)
        simulate_record!(thisstore, eco, times, interval, timestep)
    end
    return abun
end
numSpecies = 8; grd = (1,1); totalK = 1000000.0; req= 20.0; individuals=10000
eco = create_eco(numSpecies, grd, totalK, req, individuals)
abun = runsim(eco, 2000months, 100)
for n in 1:numSpecies
    if n == 1
        display(plot(abun[1,1,:,:], ylabel = "Abundance", xlabel = "Months", label="", grid = false, color = :1, linealpha = 0.1))
    else
        display(plot!(abun[n,1,:,:], label="", color = :($n), linealpha = 0.1))
    end
end
png("plots/8species_1square.png")

numSpecies = 8; grd = (4,4); totalK = 1000000.0; req= 20.0; individuals=10000
eco = create_eco(numSpecies, grd, totalK, req, individuals)
abun = runsim(eco, 2000months, 100)
sum_abun = mapslices(sum, abun, dims = 2)
for n in 1:8
    if n == 1
        display(plot(sum_abun[1,1,:,:], ylabel = "Abundance", xlabel = "Months", label="", grid = false, color = :1, linealpha = 0.1))
    else
        display(plot!(sum_abun[n,1,:,:], label="", color = :($n), linealpha = 0.1))
    end
end
png("plots/8species_16square.png")

function calc_extinction(abun::Array{Int64, 4})
    sum_abun = mapslices(sum, abun, dims = 2)[:,1,:,:]
    mat = mapslices(r-> (r[end] .== 0), sum_abun, dims = 2)[:,1,:]
    percent = mapslices(sum, mat, dims = 2)/100
    return percent
end
mean_ext = [1.0]
for i in 2.0:2:200
    numSpecies = 1; grd = (1,1); totalK = i; req= 2.0; individuals= Int64(i/req)
    eco = create_eco(numSpecies, grd, totalK, req, individuals)
    abun = runsim(eco, 2000months, 100)
    push!(mean_ext, mean(calc_extinction(abun)))
end
plot(mean_ext, ylabel = "Probability of extinction", xlabel = "Energy", label="", grid = false, color = :1)
png("plots/Energy_vs_extinction.png")

abun = runsim(2000months, 100)
plot(abun[1,1,:,:], ylabel = "Abundance", xlabel = "Months", label="", grid = false, color = :black, linealpha = 0.1)
png("plots/1species_1square.png")

numSpecies = 8; grd = (1,1); totalK = 250000.0kJ/km^2; area = 4.0km^2; req= 20.0kJ; individuals=10000
eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
abun = runsim(eco, 2000months, 100)
plot_abun(abun, 10)

for i in [1, 2, 4, 8, 16]
    numSpecies = i; grd = (1,1); totalK = 250000.0kJ/km^2; area = 4.0km^2; req= 20.0kJ; individuals=10000
    eco = create_eco(numSpecies, grd, area, totalK, req, individuals)
    abun = runsim(eco, 2000months, 10)
    if i == 1
        display(plot([1/i], [mean(abun[:,:,1500:end,:])], seriestype = :scatter, label = "", ylim = (0, 5*10^4), xlim = (0,1.1), grid = false, ylab = "Average abundance at eqm", xlab= "1 / number of species"))
    else
        display(plot!([1/i], [mean(abun[:,:,1500:end,:])], seriestype = :scatter, label = ""))
    end
end


using Profile
using ProfileView
eco = create_eco(numSpecies, grd, totalK, req, individuals)
runsim(eco, 1month, 1)  # run once to trigger compilation
Profile.clear()
@profile runsim(eco, 2000months, 100)
ProfileView.view()

using Diversity
using MyUnitful
using ClimatePref
using Simulation
#using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using Plots
pyplot()
using AxisArrays
# Set up initial parameters for ecosystem
numSpecies = 50

# Set up how much energy each species consumes
energy_vec = SolarRequirement(sample(1000.0:4000, numSpecies) .* kJ)

# Set probabilities
birth = 0.6/year
death = 0.6/year
l = 1.0
s = 0.8
boost = 1.0
timestep = 1month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s, boost)

grid = (20, 20)
area = 4.0km^2
dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
#dir = "/home/claireh/Documents/RawWorldclim"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad.array = srad.array[-10° .. 60°, 35° .. 80°]
meansrad = mean(srad.array[.!isnan.(srad.array)])
totalK = uconvert(kJ, meansrad * month * (area/(grid[1]*grid[2])))
individuals = 10000

# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-10)
movement = BirthOnlyMovement(kernel, Cylinder())

opts = rand(Normal(274.0, 10.0), numSpecies) * K
vars = rand(Uniform(0, 25/9), numSpecies) * K
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param, native)
abenv = simplehabitatAE(274.0K, grid, totalK, area)
rel = Gauss{typeof(1.0K)}()
eco = Ecosystem(sppl,abenv,rel)
plot(eco)
hab2 = HabitatCollection2(abenv.habitat, abenv2.habitat)
plot(hab2)
plot(abenv.habitat)

burnin = 1month; interval = 1month; reps = 1; times = 10years
lensim = length(0month:interval:times)
sum_abun = [sum(eco.abundances.matrix)]
for i in 1:20
    simulate!(eco, burnin, timestep)
    push!(sum_abun, sum(eco.abundances.matrix))
end
plot(sum_abun, ylabel = "Abundance", xlabel = "Months", label="", grid = false)
#png("plots/burnin.png")
plot(mapslices(sum, eco.abundances.matrix, dims = 2), seriestype = :bar,
    grid =false, label="")
plot(eco)
