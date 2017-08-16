## Code for running model on spatial grid

# Include simulation functions
using Diversity
using Simulation
using RCall
using Distributions
using StatsBase
using ProgressMeter
using DataStructures
using Unitful
using Unitful.DefaultSymbols

## Run simulation on 2 by 1 grid - investigate spatial distributions

# Set up initial parameters for ecosystem
numSpecies=8
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(collect(2:9))

# Set probabilities
birth = 0.6/month
death = 0.6/month
l = 1.0
s = 0.0
boost = 1.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s, boost)

grid = (2, 1)
area = 10.0km^2
totalK = 1000.0
individuals=100

# Create ecosystem
kernel = GaussianKernel(0.2, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)
sppl = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec, movement, param)
abenv = simplenicheAE(numNiches, grid, totalK, area)
rel = TraitRelationship(SimpleNiche)
eco = Ecosystem(sppl, abenv, rel)

plotdiv(norm_sub_alpha, eco, collect(0.0:10))

times = 10year; burnin = 1year; interval = 3month; reps = 10
# Run simulation grid
lensim = length(0month:interval:times)
#times = 1000; burnin = 500; interval = 10
#lensim = length(0:interval:times)

# Run simulations 10 times
storage = generate_storage(eco, lensim, reps)

@showprogress 1 "Computing..." for j in 1:reps
  if (j != 1) repopulate!(eco) end
  thisstore = view(storage, :, :, :, j)
  simulate!(eco, burnin, interval, timestep)
  simulate_record!(thisstore, eco, times, interval, timestep)
end
plot_abun(storage, numSpecies, grid)

ab = run_sim_spatial(eco, param, times, burnin, interval, reps)

## Look at species 1
# Convert the abundances to Integers for species 1

comp = expected_counts(ab, 1, 1)

plot_divergence(comp)


# Count how many of each abundance there are when the total is 10
# Is the distribution even?
freq_hist(ab, 1, 10, 1)
# Stop here - try for all data combined and diff values of m


## Perform same analysis for ALL SPECIES
comp = expected_counts(ab, 1)
plot_divergence(comp)

# Count how many of each abundance there were with 10 in total
freq_hist(ab, 1, 10)


## Run simulation on 5 by 5 grid - investigate spatial distributions
using RCall

# Set up initial parameters for ecosystem
numSpecies=8
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = collect(2:9)

# Set probabilities
birth = 0.6
death = 0.6
l = 1.0
s = 0.0
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, timestep, l, s]

grid = (5,5)
totalK = 1000
individuals=100

# Create ecosystem
sppl = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec)
abenv = simplenicheAE(numNiches, grid, totalK)
eco = Ecosystem(sppl, abenv, false)

# Run simulations 100 times
ab1 = run_sim_spatial(eco, param, 500, 200, 10, 50, false)
ab2 = run_sim_spatial(eco, param, 500, 200, 10, 50, false)

run1 = ab1
run2 = ab2
## Perform same analysis for ALL SPECIES
ab = cat(3, run1, run2)
tot = convert(Array{Int64,4}, ab)


total = mapslices(sum, tot[:,:,[1:50 51:100],1], 4)[:,:,:,1]

# Remove where total =0/1 because we know what the outcome should be
grid1 = tot[:,:,1:50, 1]
grid1 = grid1[total.>0]
total = total[total.>0]

# Calculate the expected distribution if things were evenly distributed over space
expected_dist = zeros(Float64, (length(total), maximum(total)+1))
for i in 1:length(total)
  expected_dist[i, 1:(total[i]+1)] = repmat([1/(total[i]+1)], total[i]+1)
end
expected = mapslices(sum, expected_dist, 1)

# Find the actual counts of the species
actual = counts(grid1+1, maximum(grid1+1))
actual = convert(Array{Float64,1}, actual)
# Trim expected values so they are the same length as the actaul
expected = expected[1:length(actual)]


@rput expected
@rput actual
@rput numSpecies; @rput move
# Calculate Kullback-Leibler divergence
KL = kldivergence(actual, expected); @rput KL
R"par(mfrow=c(1,1));plot(1:length(expected),expected, type='l',
        main = paste('Divergence =', round(KL, 2)), xlab='Abundance',
        ylab='Frequency', ylim=c(0, max(c(expected, actual))))
        abline(h=max(expected), col=1, cex=0.5, lty=3)
lines(actual, col=2)
legend('topright', legend=c('Expected', 'Observed'), col=1:2, pch='-');
abline(h=max(actual), col=2, cex=0.5, lty=3)"

# Count how many of each abundance there were with 10 in total
count_tot = grid1[total.==10]
@rput count_tot
R"hist(count_tot, breaks=c(-0.5:10.5), main=' ', xlab='Abundance')"


using RCall

# Set up initial parameters for ecosystem
numSpecies=8
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = collect(2:9)

# Set probabilities
birth = 0.6
death = 0.6
l = 1.0
s = 0.0
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, timestep, l, s]

grid = (2,2)
totalK = 1000
individuals=100

# Create ecosystem
movement = GaussianMovement(0.4, numSpecies, 10e-5)
sppl = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec, movement)
abenv = simplenicheAE(numNiches, grid, totalK, 1)
eco = Ecosystem(sppl, abenv, false)

# Run simulations 100 times
ab = run_sim_spatial(eco, param, 1000, 500, 10, 50, false)

## Perform same analysis for ALL SPECIES
tot = convert(Array{Int64,4}, ab)
total = mapslices(sum, tot, 4)[:,:,:,1]

# Remove where total =0/1 because we know what the outcome should be
grid1 = mapslices(sum, tot[:,:,:,1:2], 4)[:,:,:,1]
grid1 = grid1[total.>0]
total = total[total.>0]

# Calculate the expected distribution if things were evenly distributed over space
expected_dist = zeros(Float64, (length(total), maximum(total)+1))
for i in 1:length(total)
  expected_dist[i, 1:(total[i]+1)] = repmat([1/(total[i]+1)], total[i]+1)
end
expected = mapslices(sum, expected_dist, 1)

# Find the actual counts of the species
actual = counts(grid1+1, maximum(grid1+1))
actual = convert(Array{Float64,1}, actual)
# Trim expected values so they are the same length as the actaul
expected = expected[1:length(actual)]


@rput expected
@rput actual
@rput numSpecies; @rput move
# Calculate Kullback-Leibler divergence
KL = kldivergence(actual, expected); @rput KL
R"par(mfrow=c(1,1));plot(1:length(expected),expected, type='l',
        main = paste('Divergence =', round(KL, 2)), xlab='Abundance',
        ylab='Frequency', ylim=c(0, max(c(expected, actual))))
        abline(h=max(expected), col=1, cex=0.5, lty=3)
lines(actual, col=2)
legend('topright', legend=c('Expected', 'Observed'), col=1:2, pch='-');
abline(h=max(actual), col=2, cex=0.5, lty=3)"

# Count how many of each abundance there were with 10 in total
count_tot = grid1[total.==10]
@rput count_tot
R"hist(count_tot, breaks=c(-0.5:10.5), main=' ', xlab='Abundance')"
