## Code for running model on spatial grid
## Run simulation on 2 by 1 grid - investigate spatial distributions
using RCall

# Set up initial parameters for ecosystem
numSpecies=16
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = collect(2:17)

# Set probabilities
birth = 0.6
death = 0.6
move = 0.02
l = 1.0
s = 0.0
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, move, timestep, l, s]

grid = (2,1)
totalK = 1000
individuals=100

# Create ecosystem
sppl = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec)
abenv = MatrixAbioticEnv(numNiches, grid, totalK)
eco = Ecosystem(sppl, abenv, false)

# Run simulations 100 times
ab = run_sim_spatial(eco, param, 1000, 500, 10, 100, true)


## Look at species 1
# Convert the abundances to Integers for species 1
tot = convert(Array{Int64,3}, ab[:,1,:,:])
# Sum over both grid squares to get the total abundance
total = mapslices(sum, tot , 3)[:, :, 1]

# Calculate
expected_dist = zeros(Float64, (length(total), maximum(total)+1))
for i in 1:length(total)
  expected_dist[i, 1:(total[i]+1)] = repmat([1/(total[i]+1)], total[i]+1)
end
expected = mapslices(sum, expected_dist, 1)

# Get the abundances from grid square 1
grid1 = tot[:,:,1]
# Calculate the counts of each abundance in grid 1 and convert to float
actual = counts(grid1+1, maximum(grid1+1))
actual = convert(Array{Float64,1}, actual)
# Cut expected values to length of actual
expected = expected[1:length(actual)]


@rput expected
@rput actual
# Calculate Kullback-Leibler divergence
KL = kldivergence(actual, expected); @rput KL
R"plot(1:length(expected),expected, type='l',
        main = paste('Divergence =', round(KL, 2)), xlab='Abundance',
        ylab='Frequency', ylim=c(0, max(c(expected, actual))))
        abline(h=max(expected), col=1, cex=0.5, lty=3)
lines(actual, col=2)
legend('topright', legend=c('Expected', 'Observed'), col=1:2, pch='-');
abline(h=max(actual), col=2, cex=0.5, lty=3)"

# Count how many of each abundance there are when the total is 10
# Is the distribution even?
count_tot = grid1[total.==10]
@rput count_tot
R"hist(count_tot, breaks=c(-0.5:10.5), main=' ', xlab='Abundance')"
# Stop here - try for all data combined and diff values of m


## Perform same analysis for ALL SPECIES
tot = convert(Array{Int64,4}, ab)
total = mapslices(sum, tot, 4)[:,:,:,1]

# Remove where total =0/1 because we know what the outcome should be
grid1 = tot[:,:,:,1]
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
move = 0.1
l = 1.0
s = 0.0
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, move, timestep, l, s]

grid = (5,5)
totalK = 1000
individuals=100

# Create ecosystem
sppl = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec)
abenv = MatrixAbioticEnv(numNiches, grid, totalK)
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
move = 0.05
l = 1.0
s = 0.0
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, move, timestep, l, s]

grid = (2,2)
totalK = 1000
individuals=100

# Create ecosystem
movement = GaussianMovement(move, 0.4, numSpecies)
sppl = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec, movement)
abenv = MatrixAbioticEnv(numNiches, grid, totalK)
eco = Ecosystem(sppl, abenv, false)

# Run simulations 100 times
ab = run_sim_spatial(eco, param, 1000, 500, 10, 10, false)

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
