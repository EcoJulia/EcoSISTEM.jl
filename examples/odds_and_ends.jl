
#save("macrorun_new.jld", "storage_change", storage_change, "storage_steady", storage_steady)
storage_steady = load("macrorun_new.jld", "storage_steady")
storage_change = load("macrorun_new.jld", "storage_change")
storage = cat(3, storage_steady, storage_change)
temps = unique(eco.abenv.habitat.matrix)
hab =eco.abenv.habitat.matrix
@rput hab
im = eco.abundances.grid
@rput im; @rput temps; @rput vars

R"library(fields);par(mfrow=c(2,2))
for(i in c(1:4)){
  image.plot(temps, 1:10, im[i, , ], main =paste('Niche width =', vars[i]))
  }"
divtimes = collect(1:10:600)
alphas = zeros(400, 11, length(divtimes))
for i in 1:length(divtimes)
  met = Metacommunity(storage[:,:,divtimes[i],1], eco.spplist.types)
  alphas[:,:,i] = norm_sub_alpha(met, 0:10)[:diversity]
end
alphas = reshape(alphas, 20, 20, 11, length(divtimes))
temps = unique(eco.abenv.habitat.matrix)
hab =eco.abenv.habitat.matrix
@rput hab

@rput alphas; @rput temps; @rput vars
for i in 1:length(divtimes)
  @rput i
R"library(fields);par(mfrow=c(1,1))
jpeg(paste('plots/alphachange', i, '.jpg'), quality=100)
im = alphas[ , , 1, i]
im[is.na(im)] = 0
image.plot(1:20, 1:20, t(im),col=magma(25), breaks = c(0:25),
xlab='', ylab='')
dev.off()
"
end
gammas = zeros(400, 11, length(divtimes))
for i in 1:length(divtimes)
  met = Metacommunity(storage[:,:,divtimes[i],1], eco.spplist.types)
  gammas[:,:,i] = sub_gamma(met, 0:10)[:diversity]
end
gammas = reshape(gammas, 20, 20, 11, length(divtimes))
temps = unique(eco.abenv.habitat.matrix)
@rput gammas; @rput temps; @rput vars
for i in 1:length(divtimes)
  @rput i
R"library(fields);par(mfrow=c(1,1))
jpeg(paste('plots/gammachange', i, '.jpg'), quality=100)
im = gammas[ , , 1, i]
im[is.na(im)] = 0
image.plot(1:20, 1:20, t(im),col=magma(8500), breaks = c(0:8500),
xlab='', ylab='')
dev.off()
"
end

# Temporal change of beta
tempbetas = zeros(11, 400)
for sub in 1:400
  if sum(storage[:, sub, divtimes, 1])== 0.0
    tempbetas[:, sub] = repmat([0.0], 11)
  else
    met = Metacommunity(storage[:, sub, divtimes, 1], eco.spplist.types)
    tempbetas[:, sub] = norm_meta_beta(met, 0:10)[:diversity]
  end
end
tempbetas = reshape(tempbetas, 11, 20, 20)
@rput tempbetas;
R"library(fields);par(mfrow=c(1,1))
jpeg('plots/tempbeta.jpg', quality=100)
im = tempbetas[ 1 , , ]
im[is.na(im)] = 0
image.plot(1:20, 1:20, t(im),col=magma(20), breaks = seq(0,1.1,length.out=21),
xlab='', ylab='')
dev.off()
"

eco.abenv.habitat.change.rate = 0.01
for i in 1:300
  temps = unique(eco.abenv.habitat.matrix)
  hab =eco.abenv.habitat.matrix
  @rput hab; @rput i
  if any(i.==divtimes)
    R"library(viridis);library(fields)
    jpeg(paste('plots/tempchange', i, '.jpg'))
    image.plot(1:20, 1:20, t(hab),col=magma(30), breaks = c(-5:25),
    xlab='', ylab='')
    dev.off()"
  end
  TempChange(eco)
end

gammas = zeros(400, 11, length(divtimes))
for i in 1:length(divtimes)
  met = Metacommunity(storage[:,:,divtimes[i],1], eco.spplist.types)
  rhos[:,:,i] = norm_sub_rho(met, 0:10)[:diversity]
end
rhos = reshape(rhos, 20, 20, 11, length(divtimes))

plotdiv(norm_sub_beta, eco, 1)

datf = norm_sub_beta(eco, 1)
size(datf, 1) == length(eco.abenv.habitat.matrix) ||
  error("Metacommunity measures cannot be plotted as grid")
im = reshape(datf[:diversity], size(eco.abenv.habitat.matrix))
hab = eco.abenv.habitat.matrix
@rput im; @rput hab
R"par(mfrow=c(1,2));library(fields);
im[im==0]=NA
image.plot(im)
image.plot(hab)"


abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.1)
rel = TraitRelationship(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)

#plotdiv(norm_meta_alpha, eco, collect(0:10))

times = 10; burnin = 0; interval = 10; reps = 1
lensim = length(0:interval:times)

# Run simulations 10 times
storage = generate_storage(eco, lensim, reps)
for j in 1:reps
  if (j != 1) repopulate!(eco) end
  thisstore = @view storage[ :, :, :, j]
  simulate!(eco, burnin, interval, timestep)
  hab = gethabitat(eco).matrix
  heatmap(hab)
end

## New plots
storage_steady = load("macrorun_2steady.jld", "storage_steady")[:, :, 1:101, 1]
storage_change = load("macrorun_2steady.jld", "storage_change")
storage_newsteady = load("macrorun_2steady.jld", "storage_newsteady")[:, :, 1:101, 1]
storage = cat(3, storage_steady, storage_change, storage_newsteady)

divtimes = collect(1:10:1201)
alphas = zeros(400, 11, length(divtimes))
for i in 1:length(divtimes)
  met = Metacommunity(storage[cold,:,divtimes[i],1], UniqueTypes(length(cold)))
  alphas[:,:,i] = norm_sub_alpha(met, 0:10)[:diversity]
end
alphas = reshape(alphas, 20, 20, 11, length(divtimes))
temps = unique(eco.abenv.habitat.matrix)
hab =eco.abenv.habitat.matrix
@rput hab

@rput alphas; @rput temps; @rput vars
for i in 1:length(divtimes)
  @rput i
R"library(fields);par(mfrow=c(1,1))
jpeg(paste('plots/alphachangecold', i, '.jpg'), quality=100)
im = alphas[ , , 2, i]
im[is.na(im)] = 0
image.plot(1:20, 1:20, t(im),col=magma(9), breaks = seq(0,2, length.out=10),
xlab='', ylab='')
dev.off()
"
end
gammas = zeros(400, 11, length(divtimes))
for i in 1:length(divtimes)
  met = Metacommunity(storage[:,:,divtimes[i],1], eco.spplist.types)
  gammas[:,:,i] = sub_gamma(met, 0:10)[:diversity]
end
gammas = reshape(gammas, 20, 20, 11, length(divtimes))
temps = unique(eco.abenv.habitat.matrix)
@rput gammas; @rput temps; @rput vars
for i in 1:length(divtimes)
  @rput i
R"library(fields);par(mfrow=c(1,1))
jpeg(paste('plots/gammachange', i, '.jpg'), quality=100)
im = gammas[ , , 1, i]
im[is.na(im)] = 0
image.plot(1:20, 1:20, t(im),col=magma(5200), breaks = c(0:5200),
xlab='', ylab='')
dev.off()
"
end

tempbetas = zeros(11, 400)
for sub in 1:400
  if sum(storage[:, sub, divtimes[11:111], 1])== 0.0
    tempbetas[:, sub] = repmat([0.0], 11)
  else
    met = Metacommunity(storage[:, sub, divtimes[11:111], 1], eco.spplist.types)
    tempbetas[:, sub] = norm_meta_beta(met, 0:10)[:diversity]
  end
end
tempbetas = reshape(tempbetas, 11, 20, 20)
@rput tempbetas;
R"library(fields);par(mfrow=c(1,1))
jpeg('plots/tempbeta.jpg', quality=100)
im = tempbetas[ 1 , , ]
im[im==0] = NA
image.plot(1:20, 1:20, t(im),col=magma(20), breaks = seq(0,2,length.out=21),
xlab='', ylab='')
dev.off()
"
