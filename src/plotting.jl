using RCall

function plot_move(eco::Ecosystem, x::Int64, y::Int64, spp::Int64, plot::Bool)

  lookup = eco.lookup[spp]
  maxX = size(eco.abenv.habitat.matrix, 1) - x
  maxY = size(eco.abenv.habitat.matrix, 2) - y
  # Can't go over maximum dimension
  valid = find((lookup.x .> -x) .& (lookup.y .> -y) .&
   (lookup.x .<= maxX) .& (lookup.y .<= maxY))
  probs = lookup.p[valid]
  probs ./= sum(probs)
  xs = (lookup.x[valid] .+ x)
  ys = (lookup.y[valid] .+ y)
  A = zeros(size(eco.abenv.habitat.matrix))
  for i in eachindex(xs)
    A[xs[i], ys[i]] = probs[i]
  end
  if plot
    @rput A
    R"par(mfrow=c(1,1));library(fields);
    A[A==0]=NA
    image.plot(A)"
  else
    return A
  end
end

function plot_abun(abun::AbstractArray, numSpecies::Int64, gridSize::Int64)
  # Plot
  @rput abun
  @rput numSpecies
  @rput gridSize
  R" par(mfrow=c(gridSize,gridSize), mar=c(2, 2, 2, 2))
  for (i in 1:gridSize^2){
      for (k in 1:numSpecies){
        if (k==1) plot_fun=plot else plot_fun=lines
          plot_fun(0:100, abun[k, i, , 1], col=k, xlab='Abundance', ylab='Time', type='l',
          ylim=c(0, max(abun)))
        }
    }"
  end

function plot_divergence(expected::Vector{Float64}, actual::Vector{Float64})
  @rput expected
  @rput actual
  KL = kldivergence(actual, expected); @rput KL
  R"par(mfrow=c(1,1))
    plot(1:length(expected),expected, type='l',
          main = paste('Divergence =', round(KL, 2)), xlab='Abundance',
          ylab='Frequency', ylim=c(0, max(c(expected, actual))))
  abline(h=max(expected), col=1, cex=0.5, lty=3)
  lines(actual, col=2)
  abline(h=max(actual), col=2, cex=0.5, lty=3)
  legend('topright', legend=c('Expected', 'Observed'), col=1:2, pch='-')"
  info("Divergence = ",KL)
end

function plot_divergence(combined::Array{Array{Float64, 1}, 1})
  expected = combined[1]
  actual = combined[2]
  plot_divergence(expected, actual)
end

function freq_hist(grd::Array{Float64, 4}, sq::Int64, num::Int64)
  total = mapslices(sum, grd , length(size(grd)))
  grd = grd[:, :, :, sq]
  _freq_hist(total, grd, num)
end

function freq_hist(grd::Array{Float64, 3}, sq::Int64, num::Int64)
  total = mapslices(sum, grd , length(size(grd)))
  grd = grd[:, :, sq]
  _freq_hist(total, grd, num)
end

function freq_hist(grd::Array{Float64}, sq::Int64, num::Int64, spp::Int64)
  spp_grd = grd[:, spp, :, :]
  freq_hist(spp_grd, sq, num)
end

function _freq_hist(total::Array{Float64}, grd::Array{Float64}, num::Int64)
  grd = grd[reshape(total, size(grd)).>0]
  total = total[total.>0]
  count_tot = grd[total.== num]
  @rput count_tot
  @rput num
  R"hist(count_tot, breaks=c(-0.5:(num+0.5)), main=' ', xlab='Abundance')"
end

function plotdiv(divfun::Function, eco::Ecosystem, qs::Array{Real})
  datf = divfun(eco, qs)
  @rput datf
  R"library(ggplot2); library(cowplot)
  ggplot(datf, aes(x = q, y = diversity, col = partition_name)) + geom_line()"
end


function plotdiv(divfun::Function, eco::Ecosystem, qs::Real)
  datf = divfun(eco, qs)
  size(datf, 1) == length(eco.abenv.habitat.matrix) ||
    error("Metacommunity measures cannot be plotted as grid")
  im = reshape(datf[:diversity], size(eco.abenv.habitat.matrix))
  @rput im
  R"par(mfrow=c(1,1));library(fields);
  im[im==0]=NA
  image.plot(im)"
end
