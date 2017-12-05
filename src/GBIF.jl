using DataStructures
#using ArrayViews
using DataFrames
using JLD
using RCall

addprocs(1)
@everywhere using JuliaDB
latlon = loadndsparse("/Users/claireh/Documents/PhD/Data/GBIF/final",
       indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
       type_detect_rows =5000, chunks = 1)
save(latlon, "/Users/claireh/Documents/PhD/Data/GBIF/store/latlon")
#latlon = loadtable("/Users/claireh/Documents/PhD/Data/GBIF/CSV",
#       indexcols = [:decimallatitude, :decimallongitude],usecache = false, type_detect_rows =5000)
coords = select(latlon, :decimallatitude)

latlon = readtable(
  "/Users/claireh/Documents/PhD/Data/GBIF/0000274-170817152713382.csv",
  nrows=100000, header=true, separator= '\t')
sum(isna.(latlon[:decimallatitude]))/nrow(latlon)

coords = Array{Float64, 2}(latlon[[:decimallatitude, :decimallongitude]])
acc = mapslices(points_accuracy, coords, 2)
tab = countmap(acc)
acc_count = Vector{Int}(8)
count = 0
for i in keys(tab)
    count += 1
    acc_count[count] = tab[i]
end

resol = mapslices(res, coords, 2)
least_resol = mapslices(maximum, resol, 2)

@rput least_resol
R"
breaks=c(0, 0.1,0.5,1,2,10,50,100,150)
r.c <- cut(least_resol,breaks=breaks)
barplot(table(r.c), ylim=c(0,100000),
names=c('0 - 100m','100 - 500m', '500m - 1km', '1 - 2km' ,
'2 - 10km' ,' 10 - 50km' , '50 - 100km' , '100 - 150km'))"
Species = OrderedDict(sort(collect(countmap(latlon[:species])), by=x->x[2], rev=true))

#' Identify common species
common = collect(keys(Species))[2:11]
#' Identify rare (ish) species
#' CR category in iucn is defined as 50 individuals or less
critically_endangered = find(collect(values(Species)).> 5 .& collect(values(Species)) .< 50)
## sample ten of these randomly to match the common ones
rare = collect(keys(Species))[sample(critically_endangered, 10)]
coords_com = least_resol[vcat(map(x -> find(latlon[:species].== x), common)...),:]
coords_rare = least_resol[vcat(map(x -> find(latlon[:species].== x), rare)...),:]

@rput coords_com
@rput coords_rare
R"
breaks = c(0,0.1,0.5,1,2,10,50,100,150)
r.c.c <- cut(coords_com,breaks=breaks)
r.c.r <- cut(coords_rare,breaks=breaks)

tab.c <- table(r.c.c)
names(tab.c) <- c('0 - 100m','100 - 500m', '500m - 1km', '1 - 2km' ,
'2 - 10km' ,' 10 - 50km' , '50 - 100km' , '100 - 150km')
print(tab.c)
tab.r <- table(r.c.r)
names(tab.r) <- c('0 - 100m','100 - 500m', '500m - 1km', '1 - 2km' ,
'2 - 10km' ,' 10 - 50km' , '50 - 100km' , '100 - 150km')
print(tab.r)
dat <- cbind(rbind(data.frame(tab.r),data.frame(tab.c)))
dat <- cbind(dat, Rarity=c(rep('Rare',8),rep('Common',8)))
library(ggplot2)
library(cowplot)
g1 <- ggplot(data.frame(tab.r),aes(x=Var1, y=Freq))+
  geom_bar(stat='identity',fill='#00BFC4')+
  xlab('Accuracy')+ylab('Frequency') +
  theme(axis.text.x = element_text(size=22,color='darkgrey',
  angle = 60, vjust=0),
        axis.text.y = element_text(size=22,color='darkgrey'),
        axis.title=element_text(size=24,face='bold'),
        legend.position='none',
  plot.margin = unit(c(3,3,3,3), 'cm'))
g1

g2 <- ggplot(data.frame(tab.c),aes(x=Var1, y=Freq))+
  geom_bar(stat='identity',fill='#F8766D')+
  xlab('Accuracy')+ylab('Frequency')+
  theme(axis.text.x = element_text(size=22,color='darkgrey',
  angle = 60, vjust=0),
        axis.text.y = element_text(size=22,color='darkgrey'),
        axis.title=element_text(size=24,face='bold'),
        legend.position='none',
  plot.margin = unit(c(3,3,3,3), 'cm'))
g2


g3=ggplot(dat, aes(x=Var1, y=Freq, fill=Rarity))+
  geom_bar(stat='identity',position = position_dodge(),width=1)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend <- g_legend(g3)
## Save plot to pdf
ggsave(filename='Rare vs common.pdf',plot=gridExtra::grid.arrange(g1,g2,legend,
 ncol=3, widths=c(0.8,0.8,0.1)), scale=3)
"

function decimal_places(coord::Float64, max::Int64 = 10)
    dp = 0
    while (abs(round(coord, dp) - coord) != 0)
        dp += 1
    end
  dp
end

# Function that finds how many significant figures a coordinate has
#' Inputs:
#'
#' - A lat/long location
#'
#' Returns:
#'
#' - The number of signficant figures
function signif_figs(coord::Float64)
  ## If you round the coordinate to a specific number of sf, does it equal itself?
  coord = abs(coord)
  sf = 1
  while (abs(signif(coord, sf) - coord) != 0)
      sf += 1
  end
  sf
end
#' Function that finds to what accuracy GPS locations
#' have been written down to
#'
#' Inputs:
#'
#' - coords -- A lat/long location
#'
#' Returns:
#'
#' - The least accurate level the points were written to

function points_accuracy(coords::Vector{Float64})
  ## How many decimal places?
  dp = map(decimal_places, coords)
  ## How many significant figures?
  sf = map(signif_figs, coords)

  degrees = trunc.(coords)
  mins = trunc.(map(x -> mod(x, 60), coords*60))
  secs = map(x -> mod(x, 60), abs(coords)*3600)

    if (any(mins .== 0) & any(secs .== 0))
  ## If to any to no decimal places then is accurate to degrees
      acc="deg"
    elseif (abs(sf[1] - sf[2]) == 0 &
           abs(dp[1] - dp[2]) > 0 &
           all(sf .< 7))
    ## If the number of significant figures are the same but decimals different then put down to sf
    acc = string(sf[1],"sf")
elseif (any(mins .> 0) &&
           any(secs .== 0))
    ## If accurate to a minute
    acc = "min"
elseif (any((round(secs) - secs) .== 0))
    ## If accurate to a second
    acc = "sec"
  else
    ## Find to what decimal place the point is recorded to
    ## Or if lower than 5 decimal places then assign to lower
    acc = ifelse(maximum(dp) > 5, "lower", string(maximum(dp), "dp"))
    end
    acc
end




function res(coords::Vector{Float64}, lat_eq::Float64 = 111.32)
  lat = coords[1]
  acc = points_accuracy(coords)
  ## If degrees then find 1 degree at that lat/long
  if (acc=="deg")
  dist_lat = lat_eq
  dist_long = abs(cos(lat)*lat_eq)
  elseif (acc=="min")
    abs(cos(lat)*lat_eq)
    ## If minutes then find 1 minute at that lat/long
    dist_lat  = lat_eq/60
    dist_long = abs(cos(lat))*lat_eq/60
  elseif (acc=="sec")
    ## If seconds then find 1 second at that lat/long
    dist_lat = lat_eq/3600
    dist_long = abs(cos(lat))*lat_eq/3600
  else
    ## If other, then find to number of decimal places at that lat/long
     dp = map(decimal_places, coords)
     dist_lat = lat_eq/(10^dp[1])
     dist_long = abs(cos(lat)*lat_eq)/(10^(dp[2]))
 end
  ## Return a vector with precision of each coordinate (in km)
  resolution = [dist_lat, dist_long]
  resolution
end


#' Function to find number in last decimal place
#'
function last_decimal(coord::Float64)
  dp = decimal_places(coord, max=10)
  ll = floor(coord*(10^dp))
  (ll-floor(ll/(10))*(10))
end


TestSpecies =
