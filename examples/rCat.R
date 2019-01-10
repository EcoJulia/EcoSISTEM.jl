## Red list categorising all plant species in GBIF ##

library(foreach)
library(doParallel)
library(rCAT)

#setup parallel backend to use many processors
cl <- makeCluster(20)
registerDoParallel(cl)

# Folder name housing csv files
dirname = "final/"
files <- dir(dirname)

# Function for identifying missing data points in lat/lon of a species in
# the dataset
MissingDataFun <- function(species, dat){
  newdat = dat[dat$species == species,]
  all(is.na(newdat$decimallatitude) & is.na(newdat$decimallongitude))
}

# Parallel loop through species, calculating EOO and AOO category
ans = foreach(dat = files, .combine='rbind') %dopar% {
  # Read file and remove any with missing species names
  Genus = read.csv(paste0(dirname, dat))
  Genus = Genus[!is.na(Genus$species)&(Genus$species != ""),]
  # Use MissingDataFun to identify species that have no lat/lon coordinates
  Any_missing = unlist(sapply(unique(Genus$species),
                              MissingDataFun, dat = Genus))
  # Create a dataframe for these species and assign them as Not Evaluated
  NE_species = data.frame(taxa = unique(Genus$species)[Any_missing],
                          NOP = rep(0, sum(Any_missing)),
                          MER = rep(NA, sum(Any_missing)),
                          EOOkm2 = rep(NA, sum(Any_missing)),
                          AOO2km = rep(NA, sum(Any_missing)),
                          EOOcat = rep("NE", sum(Any_missing)),
                          AOOcat = rep("NE", sum(Any_missing)))
  # Remove rows with missing lat/lon
  GenusCoords = Genus[!is.na(Genus$decimallatitude) & !is.na(Genus$decimallongitude), ]
  # Use rCAT to calculate EOO and AOO (with 2km cell size)
  GenusCat = rCAT::ConBatch(GenusCoords$species, GenusCoords$decimallatitude,
                            GenusCoords$decimallongitude, cellsize = 2000,
                            project2gether = FALSE)
  # Add on not evaluated (NE) species
  GenusCat = rbind(GenusCat, NE_species)
  matchnames <- function(dat, x) sum(dat$species == x)
  totalrecords = sapply(GenusCat$taxa, matchnames, dat = Genus)
  return(cbind(GenusCat, totalrecords))
}
# Write the dataframe as csv and stop workers
write.csv(ans, "Conservation_status.csv")
stopCluster(cl)
