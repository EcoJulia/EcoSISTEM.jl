## Red list categorising all plant species in GBIF ##
#setup parallel backend to use many processors
addprocs(20)
using RCall
using JuliaDB
using ClimatePref
using DataFrames

gbif = load("Post1970GBIF")
Genera = searchdir("final", ".csv")
Genera = replace.(Genera, ".csv" => "")

@sync @parallel for gen in Genera
    genus = filter(g -> g.genus == gen, gbif)
    genus_dat = DataFrame(collect(genus))
    @rput genus_dat
    R"library(rCAT)
    MissingDataFun <- function(species, dat){
    newdat = dat[dat$species == species,]
    all(is.na(newdat$decimallatitude) & is.na(newdat$decimallongitude))
    }
    genus = genus_dat[!is.na(genus_dat$species)&(genus_dat$species != ''),]
    Any_missing = unlist(sapply(unique(genus$species),
                              MissingDataFun, dat = genus))
                              NE_species = data.frame(taxa = unique(Genus$species)[Any_missing],
                                                    NOP = rep(0, sum(Any_missing)),
                                                    MER = rep(NA, sum(Any_missing)),
                                                    EOOkm2 = rep(NA, sum(Any_missing)),
                                                    AOO2km = rep(NA, sum(Any_missing)),
                                                    EOOcat = rep('NE', sum(Any_missing)),
                                                    AOOcat = rep('NE', sum(Any_missing)))
    # Remove rows with missing lat/lon
    genus = genus[!is.na(genus$decimallatitude) & !is.na(genus$decimallongitude), ]
    # Use rCAT to calculate EOO and AOO (with 2km cell size)
    GenusCat = rCAT::ConBatch(genus$species, genus$decimallatitude,
                         genus$decimallongitude, cellsize = 2000,
                         project2gether = FALSE)
    # Add on not evaluated (NE) species
    GenusCat = rbind(GenusCat, NE_species)
    "
    @rget GenusCat
end
