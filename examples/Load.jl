addprocs(8)
@everywhere using JuliaDB
latlon = loadndsparse("final",
       indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
       type_detect_rows =5000,
       colparsers=Dict(:datasetkey=>String,
                       :occurrenceid=>String,
                       :gbifid=>String,
                       :locality=>String,
                       :publishingorgkey=>String,
                       :taxonkey=>String,
                       :institutioncode=>String,
                       :catalognumber=>String,
                       :recordnumber=>String))
save(latlon, "output")

addprocs(4)
@everywhere using JuliaDB
latlon = loadndsparse("test",
       indexcols = [:phylum, :class, :order, :family, :genus, :species, :scientificname],
       type_detect_rows =5000,
       colparsers=Dict(:datasetkey=>String,
                       :occurrenceid=>String,
                       :gbifid=>String,
                       :locality=>String,
                       :publishingorgkey=>String,
                       :taxonkey=>String,
                       :institutioncode=>String,
                       :catalognumber=>String,
                       :recordnumber=>String))
save(latlon, "output")
