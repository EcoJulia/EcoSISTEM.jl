using DataStructures
#using ArrayViews
using DataFrames
using JLD
dat = readtable(
  "/Users/claireh/Documents/PhD/Data/GBIF/0000274-170817152713382.csv",
  header=false)
headers = readtable(
  "/Users/claireh/Documents/PhD/Data/GBIF/headers.csv",
  header=false)
names!(dat, map(a-> convert(Symbol, a), Array(headers)[1,:]))
sum(isna.(dat[:decimallatitude]))/nrow(dat)


dat[(.!isna.(dat[:decimallatitude])).&
(.!isna.(dat[:decimallongitude])), :]
