dir = "/Users/claireh/Documents/PhD/Data/Worldclim/wc2.0_5m/"
filenames = list.files(paste(dir, "wc2.0_5m_srad", sep=""))[2:13]
for (i in seq_along(filenames)){
  r = raster(paste(dir, "wc2.0_5m_srad/", filenames[i], sep=""))
  rag = aggregate(r, fact = 9)
  writeRaster(rag, paste("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc/",
             filenames[i], sep=""), overwrite=T)
}

dir = "/Users/claireh/Documents/PhD/Data/Worldclim/wc2.0_5m/"
filenames = list.files(paste(dir, "wc2.0_5m_prec", sep=""))[2:13]
for (i in seq_along(filenames)){
  r = raster(paste(dir, "wc2.0_5m_prec/", filenames[i], sep=""))
  rag = aggregate(r, fact = 9)
  writeRaster(rag, paste("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc/wc2.0_5m_prec/",
                         filenames[i], sep=""), overwrite=T)
}