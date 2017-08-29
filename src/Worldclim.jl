using FileIO
using Images
using ImageView
using RasterIO

img = RasterIO.openraster(Pkg.dir("Simulation", "examples/data/Worldclim",
 "wc2.0_bio_30s_01.tif"))
imshow(img)
