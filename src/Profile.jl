filename = pop!(ARGS)
using Coverage
using JLD



function saveprofile(name::String)
 allocation = analyze_malloc(Pkg.dir("Simulation"))
 filename = string(name, ".jld")
 save(Pkg.dir("Simulation", filename), "allocation", allocation)
end

saveprofile(filename)
