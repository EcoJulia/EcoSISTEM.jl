
using Distributions
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
using JLD
using ClimatePref
using Simulation
using Diversity
using TOML
using JuliaDB
# Dictionary object of types
typedict = Dict("TempBin" => TempBin, "RainBin" => RainBin,
"GaussianKernel" => GaussianKernel, "Worldclim" => Worldclim,
"ERA-interim" => ERA,"CERA-20C" => CERA, "Trapeze" => Trapeze, "Uniform" => Unif)
# Dictionary object of functions
funcdict = Dict("norm_sub_alpha" => norm_sub_alpha)

# Dictionary object of units
unitdict = Dict("month" => month, "year" => year, "km" => km, "mm" => mm,
"day^-1kJm^-2" => day^-1*kJ*m^-2, "deg" => °)

function readoutput(file::String, name::String)
    filenames = searchdir(file, name)
    joinname = joinpath(file, filenames[1])
    withoutend = split(joinname, ".jld")[1]
    withoutbegin = split(withoutend, file)[2]
    mat = JLD.load(joinname, String(withoutbegin))
    for i in eachindex(filenames)[2:end]
        joinname = joinpath(file, filenames[i])
        withoutend = split(joinname, ".jld")[1]
        newmat = JLD.load(joinname, String(withoutend))
        mat = cat(3, mat, newmat)
    end
    return mat
end
function runTOML(file::String, eco::Ecosystem)
    fulldict = TOML.parsefile(file)
    params = fulldict["params"]
    datadump = fulldict["datadump"]

    unit = unitdict[params["unit"]]
    burnin = params["burnin"] * unit
    times = params["times"] * unit
    interval = params["interval"] * unit
    timestep = params["timestep"] * unit

    dumpinterval = datadump["interval"] * unit
    outfile = datadump["outfile"]
    divides = Int(times/dumpinterval)

    if haskey(fulldict, "measure")
        measure = fulldict["measure"]
        qs = measure["qs"]
        divfun = funcdict[measure["measures"]]
        simulate!(eco, burnin, timestep)
        resettime!(eco)
        for i in 1:divides
            lensim = ifelse(i == 1, length(0month:interval:dumpinterval),
                length(timestep:interval:dumpinterval))
            abun = generate_storage(eco, 1, lensim, 1)
            simulate_record_diversity!(abun, eco, dumpinterval, interval, timestep,
                divfun, qs)
            JLD.save(string(outfile, "Run", i, ".jld"), string("Run", i), abun)
        end
    else
        simulate!(eco, burnin, timestep)
        resettime!(eco)
        for i in 1:divides
            lensim = ifelse(i == 1, length(0month:interval:dumpinterval),
                length(timestep:interval:dumpinterval))
            abun = generate_storage(eco, lensim, 1)
            simulate_record!(abun, eco, dumpinterval, interval, timestep)
        end
    end
    numSpecies = size(eco.abundances.matrix, 1)
    print(string("Model run:", "\n", times, "\n",
        numSpecies, " species", "\n", burnin, " burnin",
        "\n", interval, " recording interval", "\n",
        timestep, " timestep", "\n", "Output: ", outfile))
end

function readTOML(file::String)
    fulldict = TOML.parsefile(file)
    dir = fulldict["dir"]
    # Traits
    sppl = TOML_sppl(fulldict)
    abenv = TOML_abenv(fulldict)
    locations = TOML_locs(fulldict, sppl)
    rel = TOML_rel(fulldict)
    eco = Ecosystem(locations, sppl, abenv, rel)
    return eco
end

function TOML_sppl(fulldict::Dict)
    dir = fulldict["dir"]
    species = fulldict["species"]
    requirements = fulldict["requirements"]
    movement = fulldict["movement"]
    # Traits
    if typeof(species["names"]) == Vector{String}
        names = species["names"]
    else
        numspp = species["names"]
        gbif = JuliaDB.load(joinpath(dir, species["file"]))
        chunk = Int(round(length(gbif)/numspp))
        names = Vector{String}(numspp)
        for i in 1:numspp
            sel = rand(1:1000, 1)
            names[i] = rows(gbif, :species)[sel][1]
        end
    end
    numSpecies = length(names)
    traitcol = TOML_trait(fulldict, names)
    # Requirements
    req = TOML_requirement(requirements, numSpecies)
    movements = TOML_move(movement, numSpecies)
    if species["allnative"]
        native = repmat([true], numSpecies)
    else
        error("Need method for individual native status")
    end
    params = species["params"]
    abun = Multinomial(params["individuals"], numSpecies)
    param = EqualPop(params["birth"]/month, params["death"]/month, params["long"], params["surv"], params["boost"])
    sppl = SpeciesList(numSpecies, traitcol, abun, req,
        movements, param, native)
    sppl.names = names
    return sppl
end
function TOML_trait(fulldict::Dict, names::Vector{String})
    dir = fulldict["dir"]
    traits = fulldict["traits"]
    species = fulldict["species"]
    numTraits = length(traits)

    if haskey(traits, "temperature")
        tempdict = traits["temperature"]
        filename = joinpath(dir, tempdict["file"])
        Temp = JLD.load(filename, tempdict["name"])
         tp = typedict[tempdict["type"]]
        trait1 = tp(Array(transpose(Temp[names,:])))
    end
    if haskey(traits, "rainfall")
        raindict = traits["rainfall"]
        filename = joinpath(dir,raindict["file"])
        Rain = JLD.load(filename, raindict["name"])
        tp = typedict[raindict["type"]]
        trait2 = tp(Array(transpose(Rain[names,:])))
    end
    if numTraits == 2
        traitcol = TraitCollection2(trait1, trait2)
    else
        traitcol = ifelse(haskey(traits, "temperature"), traits1, traits2)
    end
    return traitcol
end
function TOML_requirement(requirements::Dict, numSpecies::Int64)
    numReqs = length(setdiff(keys(requirements), ["allsame"]))
    if haskey(requirements, "solar")
        if requirements["allsame"]
            unit = unitdict[requirements["solar"]["unit"]]
            req1 = SolarRequirement(repmat([requirements["solar"]["req"] * unit],
                numSpecies))
        else
            error("Need method for individual requirements")
        end
    end
    if haskey(requirements, "water")
        if requirements["allsame"]
            unit = unitdict[requirements["water"]["unit"]]
            req2 = WaterRequirement(repmat([requirements["water"]["req"] * unit],
                numSpecies))
        else
            error("Need method for individual requirements")
        end
    end
    if numReqs == 2
        req =  ReqCollection2(req1, req2)
    else
        req = ifelse(haskey(requirements, "solar"), req1, req2)
    end
    return req
end
function TOML_move(movement::Dict, numSpecies::Int64)
    # Movement
    kernel = typedict[movement["type"]]
    distance = movement["distance"] * unitdict[movement["unit"]]
    thresh = movement["threshold"]
    movements = BirthOnlyMovement(kernel(distance, numSpecies, thresh))
    return movements
end

function TOML_abenv(fulldict::Dict)

    dir = fulldict["dir"]
    climate = fulldict["climate"]
    budget = fulldict["budget"]

    grid = fulldict["grid"]
    unit = unitdict[grid["unit"]]
    minX = grid["minX"] * unit
    maxX = grid["maxX"] * unit
    minY = grid["minY"] * unit
    maxY = grid["maxY"] * unit

    numHabs = length(climate)
    numEnergy = length(budget)

    # Climate
    if haskey(climate, "temperature")
        tempdict = climate["temperature"]
        filename = joinpath(dir, tempdict["file"])
        tp = typedict[tempdict["type"]]
        if tp == ERA
            hab1 = extractERA(filename, tempdict["name"],
                collect(1.0month:1month:tempdict["years"]*year))
            hab1.array = hab1.array[:, -89.25° .. 90°][minX .. maxX, minY .. maxY]
        elseif tp == CERA
            hab1 = extractCERA(joinpath(fulldict["dir"], tempdict["file"]), tempdict["search"], tempdict["name"])
            step = hab1.array.axes[1][2] - hab1.array.axes[1][1]
            hab1.array = hab1.array[:, -89.25° .. 90°][minX .. maxX, minY - step .. maxY]
        else
            hab1 = extractworldclim(joinpath(filename, tempdict["name"]))
            hab1.array = hab1.array[minX .. maxX, minY .. maxY]
        end
        if haskey(tempdict, "upscale")
            hab1 = upresolution(hab1, tempdict["upscale"])
        end

    end
    if haskey(climate, "rainfall")
        raindict = climate["rainfall"]
        filename = joinpath(dir, raindict["file"])
        tp = typedict[raindict["type"]]
        if tp == ERA
            hab2 = extractERA(filename, raindict["name"],
            collect(1.0month:1month:raindict["years"]*year))
            hab2.array = hab2.array[:, -89.25° .. 90°][minX .. maxX, minY .. maxY]
        elseif tp == CERA
            hab2 = extractCERA(joinpath(fulldict["dir"], tempdict["file"]), tempdict["search"], raindict["name"])
            hab2.array = hab2.array[:, -89.25° .. 90°][minX .. maxX, minY - step .. maxY]
        else
            hab2 = extractworldclim(joinpath(filename, raindict["name"]))
            hab2.array = hab2.array[minX .. maxX, minY .. maxY]
        end
    if haskey(raindict, "upscale")
        hab2 = upresolution(hab2, raindict["upscale"])
    end

    end
    # Budget
    if haskey(budget, "solar")
        solar = budget["solar"]
        filename = joinpath(dir, solar["file"])
        srad = extractworldclim(joinpath(filename, solar["name"]))
        srad.array = srad.array[minX .. maxX, minY .. maxY]
        if haskey(solar, "upscale")
            srad = upresolution(srad, solar["upscale"])
        end
        budget1 = SolarBudget(convert(Array{typeof(2.0*day^-1*kJ*m^-2),3},
            srad.array), 1)
    end
    if haskey(budget, "water")
        water = budget["water"]
        if water["fromclimate"]
            budget2 = WaterBudget(Array{typeof(1.0*mm), 3}(hab2.array), 1)
        else
            filename = joinpath(dir, water["file"])
            water = extractworldclim(joinpath(filename, water["name"]))
                        water.array = water.array[minX .. maxX, minY .. maxY]
            if haskey(budget["water"], "upscale")
                water = upresolution(water, budget["water"]["upscale"])
            end
            budget2 = WaterBudget(convert(Array{typeof(2.0*mm),3},
                water.array), 1)
        end
    end
    # Active
    active = fulldict["active"]
    filename = joinpath(dir, active["file"])
    world = extractfile(filename)[:, -89.25° .. 90°]
    world = world[minX .. maxX, minY .. maxY]
    if haskey(active, "upscale")
        world = upresolution(world, active["upscale"])
    end
    activegrid = Array{Bool, 2}(.!isnan.(ustrip.(world)))


    if numEnergy == 2
        bud = BudgetCollection2(budget1, budget2)
    else
        bud = ifelse(haskey(budget, "solar"), budget1, budget2)
    end
    if numHabs == 2
        abenv1 = eraAE(hab1, budget1, activegrid)
        abenv2 = worldclimAE(hab2, budget1, activegrid)
        hab = HabitatCollection2(abenv1.habitat, abenv2.habitat)
    else
        hab = ifelse(haskey(climate, "temperature"),
            eraAE(hab1, budget1, activegrid),
            worldclimAE(hab2, budget1, activegrid))
    end


    abenv = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, abenv1.active,
            bud, abenv1.names)
end

function TOML_locs(fulldict::Dict, sppl::SpeciesList)
    dir = fulldict["dir"]
    grid = fulldict["grid"]
    locs = fulldict["locs"]
    names = sppl.names
    refgrid = Array{Int64,2}(Tuple(grid["gridsize"]))
    unit = unitdict[grid["unit"]]
    x = (grid["minX"]*unit):(grid["resolution"]*unit):(grid["maxX"]*unit)
    y =  (grid["minY"]*unit-grid["resolution"]*unit):(grid["resolution"]*unit):(grid["maxY"]*unit)
    refgrid =  AxisArray(refgrid, Axis{:longitude}(x), Axis{:latitude}(y))
    refgrid[1:end] = 1:length(refgrid)
    ref = Reference(refgrid)
    tab = JuliaDB.load(joinpath(dir,locs["file"]))
    filtab = filter(t -> t.species in names, tab)
    ctab = collect(filtab)
    lon = select(ctab, :decimallongitude)
    lat = select(ctab, :decimallatitude)
    locs = find((lon .> grid["minX"]) .& (lon .< grid["maxX"]) .&
        (lat .> grid["minY"]) .&  (lat .< grid["maxY"]))
    vals = extractvalues(lon[locs] * °, lat[locs] * °, ref)
    locations = table(collect(select(ctab, :species))[locs, :], vals)
    return locations
end
function TOML_rel(fulldict::Dict)
    rel = fulldict["relationship"]
    numRel = length(rel) - 1
    if haskey(rel, "temperature")
        tp = typedict[rel["temperature"]["type"]]
        rel1 = tp{eltype(1.0°C)}()
    end
    if haskey(rel, "rainfall")
        tp = typedict[rel["rainfall"]["type"]]
        rel2 = tp{eltype(1.0mm)}()
    end
    if numRel == 2 && rel["multiplicative"]
        rel = multiplicativeTR2(rel1, rel2)
    else
        rel = ifelse(haskey(rel, "temperature"), rel1, rel2)
    end
    return rel
end
