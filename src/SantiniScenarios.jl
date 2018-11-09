using StatsBase
function trait_populate!(ml::GridLandscape, spplist::SpeciesList,
                       abenv::AbstractAbiotic)
  # Calculate size of habitat
  dim = _getdimension(abenv.habitat)
  len = dim[1] * dim[2]
  grid = collect(1:len)
  # Set up copy of budget
      b = reshape(copy(_getbudget(abenv.budget)), size(grid))
      units = unit(b[1])
      activity = reshape(copy(abenv.active), size(grid))
      b[.!activity] = 0.0 * units
      # Loop through species
      for i in eachindex(spplist.abun)
        # Get abundance of species
        if spplist.native[i]
        abun = spplist.abun[i]
        options = unique(abenv.habitat.matrix)
        pref = spplist.traits.val[i]
        if (pref ∉ options)
            wv = rand(Beta(2,2), len)
        else
            wv= Vector{Float64}(len)
            wv[reshape(abenv.habitat.matrix, len).==pref]= 1.0
            wv[reshape(abenv.habitat.matrix, len).!=pref]= 0.0
        end
        # Randomly choose position on grid (weighted)
        pos = sample(grid[b .> (0 * units)], weights(wv), abun, replace=true)
        # Add individual to this location
        for j in pos
            ml.matrix[i, j] = ml.matrix[i, j] + 1
        end
      end
  end
end

function trait_repopulate!(eco::Ecosystem)
  eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
  eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun), length(eco.spplist.abun)))
  trait_populate!(eco.abundances, eco.spplist, eco.abenv)
end
