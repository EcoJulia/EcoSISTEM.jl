"""
    GridLandscape

Ecosystem abundances housed in the landscape. These are represented in both 2
dimensions (for computational efficiency in simulations) and 3 dimensions (to
represent species, their abundances and position in the grid).

"""
mutable struct GridLandscape
  matrix::Matrix{Int64}
  grid::Array{Int64, 3}

  function GridLandscape(abun::Matrix{Float64}, dimension::Tuple)
    return new(abun, reshape(abun, dimension))
  end
end

"""
    emptygridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)

Function to create an empty GridLandscape given a GridAbioticEnv and a
SpeciesList.
"""
function emptygridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)
  mat = zeros(counttypes(spplist, true), countsubcommunities(gae))

  dimension = (counttypes(spplist, true), size(gae.habitat.matrix)...)
  return GridLandscape(mat, dimension)
end
