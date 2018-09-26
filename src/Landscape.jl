using Missings
"""
    GridLandscape

Ecosystem abundances housed in the landscape. These are represented in both 2
dimensions (for computational efficiency in simulations) and 3 dimensions (to
represent species, their abundances and position in the grid).

"""
mutable struct GridLandscape
  matrix::Matrix{Int64}
  grid::Array{Int64, 3}
  seed::MersenneTwister

  function GridLandscape(abun::Matrix{Int64}, dimension::Tuple)
    a = abun
    return new(a, reshape(a, dimension), Base.GLOBAL_RNG)
  end
end

mutable struct CachedGridLandscape
  matrix::Vector{Union{GridLandscape, Missing}}
  outputfolder::String
  saveinterval::Unitful.Time

  function CachedGridLandscape(dim::Int64, file::String, interval::Unitful.Time)
    v = Vector{Union{GridLandscape, Missing}}(undef, dim)
    fill!(v, missing)
    return new(v, file, interval)
  end
end

"""
    emptygridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)

Function to create an empty GridLandscape given a GridAbioticEnv and a
SpeciesList.
"""
function emptygridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)
  mat = zeros(Int64, counttypes(spplist, true), countsubcommunities(gae))

  dimension = (counttypes(spplist, true), _getdimension(gae.habitat)...)
  return GridLandscape(mat, dimension)
end
