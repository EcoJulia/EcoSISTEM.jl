using Diversity
import Diversity.sumovertypes, Diversity.sumoversubcommunities, Diversity.powermean

#function sumovertypes{A, H}(ml::MatrixLandscape{A, H}, vals::A)
#  mapslices(sum, vals, 1)::A
#end

#function sumoversubcommunities{A, H}(ml::MatrixLandscape{A, H}, vals::A)
#  mapslices(sum, vals, (2,3))::A
#end

#function powermean{S <: AbstractFloat}(values::AbstractArray{S,3}, orders, weights::AbstractArray{S,3})
#  arr = Matrix{Float64}(size(values, 2), size(values, 3))
#  for (i = 1:size(arr, 1))
#    for (j = 1:size(arr, 2))
#      arr[i, j] = powermean(values[:, i, j], orders, weights[:, i, j])
#    end
#  end
#  arr
#end
