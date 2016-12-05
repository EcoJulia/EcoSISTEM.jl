using Distributions
using DataFrames

function jmulti(n::Int64, p::AbstractArray)
  rand(Multinomial(n, p))
end
function jmulti(n::Int64, p::Real)
  rand(Multinomial(n, repmat([p], n)))
end

function neighbours(x, y, dim)

function create_habitat(dims::Tuple, corr, prop)
  mat=zeros(10, 10)
  rand!(mat,[0 1])
end

mat=ones(10, 10)

dims=size(mat)
squares=dims[1]*dims[2]
function populate(species::Int64, individuals::Int64, habitat::Habitats)
  dim=size(habitat.matrix)
  P=zeros(Int64,species,dim[1],dim[2])
  abun_vec=rand(Multinomial(individuals, species))
  for i in eachindex(abun_vec)
    abun=abun_vec[i]
      while abun>0
      x=rand(1:dim[1])
      y=rand(1:dim[1])
      P[i,x,y]=P[i,x,y]+1
      abun=abun-1
    end
  end
  MatrixLandscape(P, habitat)
end
