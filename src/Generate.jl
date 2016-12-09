using Distributions

function jmulti(n::Int64, p::AbstractArray)
  rand(Multinomial(n, p))
end
function jmulti(n::Int64, p::Real)
  rand(Multinomial(n, repmat([p], n)))
end



#function create_habitat(dims::Tuple, corr, prop)
#  mat=zeros(10, 10)
#  rand!(mat,[0 1])
#end



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


function SR(ecosystem::Ecosystem, individuals::Int64)
  sz=size(ecosystem.partition.abundances,2,3)
  ms=mapslices(sum, ecosystem.partition.abundances, 1)#*individuals
  reshape(ms, sz)
end
function SR(ecosystem::AbstractHabitat)
 sz=size(ecosystem.partition.abundances,2,3)
 ms=mapslices(sum, ecosystem.partition.abundances, 1)
 reshape(ms, sz)
end

function create_habitat(dim, types, prop)
  sample(types, WeightVec(prop), dim)
end

function populate(species::Int64, individuals::Int64, habitat::Niches, traits::Vector)
  dim=size(habitat.matrix)
  P=zeros(Int64,species,dim[1],dim[2])
  abun_vec=rand(Multinomial(individuals, species))
  for i in eachindex(abun_vec)
    abun=abun_vec[i]
    pref=traits[i]
    wv= Vector{Float64}(100)
    wv[find(reshape(habitat.matrix, (100,1)).==pref)]= 0.9
    wv[find(reshape(habitat.matrix, (100,1)).!=pref)]= 0.1
      while abun>0
      pos=sample(1:(dim[1]*dim[2]), weights(wv))
      P[i,pos]=P[i,pos]+1
      abun=abun-1
    end
  end
  MatrixLandscape(P, habitat)
end
