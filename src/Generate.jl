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
function get_neighbours(mat::Matrix, y::Int64, x::Int64, chess::Int64=4)
  dims=size(mat)
  if chess==4
    neighbour_vec=[x y-1; x y+1; x-1 y; x+1 y]
  elseif chess==8
    neighbour_vec=[x y-1; x y+1; x-1 y; x+1 y; x-1 y-1; x-1 y+1; x+1 y-1; x+1 y+1]
  else
    error("Can only calculate neighbours in 4 or 8 directions")
  end
  remove=vcat(mapslices(all, [neighbour_vec.>=1 neighbour_vec[:,1].<=dims[1] neighbour_vec[:,2].<=dims[2]], 2)...)
  neighbour_vec=neighbour_vec[remove,:]
  neighbour_vec
end
function update!(eco::Ecosystem,  birth::Float64, death::Float64, move::Float64, timestep::Real)
  abun=map(i->sum(eco.partition.abundances[i,:,:]), 1:size(eco.partition.abundances,1))
  dims=size(eco.partition.abundances)[2:3]
  spp=size(eco.partition.abundances,1)
  birthrate=birth*timestep
  deathrate=death*timestep
  moverate=move*timestep
  for x in 1:dims[1]
    for y in 1:dims[2]
      square=eco.partition.abundances[:,x, y]
      for j in 1:spp
        if eco.traits.traits[j]!=eco.partition.habitat.matrix[x,y]
          birthrate=birthrate*0.8
          deathrate=deathrate
        end
        births= jbinom(1,j, birthrate)[1]
        deaths= jbinom(1,j, deathrate)[1]
        moves=jbinom(1, j, moverate)[1]

        neighbours=get_neighbours(eco.partition.habitat.matrix,y, x, 4)
        while moves > 0
          choose=sample(1:size(neighbours,1))
          destination=neighbours[choose,:]
          eco.partition.abundances[j, destination[1], destination[2]]=eco.partition.abundances[j, destination[1], destination[2]] + 1
          moves = moves - 1
        end
        eco.partition.abundances[j,x, y]= eco.partition.abundances[j,x, y] + births - deaths - moves
      end
    end
  end
end
