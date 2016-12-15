# Function to randomly populate a habitat matrix
function populate(species::Int64, individuals::Int64, habitat::Habitats)
  # Calculate size of habitat
  dim=size(habitat.matrix)
  # Create empty population matrix of correct dimensions
  P=zeros(Int64,species,dim[1],dim[2])
  # Randomly choose abundances for each species from Multinomial
  abun_vec=rand(Multinomial(individuals, species))
  # Loop through species
  for i in eachindex(abun_vec)
    # Get abundance of species
    abun=abun_vec[i]
      # Loop through individuals
      while abun>0
      # Randomly choose position on grid
      x=rand(1:dim[1])
      y=rand(1:dim[1])
      # Add individual to this location
      P[i,x,y]=P[i,x,y]+1
      abun=abun-1
    end
  end
  # Create MatrixLandscape from P matrix and habitat
  MatrixLandscape(P, habitat)
end

# Function to calculate species richness from an Ecosystem object
function SR(ecosystem::Ecosystem)
  sz=size(ecosystem.partition.abundances,2,3)
  ms=map(x-> sum(ecosystem.partition.abundances[:,x].>0), 1:(sz[1]*sz[2]))#*individuals
  reshape(ms, sz)
end

# Function to create a habitat from a discrete set of types
function create_habitat(dim, types, prop)
  # Weighted sample from the types in the correct dimension
  sample(types, WeightVec(prop), dim)
end



# Function to populate a Niche habitat
function populate(species::Int64, individuals::Int64, habitat::Niches, traits::Vector)
  # Calculate size of habitat
  dim=size(habitat.matrix)
  # Create empty population matrix of correct dimensions
  P=zeros(Int64,species,dim[1],dim[2])
  # Randomly choose abundances for each species from Multinomial
  abun_vec=rand(Multinomial(individuals, species))
  # Loop through species
  for i in eachindex(abun_vec)
    # Get abundance of species
    abun=abun_vec[i]
    # Get species trait
    pref=traits[i]
    # Calculate weighting, giving preference to squares that match with trait
    wv= Vector{Float64}(100)
    wv[find(reshape(habitat.matrix, (100,1)).==pref)]= 0.9
    wv[find(reshape(habitat.matrix, (100,1)).!=pref)]= 0.1
    # Loop through individuals
      while abun>0
        # Randomly choose position on grid (weighted)
      pos=sample(1:(dim[1]*dim[2]), weights(wv))
      # Add individual to this location
      P[i,pos]=P[i,pos]+1
      abun=abun-1
    end
  end
  # Create MatrixLandscape from P matrix and habitat
  MatrixLandscape(P, habitat)
end


# Function to get the neighbours of a grid square in a matrix in 4 or 8 directions
function get_neighbours(mat::Matrix, y::Int64, x::Int64, chess::Int64=4)
  # Calculate dimensions
  dims=size(mat)
  # Include 4 directions
  if chess==4
    neighbour_vec=[x y-1; x y+1; x-1 y; x+1 y]
  # Include 8 directions
  elseif chess==8
    neighbour_vec=[x y-1; x y+1; x-1 y; x+1 y; x-1 y-1; x-1 y+1; x+1 y-1; x+1 y+1]
  else
    # Give error if other number chosen than 4 or 8
    error("Can only calculate neighbours in 4 or 8 directions")
  end
  # Remove answers outside of the dimensions of the matrix
  remove=vcat(mapslices(all, [neighbour_vec.>=1 neighbour_vec[:,1].<=dims[1] neighbour_vec[:,2].<=dims[2]], 2)...)
  neighbour_vec=neighbour_vec[remove,:]
  neighbour_vec
end

# Function to update a Ecosystem after one timestep- stochastic birth, death and movement
function update!(eco::Ecosystem,  birth::Float64, death::Float64, move::Float64, timestep::Real)
  # Calculate abundance in overall grid (to be implemented later?)
  #abun=map(i->sum(eco.partition.abundances[i,:,:]), 1:size(eco.partition.abundances,1))
  # Calculate dimenions of habitat and number of species
  dims=size(eco.partition.abundances)[2:3]
  spp=size(eco.partition.abundances,1)

  # Calculate effective rates
  birthrate=birth*timestep
  deathrate=death*timestep
  moverate=move*timestep

  # Loop through grid squares
  for x in 1:dims[1]
    for y in 1:dims[2]

      square=eco.partition.abundances[:,x, y]

      # Loop through species in chosen square
      for j in 1:spp
        # If zero abundance then go extinct
        if square[j] == 0
          birthrate = 0
          deathrate = 0
          moverate = 0
        # If traits are same as habitat type then give birth "boost"
        elseif eco.traits.traits[j]!=eco.partition.habitat.matrix[x,y]
          birthrate=birthrate*0.8
          deathrate=deathrate
        end
        # Throw error if rates exceed 1
        birthrate<=1 && deathrate<=1 && moverate<=1 || error("rates larger than one in binomial draw")

        # Calculate births, deaths and movements
        births= jbinom(1,square[j], birthrate)[1]
        deaths= jbinom(1,square[j], deathrate)[1]
        moves=jbinom(1, square[j], moverate)[1]

        # Find neighbours of grid square
        neighbours=get_neighbours(eco.partition.habitat.matrix,y, x, 8)
        # Loop through number of movements
        while moves > 0
          # Randomly sample one of the neighbours
          choose=sample(1:size(neighbours,1))
          destination=neighbours[choose,:]
          # Add one to this neighbour
          eco.partition.abundances[j, destination[1], destination[2]]=eco.partition.abundances[j, destination[1], destination[2]] + 1
          moves = moves - 1
        end
        # Update population
        eco.partition.abundances[j,x, y]= eco.partition.abundances[j,x, y] + births - deaths - moves
      end
    end
  end
end
