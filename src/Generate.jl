using Distributions
using DataFrames
function neighbours(x, y, dim)

function create_habitat(dims::Tuple, corr, prop)
  mat=zeros(10, 10)
  rand!(mat,[0 1])

end
P=DataFrame(Ind=assign[:,1],Species=assign[:,2],Loc=assign[:,3])
df[df[:A] % 2 .== 0, :]
P[P[:Loc] .== 1, :]
mat=ones(10, 10)
function jmulti(n::Int64, p::AbstractArray)
  rand(Multinomial(n, p))
end
function jmulti(n::Int64, p::Real)
  rand(Multinomial(n, repmat([p], n)))
end

dims=size(mat)
squares=dims[1]*dims[2]
assign=zeros(1000,3)
assign[:,1]=repmat([1],1000)
for i in 1:size(assign,1)
  assign[i,2]=rand(1.0:50,1)[1]
  assign[i,3]=rand(1:squares)
end
P=zeros(Int64,50,2,100)
map(a-> P[:,1,a]=collect(1:50), 1:100)
indvs=1000
while indvs>0
  sp=rand(1:50,1)
  loc=rand(1:100)
  P[sp,2,loc]=P[sp,2,loc]+1
  indvs=indvs-1
end
part=MatrixLandscape(P, Habitats(mat))


P=zeros(Int64,50,10,10)
abun_vec=rand(Multinomial(1000, 50))
for i in eachindex(abun_vec)
  abun=abun_vec[i]
    while abun>0
    x=rand(1:10)
    y=rand(1:10)
    P[i,x,y]=P[i,x,y]+1
    abun=abun-1
  end
end
part=MatrixLandscape(P, Habitats(mat))

for i in 1:100
 P.sub=P[P[:Loc] .== 1, :]
 =map(a->sum(assign[:,2].==a), collect(1.0:50))
end

P=map(a->sum(assign[:,2].==a), collect(1.0:50))

P=hcat(collect(1:50), P)


using Gadfly

display(plot(P, Geom.bar))


function populate(species::Int64, individuals::Int64, mat::AbstractArray, prob)
  dims=size(mat)
  pop=zeros(dims)
  # How many squares to fill?
  squares=dims[1]*dims[2]
  # What are the abundances of each species?
  abun_vec=jmulti(individuals, repmat([prob], species))
  for i in 1:length(abun_vec)
    abun=abun_vec[i]
    while abun>0
    random_ij= mapslices(rand,hcat(1:dims[1], 1:dims[2]), 1)
    pop[random_ij[1], random_ij[2]]=pop[random_ij[1], random_ij[2]]+1
    abun=abun-1
  end

end
