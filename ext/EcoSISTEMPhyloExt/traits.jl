# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Trees and traits -------------------------------------------------------------------------------
#
# The phylogenetic operations proper: rooting, assigning and reading traits, and evolving them
# along branches. Every one is the sole method for a hook the parent declares, and every docstring
# stays on that parent stub — `@autodocs` cannot see inside an extension.

function EcoSISTEM.reroot!(tree::AbstractTree, node::String)
    root = getroot(tree)
    deletenode!(tree, node)
    createnode!(tree, node)
    createnode!(tree, "NewRoot")
    createbranch!(tree, "NewRoot", root)
    return createbranch!(tree, "NewRoot", node)
end

function EcoSISTEM.resettraits!(tree::AbstractTree)
    nodes = getnodenames(tree)
    for i in nodes
        setnodedata!(tree, i, DataFrame())
    end
end

function EcoSISTEM.pair(vec)
    # Calc number of pairs
    npairs = length(vec) - 1
    # Create empty array
    newvec = Matrix{String}(undef, npairs, 2)
    # Split into pairs
    for i in collect(1:npairs)
        newvec[i, :] = vec[i:(i + 1)]
    end

    return newvec
end

function EcoSISTEM.root_to_tips(tree)
    tips = collect(nodenamefilter(isleaf, tree))
    paths = map(tips) do tps
        return reverse(nodehistory(tree, tps))
    end

    return paths
end

function EcoSISTEM.arenoderecordsempty(tree::AbstractTree,
                                       nodes::Vector{String})
    map(nodes) do nod
        return isempty(getnodedata(tree, nod))
    end
end

function EcoSISTEM.assigntraits!(tree::AbstractTree,
                                 switch_rate::Vector{Float64},
                                 traits::DataFrame)
    # Check if tree already assigned
    check = arenoderecordsempty(tree, collect(getnodenames(tree)))
    all(check) || error("Some nodes already assigned traits")

    # Calculate all branch paths from root to tips
    tips = collect(nodenamefilter(isleaf, tree))
    root = first(collect(nodenamefilter(isroot, tree)))

    paths = root_to_tips(tree)
    samp = DataFrame(hcat([rand(col) for col in eachcol(traits)]),
                     names(traits))
    # Assign first node a trait randomly
    setnodedata!(tree, root, samp)

    # Loop through all paths
    for i in paths
        # Split path into pairs of nodes
        pairs = pair(i)
        # Test if any branches have been assigned already
        unassigned = arenoderecordsempty(tree, i[i .!= root])
        pairs = pairs[unassigned, :]

        # Calculate how long the path is already
        len = distance(tree, first(pairs), last(pairs))
        # Calculate time to next switch
        # Draw switch times from exponential distribution
        # Stop when they are larger than the length of the path
        alltimes = map(switch_rate) do swt
            times = Array{Float64}(undef, 0)
            while (sum(times) < len)
                time_switch = rand(Exponential(swt * len))
                append!(times, time_switch)
            end
            return times
        end

        # Sum up the event times cumulatively
        cum_times = cumsum.(alltimes)

        # Run through the branches for the path, assigning a trait
        for j in axes(pairs, 1)
            sel_pair = pairs[j, :]
            # Calculate the length of the branch
            branch_len = distance(tree, first(sel_pair), last(sel_pair))

            # Calculate switch times
            num_switches = map(cum_times) do cum
                return sum(cum .< branch_len)
            end

            # Find trait of last node
            last_label = getnodedata(tree, first(sel_pair))

            map(num_switches) do number
                # If there are no switches, give same trait as previous node
                if number == 0
                    set_node = last(sel_pair)
                    setnodedata!(tree, set_node, last_label)
                else
                    # Else loop through for the required number of switches, sampling from list of traits
                    while number > 0
                        set_node = last(sel_pair)
                        newtrait = map(names(traits)) do trt
                            col = traits[!, trt]
                            return sample(col[col .!= last_label[:, trt]])
                        end
                        newtrait = DataFrame(hcat(newtrait), names(traits))
                        setnodedata!(tree, set_node, newtrait)
                        last_label = getnodedata(tree, set_node)
                        number = number - 1
                    end
                end
            end
        end
    end
end

function EcoSISTEM.assigntraits!(tree::AbstractTree, switch_rate::Float64,
                                 traits::DataFrame)
    return assigntraits!(tree, [switch_rate], traits)
end

function EcoSISTEM.brownian_motion(T::Real, σ²::Float64, start::Float64,
                                   lab::String = "")
    t = 0:T  # time
    # first, simulate a set of random deviates
    x = rand(Normal(0, sqrt(σ²)), length(t))
    # now compute their cumulative sum
    x = cumsum(append!([start], x))

    return x
end

function EcoSISTEM.assigntraits!(tree::AbstractTree, traits::DataFrame)

    # Warning for nodes that have already been assigned traits
    check = arenoderecordsempty(tree, collect(getnodenames(tree)))
    all(check) || @warn "Some nodes already assigned traits"

    # Check that the length of the starting values and variances are the same
    length(traits[!, :start]) == length(traits[!, :σ²]) ||
        error("Start values and variance must have same number of traits")

    # Find all names of nodes
    names = collect(getnodenames(tree))

    # Sort by distance from root
    root = first(collect(NodeNameIterator(tree, isroot)))
    dist = map(names) do node
        return distance(tree, root, node)
    end
    names = names[sortperm(dist)]

    # Loop through nodes in order of appearance
    for i in names
        # Set start value of the Brownian motion to root
        if isroot(tree, i)
            setnodedata!(tree, i, traits)
        else
            # Get value information for parent node
            pnt = getparent(tree, i)
            srt = getnodedata(tree, pnt)[!, :start]
            # Find length of path between parent and child node
            path = first(branchroute(tree, pnt, i))
            ln = getlength(tree, getbranch(tree, path))

            # Run the Brownian-motion model on each trait and set record
            newtrait = map(srt, traits[!, :σ²]) do start, sig
                return last(brownian_motion(ln, sig, start))
            end
            newdat = DataFrame(start = newtrait, σ² = traits[!, :σ²])
            setnodedata!(tree, i, newdat)
        end
    end
end

function EcoSISTEM.gettraits(tree::AbstractTree, tips::Bool = true)
    check = .!arenoderecordsempty(tree, getnodes(tree))
    all(check) || error("All node records empty")
    if tips
        nodes = nodenamefilter(isleaf, tree)
    else
        nodes = getnodenames(tree)
    end
    df = vcat(map(n -> getnodedata(tree, n), nodes)...)
    df[!, :species] .= nodes
    return df
end

function EcoSISTEM.discrete_evolve(numTraits::Int64, tree::BinaryTree,
                                   switch_rate = 0.5)
    # Create traits and assign to tips
    traits = DataFrame(trait1 = 1:numTraits)
    assigntraits!(tree, 0.5, traits)
    # Get traits from tree
    # `penalty = 0.5` pins the released behaviour — see the note in `species_list.jl`: this is a
    # niche label, so being outside it is a disadvantage rather than an exclusion.
    return SimpleCategoricalTolerance(Array(gettraits(tree, true)[:, 1]),
                                      axis = EcoSISTEM.TypologyAxis,
                                      penalty = 0.5)
end

function EcoSISTEM.continuous_evolve(val::Union{Float64,
                                                Unitful.Quantity{Float64}},
                                     var::Union{Float64,
                                                Unitful.Quantity{Float64}},
                                     tree::BinaryTree)
    # Create traits and assign to tips
    numspecies = length(getleafnames(tree))
    traits = DataFrame(start = ustrip(val), σ² = ustrip(var))
    assigntraits!(tree, traits)
    # Get traits from tree
    newtraits = gettraits(tree, true)
    newtraits[!, :start] = newtraits[!, :start] .* unit(val)
    return NicheTolerance(EcoSISTEM.NicheAxis, Normal, newtraits[!, :start],
                          newtraits[!, :σ²])
end
