# SPDX-License-Identifier: LGPL-3.0-or-later

# The tree operations a pool's phylogeny needs, the methods of the parent's hooks in
# `src/extensions.jl`. A pool keeps its own copy of the tree it is given, its tips named for its
# species, so each operation returns a new tree and leaves the one it was given alone.

function EcoSISTEM._tipnames(tree::AbstractTree)
    return String[tip for tip in getleafnames(tree)]
end

# Built afresh as a `RootedTree` of the same shape and branch lengths, since Phylo's `renamenode!`
# is optional and returns `false` on the `BinaryTree` that `rand(Ultrametric(...))` draws, and
# `droptips!` throws there too (Phylo 0.5.4: its node iterator declares `BinaryNode`s and yields
# names). The pool's tree is pruned and grafted later, so it stays a `RootedTree`. The kept tips are
# made first, under their new names, and every other node under a name of Phylo's choosing, so no
# two can clash; the rest are then pruned. Node data is not carried: the pool keeps a tree for its
# shape alone.
function EcoSISTEM._renamedtree(tree::AbstractTree, rename::AbstractDict)
    out = RootedTree()
    tips = getleafnames(tree)
    named = Dict{String, String}()
    for tip in tips
        haskey(rename, tip) &&
            (named[tip] = getnodename(out, createnode!(out, rename[tip])))
    end
    for node in getnodenames(tree)
        haskey(named, node) ||
            (named[node] = getnodename(out, createnode!(out)))
    end
    for branch in getbranches(tree)
        createbranch!(out, named[getnodename(tree, src(tree, branch))],
                      named[getnodename(tree, dst(tree, branch))],
                      getlength(tree, branch))
    end
    gone = [named[tip] for tip in tips if !haskey(rename, tip)]
    isempty(gone) || droptips!(out, gone)
    return out
end

function EcoSISTEM._prunedtree(tree::AbstractTree, names)
    out = deepcopy(tree)
    droptips!(out, collect(String, names))
    return out
end

# The new tip and the original hang from a new node where the original was, each on a branch of no
# length, so the two are as close as two tips can be and every other distance is unchanged.
function EcoSISTEM._graftedtree(tree::AbstractTree, original::AbstractString,
                                copy::AbstractString)
    out = deepcopy(tree)
    fork = createnode!(out)
    if hasinbound(out, original)
        branch = getinbound(out, original)
        parent = src(out, branch)
        len = getlength(out, branch)
        deletebranch!(out, branch)
        createbranch!(out, parent, fork, len)
    end
    createbranch!(out, fork, original, 0.0)
    createnode!(out, copy)
    createbranch!(out, fork, copy, 0.0)
    return out
end

function EcoSISTEM._phylotypes(tree::AbstractTree,
                               names::AbstractVector{<:AbstractString})
    return PhyloBranches(tree, names)
end
