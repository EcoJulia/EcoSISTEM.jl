# SPDX-License-Identifier: LGPL-3.0-or-later
#
# THE EXTENSION BOUNDARY, grouped by extension so that "what does this extension owe us, and how do
# we know it is there?" is one place to look. Two kinds of thing live here:
#
#   - the HOOKS an extension supplies, each declared method-less. The parent states the generic and
#     its documentation; the sole method arrives with the package named in the section heading.
#   - the parent-side code that DETECTS an extension, at the end of the MPI section. It is real code
#     rather than a declaration, and it is here so that `Base.get_extension` is asked in one file.
#
# Two rules decide what that looks like:
#
#   - An extension ADDS a method, never overwrites one, so a method-less `function f end` here plus
#     the sole method in `ext/` is what avoids a precompilation warning.
#   - A docstring cannot live in an extension: `docs/src/api.md` is an `@autodocs` block over this
#     package's own modules and cannot see inside one. So every hook's documentation stays here, on
#     the stub, and describes what the extension actually implements. The same problem applies to a
#     type, which cannot be stubbed at all: the parent declares an abstract supertype carrying the
#     documentation and the extension supplies the concrete subtype under the same name.
#
# This mirrors `BaseInterface.jl`/`DiversityInterface.jl`/`EcoBaseInterface.jl`, which hold methods
# we supply on someone else's generic; here are the generics we declare for someone else to supply.
#
# Three things that look like they belong here and do not:
#
#   - `canonicalunit`, `bounds`, `supplytype` and `demandtype` (`nicheaxis_macro.jl`) are forward
#     declarations `@nicheaxis` needs, not extension hooks: it emits its methods fully qualified, and
#     a qualified definition needs the binding to exist already. Moving them breaks every axis.
#   - `_preferredcode` (`Climate.jl`) has its fallback method immediately below it, so declaration
#     and fallback stay together.
#   - `MPIEcosystem` and `MPIGridLandscape` are types, declared in `MPIEcosystem.jl` with the rest.
#     Only their hooks are here.

# ---------------------------------------------------------------------------
# EcoSISTEMPhyloExt - requires Phylo
# ---------------------------------------------------------------------------
# Everything to do with phylogenies. These are Phylo's concepts rather than EcoSISTEM's, which is why
# they are grouped here instead of sitting beside `SpeciesList` and the tolerances they were declared
# next to.

"""
    reroot!(tree::AbstractTree, node::String)

Reroot a phylogenetic tree by removing and recreating `node`, then attaching a
new root node `"NewRoot"` above the original root.
"""
function reroot! end

"""
    resettraits!(tree::AbstractTree)

Clear all node data records in the tree, resetting every node to an empty
`DataFrame`.
"""
function resettraits! end

"""
    assigntraits!(tree::AbstractTree, switch_rate::Vector{Float64}, traits::DataFrame)
    assigntraits!(tree::AbstractTree, switch_rate::Float64, traits::DataFrame)
    assigntraits!(tree::AbstractTree, traits::DataFrame)

Evolve functional traits through a phylogenetic tree, writing the result onto the
tree's node records.

The first two forms evolve **categorical** traits, switching between the values in
each column of `traits` at `switch_rate` - one rate per trait, or a single rate
shared by all of them. The third evolves **continuous** traits by Brownian motion,
and takes a `traits` with a `start` column of initial values and a `σ²` column of
rates, one row per trait.

# Arguments

  - `tree`: the phylogeny to evolve the traits along.
  - `switch_rate`: the rate of categorical trait change along a branch, per trait
    or shared by all of them. Omitted for the Brownian form.
  - `traits`: the traits to evolve - the categorical values per trait, or the
    `start`/`σ²` pairs the Brownian form needs.
"""
function assigntraits! end

"""
    gettraits(tree::AbstractTree, tips::Bool=true)

Retrieve functional traits assigned to a phylogenetic tree, either just tips or
all nodes.

"""
function gettraits end

"""
    discrete_evolve(numTraits::Int64, tree::BinaryTree)

Evolve a discrete switching trait along a BinaryTree, `tree`. Takes in a number
of traits, `numTraits` to be switched between and rate to switch between traits,
`switch_rate` with default value of 0.5.
"""
function discrete_evolve end

"""
    continuous_evolve(val::Union{Float64, Unitful.Quantity{Float64}}, var::Union{Float64, Unitful.Quantity{Float64}}, tree::BinaryTree)

Evolve a continuous trait along a BinaryTree, `tree` via Brownian motion. Takes
in a starting value, `val` and a variance, `var`.
"""
function continuous_evolve end

# Declared abstract because a struct cannot be stubbed the way a function can: the extension supplies
# the sole concrete subtype under the same name, as `MPIEcosystem` does. They must stay `public`
# rather than exported, or `using EcoSISTEM` inside the extension would import the bare names and the
# concrete definitions could not be written.
"""
    Brownian

A fitted Brownian-motion model of trait evolution, as returned by
[`fitbrownian`](@ref). Under it a trait wanders at random along every branch,
so expected divergence grows with the time two lineages have been apart.

# Fields of the concrete type `EcoSISTEMPhyloExt.Brownian`

  - `optimum`: the maximum-likelihood parameters, `[σ², z̄₀]` - the diffusion
    rate and the inferred root state.
  - `se`: their standard errors, from the Hessian; `show` prints each parameter
    with a ±2 SE interval.
  - `H`: the Hessian of the negative log-likelihood at the optimum.
  - `LL`: the log-likelihood there, for comparison against a competing model
    fitted to the same trait.
"""
abstract type Brownian end

"""
    varcovar(tree::AbstractTree)

Compute the phylogenetic variance-covariance matrix from the branch lengths of
`tree`. The diagonal entries are the root-to-tip distances and the off-diagonal
entries are the shared root-to-ancestor distances between pairs of tips.

# Arguments

  - `tree`: the phylogeny, whose **leaves** are the taxa the matrix is over. Its
    order fixes the row/column order, so a trait vector passed to
    [`fitbrownian`](@ref) must be in that same order.

Two tips sharing a longer ancestry share more of their history, so the matrix
is what turns "related species resemble each other" into a covariance a model
can be fitted against.
"""
function varcovar end

"""
    fitbrownian(tree::AbstractTree, traits::Vector{F}) where F <: AbstractFloat

Fit a Brownian motion model of trait evolution to `traits` measured at the tips
of `tree`, returning a [`Brownian`](@ref).

# Arguments

  - `tree`: the phylogeny relating the taxa. Its [`varcovar`](@ref) matrix is
    what carries the relatedness into the likelihood.
  - `traits`: one measured value per **tip**, in the tree's own leaf order. The two are matched by
    position and not by name, so a mismatch is silent rather than an error.
"""
function fitbrownian end

# The building blocks the trait-evolution functions above are written in terms of. Private to the
# package but declared here for the same reason as the rest: the methods are in the extension.
"""
    pair(vec)

Split a vector into consecutive overlapping pairs, returning a matrix with one
pair per row. Used to traverse root-to-tip paths branch by branch.
"""
function pair end

"""
    root_to_tips(tree)

Return all root-to-tip paths in a phylogenetic tree as a vector of node-name
vectors, one per tip.
"""
function root_to_tips end

"""
    arenoderecordsempty(tree::AbstractTree, nodes::Vector{String})

Check whether the node data records are empty for each node in `nodes`. Returns
a vector of `Bool` values, one per node.
"""
function arenoderecordsempty end

"""
    brownian_motion(T::Real, σ²::Float64, start::Float64, lab::String = "")

Evolve a Real value through Brownian motion, with a starting value,
 `start`, and rate, `σ²`.

"""
function brownian_motion end

# ---------------------------------------------------------------------------
# EcoSISTEMMPIExt - requires MPI
# ---------------------------------------------------------------------------
# The distributed hot loop's plumbing. The two abstract types it fills in are in
# `MPIEcosystem.jl`; what follows is what the extension must implement.

# Whether this process should build a distributed `MPIEcosystem`. The MPI extension supplies the sole
# method, `MPI.Initialized() && MPI.Comm_size(MPI.COMM_WORLD) > 1`; with the extension absent there is
# no method at all, which is why `_usempi` below checks for it before asking.
function _should_mpi end

"""
    gatherabundance(eco::MPIEcosystem)

Gather the full abundances matrix onto the root node.
"""
function gatherabundance end

"""
    gatherdiversity(eco::MPIEcosystem, divmeasure::F, q) where F <: Function

Gather diversity calculated by `divmeasure` at value `q` from all MPI nodes onto
the root node (rank 0), combining subcommunity diversity values using a power
mean weighted by total abundances across nodes.
"""
function gatherdiversity end

export gatherabundance, gatherdiversity

"""
    synchronise_from_rows!(mpigrid)

Copy the abundance matrix from its **species-partitioned** view into its
cell-partitioned one, so the next phase of the timestep sees the data in the
split it needs.

An [`MPIGridLandscape`](@ref) holds one matrix in two layouts because the
simulation needs different splits at different moments: demographics is per
species, dispersal is per cell. This is one half of moving between them.
"""
function synchronise_from_rows! end

"""
    synchronise_from_cols!(mpigrid)

Copy the abundance matrix from its **cell-partitioned** view back into the
species-partitioned one - the mirror of [`synchronise_from_rows!`](@ref).
"""
function synchronise_from_cols! end

# Every call site is inside `EcoSISTEMMPIExt` - these are hooks the extension implements, not names
# a user writes. Deliberately not symmetric with `MPIEcosystem` and `MPIGridLandscape`, which stay
# exported because a user names those types directly.
public empty_landscape, synchronise_from_rows!, synchronise_from_cols!

# --- Detecting it, on the parent's side ------------------------------------
# Everything below is implemented here rather than declared. It lives with the hooks because it is
# what asks whether they have arrived: both `Base.get_extension` calls in the package are here, and
# `build_ecosystem` reaches all of this through a single call to `_resolvedistributed`, naming MPI
# nowhere itself.

# Is a distributed run both possible and wanted? The extension must be loaded (so `_should_mpi` has a
# method at all) and that method must say yes.
function _usempi()
    return !isnothing(Base.get_extension(@__MODULE__, :EcoSISTEMMPIExt)) &&
           _should_mpi()
end

# Launcher env-vars set by Open MPI / MPICH / PMIx / Slurm - a cheap "was I started under an MPI
# launcher?" check that needs no MPI dependency. Used only to warn when a run looks distributed but
# MPI wasn't initialised (which would otherwise silently build one independent serial ecosystem per
# rank).
function _launchedundermpi()
    return any(k -> haskey(ENV, k),
               ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_RANK",
                "SLURM_STEP_NUM_TASKS"))
end

# Resolve the `distributed` kwarg of `build_ecosystem` to a Bool: `:auto` builds an `MPIEcosystem`
# only for a live multi-rank MPI session (`_usempi()`), `true` forces it (erroring if the MPI
# extension isn't loaded), `false` forces a serial `Ecosystem`.
function _resolvedistributed(distributed)
    distributed === false && return false
    if distributed === true
        isnothing(Base.get_extension(@__MODULE__, :EcoSISTEMMPIExt)) &&
            error("build_ecosystem(distributed = true) needs the MPI extension: run `using MPI` (and MPI.Init()) first.")
        return true
    end
    distributed === :auto ||
        error("`distributed` must be true, false or :auto; got $(repr(distributed)).")
    if _launchedundermpi() && !_usempi()
        @warn "Process looks MPI-launched but MPI is not initialised with >1 rank; building a serial Ecosystem. Call MPI.Init() before build_ecosystem for a distributed run."
    end
    return _usempi()
end

# ---------------------------------------------------------------------------
# EcoSISTEMERAExt - requires PyCall
# ---------------------------------------------------------------------------

"""
    retrieve_era5(args...; kwargs...)

Download ERA5 reanalysis data from the Copernicus Climate Data Store.

The download goes through the CDS Python client, so the method appears only once `PyCall` is loaded.

# Arguments

  - `args...`, `kwargs...`: forwarded unchanged to the extension's method, which documents them -
    load `PyCall` and consult `retrieve_era5` again to see the live signature.

The result is an ERA5 netCDF that [`readfile`](@ref)/[`ERA`](@ref) read like any other raster.
"""
function retrieve_era5 end

# `public` rather than exported because it is on its way out: the package will move to `CDSAPI.jl`,
# JuliaClimate's native interface to the same Copernicus service, which removes the `PyCall`
# dependency this hook exists to bridge. Supported until then.
public retrieve_era5

# ---------------------------------------------------------------------------
# EcoSISTEMRasterDataSourcesExt - requires RasterDataSources
# ---------------------------------------------------------------------------

"""
    compress_landcover(landcover::ClimateRaster{<:EarthEnv{<:LandCover}})

Collapse EarthEnv's twelve per-class cover fractions into a single layer of **class codes**, taking
each cell's dominant class.

# Arguments

  - `landcover`: the twelve bands as one multi-band [`ClimateRaster`](@ref), which is what
    `SourceSpec(EarthEnv{LandCover})` (no code) reads. The class code is the band's **position**, so
    all twelve must be present and in order; look one up by name with
    [`landcoverclass`](@ref EcoSISTEM.landcoverclass) rather than counting.

The result holds land-cover *types* (`1`-`12`, addressable by name with
[`landcoverclass`](@ref EcoSISTEM.landcoverclass)), which is a different quantity from its inputs:
each input band is a fraction of the cell covered by one class - a [`SurfaceArea`](@ref) - while the
output names which class won. Declaring the spec's `axis = LandCoverTypology` is what says so -
a `TypologyAxis` holds class labels - and it is why nothing downstream averages between class codes:
resampling takes the nearest class rather than interpolating, and a regime built from it is a
`CategoricalRegime`.

**Its niche axis is declared where it is used, not here.** A raster carries no axis of its own -
only a layer code the shipped catalogue can resolve - and no EarthEnv code means "typology", so say
so on the spec:

```julia
ConstructedSpec(EarthEnv{LandCover}, axis = LandCoverTypology) do lc
    compress_landcover(lc)
end
```

Prefer to let the combine run **after** the layers are sampled, which is the default: interpolating
the *percentages* and taking the argmax afterwards keeps each cell's sub-cell mix and cannot fabricate
a class, whereas collapsing first and resampling the codes invents classes that occur nowhere in the
data.
"""
function compress_landcover end

export compress_landcover

# The only public, non-deprecated name whose implementation lives in that extension, and so the one
# place a user meets the dependency without asking for climate data: `SetLandCover(:open_water)`
# reaches it, while `SetLandCover(7)` does not.
"""
    landcoverclass(name::Symbol)

The raw numeric code of the shipped `EarthEnv` land-cover class `name`, looked up by name in the
shipped table (never a hardcoded number). A building block for land-cover masks - e.g. an argmax of
all class bands, excluding open water:

```julia
ConstructedSpec(EarthEnv{LandCover}) do lc
    compress_landcover(lc) .!= landcoverclass(:open_water)
end
```

A raster broadcasts and yields a raster, so a combine names no array type: neither `.array` on the
way in nor a wrapper on the way out.
"""
function landcoverclass end

"""
    sourcecrs(T::Type{<:RasterDataSource}, layers = RasterDataSources.layers(T); kw...)

Return the coordinate reference system of a source's files without reading any of their data.

The file is opened lazily, so only its header is touched - which is what lets a [`StudyArea`](@ref)
settle its target CRS *before* deciding how much of each layer to read, rather than reading every
layer whole just to discover where it is. Returns `nothing` for a file that declares no CRS.

`cut`/`scale`/`fn` are accepted and ignored, so a [`SourceSpec`](@ref)'s stored read keywords can be
splatted in unchanged; the aggregation `scale` in particular cannot affect a CRS.

# Arguments

  - `T`: the `RasterDataSources` dataset type.
  - `layers`: which layers to inspect, defaulting to all of the dataset's. Only the first file that
    declares a CRS is needed, so this is rarely worth narrowing.
  - `kw...`: read keywords, all ignored - see above.
"""
function sourcecrs end

# ---------------------------------------------------------------------------
# EcoSISTEMDataPipelineExt - requires DataPipeline
# ---------------------------------------------------------------------------

# DataPipeline extension
"""
    EcoSISTEM.unziptemp(path::String)

Helper function for the FAIR Data Pipeline to unzip files that are stored as
zips to a temporary folder.
"""
function unziptemp end
