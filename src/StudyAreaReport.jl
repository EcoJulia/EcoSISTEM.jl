# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The report itself: how the grid was decided, what each layer contributed, what went wrong on the
# way, and the cache of reads that got there.

using DimensionalData

using Unitful

# What determines the bytes a source read produces: the dataset, the layer code and the read keywords.
# Deliberately not the spec object — two different specs asking for the same layer should share one
# read, a `ConstructedSpec` closes over a function that can never compare equal, and `SourceSpec`
# equality would stop holding the moment a `readkw` held a range or an array. Keying on the read also
# guarantees the cache holds nothing grid-dependent, since a read identity says nothing about any
# target grid.
struct ReadKey
    source::Type
    code::Union{CODE_TYPE, Vector{CODE_TYPE}}
    readkw::NamedTuple
end

"""
    LayerCache()

Hold the raw source reads made while deciding a [`StudyArea`](@ref), so that refining an area does
not read the same layer twice.

Only **raw** reads are cached, never anything reprojected onto a grid, which would be wrong to reuse
for a different area. Downloads are already cached on disk, so what this saves is decoding rather
than fetching.

# Fields

  - `reads`: the reads so far, keyed by what determines their content.
"""
struct LayerCache
    reads::Dict{ReadKey, Any}

    LayerCache() = new(Dict{ReadKey, Any}())
end

"""
    Problem

One thing worth telling the user about a study area.

# Fields

  - `severity`: a [`ProblemNotice`](@ref) or a [`ProblemWarning`](@ref).
  - `code`: names the kind of problem, so that it can be tested for and filtered on without matching
    against the prose.
  - `message`: the prose itself.
"""
struct Problem
    severity::AbstractProblemSeverity
    code::Symbol
    message::String

    function Problem(severity::AbstractProblemSeverity, code::Symbol,
                     message::AbstractString)
        return new(severity, code, String(message))
    end
end

"""
    LayerPlan

One layer's place in a study area: where it is, and what the chosen grid costs it.

# Fields

  - `name`: the layer's name.
  - `crs`, `step`, `bounds`: where the layer is, in its **own** coordinates. Stated that way so that
    the three agree with each other and are exact — an extent re-expressed in a distant CRS need be
    neither.
  - `kind`: what the target grid does to it, as an [`AbstractLayerFate`](@ref), which carries the
    aggregation factor or the reason for resampling along with it.
"""
struct LayerPlan
    name::Symbol
    crs::Any
    step::Any
    bounds::Any
    kind::AbstractLayerFate
end

"""
    StudyAreaReport

Everything decided about a study area, and what it costs.

Returned by [`investigate_study_area`](@ref) and held by every [`StudyArea`](@ref). The two share one
analysis, so a report can never describe a grid other than the one that would be built.

# Fields

  - `crs`, `cellsize`, `align`: the grid that was chosen.
  - `crssource`, `cellsizesource`: where each of those decisions came from — an
    [`AbstractDecisionSource`](@ref), [`GivenByUser`](@ref) if you supplied it, otherwise how it was
    derived. This is what lets the constructor announce everything it guessed.
  - `active`: **one array answering two questions**. Its coordinates say where the cells are — units,
    CRS, span, locus — and its values say which of them are simulated.
  - `simulate_safely`: the resolved rule `active` was decided by.
  - `layers`: a [`LayerPlan`](@ref) each, saying what the chosen grid costs that layer.
  - `footprint`: how much memory a run on this grid needs.
  - `problems`: the [`Problem`](@ref)s found on the way.
  - `specs`, `constraints`: enough to use this report as the base of another — the specs it was
    decided from, and the constraints **as given** rather than as resolved.
  - `cache`: the [`LayerCache`](@ref) of reads, or `nothing` on an as-built report, where it has been
    discarded.
  - `stage`: which of the two kinds of report this is — see [`AbstractReportStage`](@ref).
"""
struct StudyAreaReport
    crs::Any
    crssource::AbstractDecisionSource
    cellsize::Any
    cellsizesource::AbstractDecisionSource
    align::Union{Symbol, Nothing}
    # The eltype is not load-bearing, which is what lets one array carry both the geometry and the
    # mask: `_sampledata` hands this straight to `Rasters.resample(src, to = ...)`, which reads the
    # lookups and ignores the payload — `Float64`, `Matrix{Bool}` and `BitMatrix` targets all give
    # `isequal`-identical output. Rasters is fussy about the lookups themselves, `Regular` against
    # `Irregular`, which a `Bool` payload leaves untouched.
    #
    # Typed abstractly on purpose: a concrete `DimArray` type spells out every lookup parameter, and
    # the no-`Any` rule is about hot-path outputs, where a report never appears. The annotation is
    # here as a guard against storing a bare `Matrix`.
    active::DimensionalData.AbstractDimArray{Bool, 2}
    # The resolved rule rather than the constraint as given, so that the builder can apply the same
    # rule to a layer that was never named here, and the announcement can say which rule made the
    # grid.
    simulate_safely::Bool
    layers::Vector{LayerPlan}
    footprint::Any
    problems::Vector{Problem}
    # As given, not as resolved: re-deriving has to reach the same answer by the same route, and
    # recording a derived value as though it had been given would misreport its provenance.
    specs::NamedTuple
    constraints::NamedTuple
    # `nothing` means **discarded**, not "read nothing", and the `Union` exists to keep those two
    # apart: a synthetic area legitimately reports an empty cache, so emptiness alone could not say
    # which had happened. Discarded because the reads are consumed inputs — keeping them would pin
    # every raster a build touched for the life of the run, on every MPI rank, and nothing ever clears
    # the cache.
    cache::Union{LayerCache, Nothing}
    # Set internally at the two construction sites and by `_refinedreport`, and deliberately not a
    # constructor keyword: nothing outside the package should be able to claim a report is as-built.
    stage::AbstractReportStage
end

function Base.show(io::IO, c::LayerCache)
    return print(io,
                 "LayerCache($(length(c)) read$(length(c) == 1 ? "" : "s"))")
end

function Base.show(io::IO, ::MIME"text/plain", r::StudyAreaReport)
    ny, nx = Base.size(r.active)
    # The stage goes in the header because it changes what every number below means: on an `AsBuilt`
    # report `active` is what was built, while `specs` and `constraints` still say what was asked for.
    println(io, "StudyAreaReport ($(_stagephrase(r.stage)))")
    println(io,
            "  grid      $(ny) × $(nx) cells of $(_cellsizetext(r.cellsize)) ",
            "($(_sourcephrase(r.cellsizesource)))")
    println(io, "  crs       $(_crsname(r.crs)) ",
            "($(_sourcephrase(r.crssource)))")
    println(io, "  aligned   ",
            isnothing(r.align) ? "nothing — no layer is in the target CRS" :
            ":$(r.align)")
    println(io, "  active    $(count(r.active)) of $(length(r.active)) cells ",
            "($(round(100 * count(r.active) / length(r.active), digits = 1))%)")
    println(io, "  cells     ",
            r.simulate_safely ?
            "only those every layer covers completely (`simulate_safely`)" :
            "any whose centre has data, including partly covered ones (`simulate_safely = false`)")
    println(io, "  memory    $(_footprintmessage(r.footprint))")
    # Only worth a line when it is `nothing`: an empty cache is ordinary, since a synthetic area reads
    # nothing, but a discarded one is a fact about this report that nothing else records.
    isnothing(r.cache) &&
        println(io,
                "  reads     discarded — the layers were read when this was built, not kept")
    if !isempty(r.layers)
        println(io, "\n  layers")
        for p in r.layers
            println(io, _layerline(p))
        end
    end
    if isempty(r.problems)
        print(io, "\n  no problems found")
    else
        println(io, "\n  problems")
        for p in r.problems
            println(io,
                    "  [$(_severitytag(p.severity))] $(p.code): $(p.message)")
        end
    end
    return nothing
end

Base.show(io::IO, r::StudyAreaReport) = show(io, MIME"text/plain"(), r)

# ══ Functions ══════════════════════════════════════════════════════════════════════════════════

# How a severity is tagged in the report's own listing.
_severitytag(::ProblemNotice) = "info"

_severitytag(::ProblemWarning) = "warn"

# Raise a problem at its own severity — one method each, rather than a conditional on the value.
_report(::ProblemNotice, message) = @info message

_report(::ProblemWarning, message) = @warn message
