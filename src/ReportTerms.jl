# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The rest of a study-area report's vocabulary: what became of each layer, how far the report has
# got, and how loudly a problem with it should speak.

using Unitful

"""
    AbstractLayerFate

What the chosen study grid costs one layer: [`LayerKeptExactly`](@ref), [`LayerAggregated`](@ref) or
[`LayerResampled`](@ref). Each carries exactly the detail its own case has - the aggregation factor, or
the reason for resampling - so a field can never be set for a case it means nothing to.
"""
abstract type AbstractLayerFate end

"""    LayerKeptExactly <: AbstractLayerFate - the layer is copied onto the grid, cell for cell. """
struct LayerKeptExactly <: AbstractLayerFate end

"""
    LayerAggregated(factor)

The layer is combined `factor` source cells to a target cell - still exact, no interpolation.
"""
struct LayerAggregated <: AbstractLayerFate
    factor::Int
end

"""
    LayerResampled(reason)

The layer is interpolated onto the grid, introducing uncertainty - the thing the whole alignment
story exists to avoid. `reason` says why it was unavoidable.
"""
struct LayerResampled <: AbstractLayerFate
    reason::String
    LayerResampled(reason::AbstractString) = new(String(reason))
end

"""
    AbstractReportStage

Which of the two things a [`StudyAreaReport`](@ref) describes: an area as **proposed**
([`AsInvestigated`](@ref)) or as **built** ([`AsBuilt`](@ref)). Both are the same type with the same
fields, so this is what says which question the numbers answer.

It matters because building can narrow the grid: a report's `specs` and `constraints` always say
what was *asked for*, while a [`GridHabitat`](@ref)'s layers may have removed cells the study area
still lists.

It is therefore also what makes a built report safe to reuse. Any value carrying a report can be a
[`StudyArea`](@ref)'s `base`, and re-deriving from one would rebuild the *proposed* grid and discard
that narrowing - so an `AsBuilt` base given no other keyword is copied verbatim instead, which is how
[`build_habitat`](@ref) reseeds a habitat from the grid it is actually on.
"""
abstract type AbstractReportStage end

"""    AsInvestigated <: AbstractReportStage - an area as proposed: nothing has been built on it yet. """
struct AsInvestigated <: AbstractReportStage end

"""    AsBuilt <: AbstractReportStage - an area as built: this describes a `GridHabitat` that exists. """
struct AsBuilt <: AbstractReportStage end

"""
    AbstractProblemSeverity

How much a [`Problem`](@ref) matters: a [`ProblemNotice`](@ref) or a [`ProblemWarning`](@ref). A
type rather than a symbol, so an unrecognised severity is refused by the signature where it is
written rather than by a check inside the constructor.
"""
abstract type AbstractProblemSeverity end

"""    ProblemNotice <: AbstractProblemSeverity - something was guessed, or is lossy but expected. """
struct ProblemNotice <: AbstractProblemSeverity end

"""    ProblemWarning <: AbstractProblemSeverity - it will work, but is probably not what you want. """
struct ProblemWarning <: AbstractProblemSeverity end

# --- Saying it in prose ------------------------------------------------------
# How each term reads in a report a person is shown, given the term itself.

# A layer's fate, as it appears in the report's per-layer listing.
_fatephrase(::LayerKeptExactly) = "kept exactly"

_fatephrase(fate::LayerAggregated) = "aggregated $(fate.factor)× (exact)"

_fatephrase(fate::LayerResampled) = "RESAMPLED - $(fate.reason)"

# A report's stage, as it appears in its header line.
_stagephrase(::AsInvestigated) = "as investigated"

_stagephrase(::AsBuilt) = "as built"
