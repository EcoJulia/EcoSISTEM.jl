# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What an assembled model was built from: one question, asked of anything that holds a study area -
# its report, the area itself, a habitat, an ecosystem - each answering by walking down to the
# report, which holds the records of what was read, and adding what it knows itself: the grid, and
# for an ecosystem the run, as the `GridContext` of a `Provenance`.

"""
    provenance(report::StudyAreaReport)
    provenance(area::StudyArea)
    provenance(habitat::GridHabitat)
    provenance(eco::AbstractEcosystem)

Return a [`Provenance`](@ref): what this package version was, every input the object was built
from, and as its [`GridContext`](@ref) the grid it sits on and, for an ecosystem, the run. The
same answer comes from each object holding the same study area, and on every MPI rank, since it
is taken from the inputs alone.

An input is recorded by its provenance record where it has one - every file this package fetched -
and otherwise by its file name, with the catalogue's facts about its dataset where it is a layer of
one. The files read for the layers and those of the shape or mask that cut the grid are listed
alike; a synthetic layer reads nothing and adds nothing. An ecosystem adds the records given as
`provenance` to its species list and to [`build_ecosystem`](@ref), and those of every
[`Intervention`](@ref) that has acted on it.

# Arguments

  - `report`, `area`, `habitat`, `eco`: what to report on.
"""
function provenance(report::StudyAreaReport)
    within = get(report.constraints, :within, nothing)
    return Provenance(Provenancer.softwarerecord(EcoSISTEM),
                      _uniqueinputs(vcat(_reportinputs(report),
                                         _specinputs(within))),
                      GridContext(_gridrecord(report), nothing))
end

provenance(area::StudyArea) = provenance(area.report)

provenance(habitat::GridHabitat) = provenance(habitat.area)

function provenance(eco::AbstractEcosystem)
    built = provenance(eco.habitat)
    inputs = _uniqueinputs(vcat(built.inputs, eco.spplist.inputs, eco.inputs))
    return Provenance(built.software, inputs,
                      GridContext(built.context.grid, _runrecord(eco)))
end

"""
    write_provenance(path::AbstractString, obj)

Write what a study area, a habitat or an ecosystem was built from to `path` as TOML, and return the
path: the software, the grid, the run where there is one, and every input with its DOI, licence,
version, citation, file name and checksum - the file a paper's repository commits beside its
figures, so a reader can see what made them. The object's [`provenance`](@ref) is written as a
[`Provenance`](@ref) is.

An input is written as its record holds it, named by its file and never by a directory, so nothing
in the file is particular to the machine that wrote it. Quantities, the coordinate reference system
and the grid's extent are written as text, and fetch times and an epoch as TOML dates.

# Arguments

  - `path`: the file to write; its directory is created if it does not exist.
  - `obj`: a study area report, a `StudyArea`, a `GridHabitat` or an ecosystem, whose
    [`provenance`](@ref) is written, or a recorder, whose run's provenance as of its last write is.
"""
function write_provenance(path::AbstractString,
                          obj::Union{StudyAreaReport, StudyArea, GridHabitat,
                                     AbstractEcosystem})
    return write_provenance(path, provenance(obj))
end

function write_provenance(path::AbstractString, recorder::AbstractRecorder)
    record = provenance(recorder)
    isnothing(record) &&
        error("this `$(nameof(typeof(recorder)))` has not recorded anything yet, so there is no " *
              "run to write the provenance of.")
    return write_provenance(path, record)
end

# Where the cells are: the grid's coordinate reference system, its cell size, the extent its cells
# cover, and how many cells it has and how many are simulated.
function _gridrecord(report::StudyAreaReport)
    return (crs = report.crs, cellsize = report.cellsize,
            extent = Extents.extent(report.active),
            cells = length(report.active), active = count(report.active))
end

# What a run adds: the seed its streams are drawn from, the date its elapsed time counts from, and
# how far it has run.
function _runrecord(eco::AbstractEcosystem)
    return (seed = eco.seed, epoch = eco.epoch, calendar = eco.calendar,
            elapsed = eco.elapsed, placement = _placementrecord(eco.placement))
end

# How a run was placed, as its provenance says it: the placement as written, a function by its
# name, or `nothing` where the ecosystem records none.
_placementrecord(::Nothing) = nothing

_placementrecord(placement) = _placementphrase(placement)

# What a provenance listing adds about the run calendar: nothing for exact dates, and that months
# were counted as mean months otherwise.
_calendarnote(::ExactDates) = ""

_calendarnote(::MeanMonths) = ", months of 30.44 d"

# The files a `within` constraint names, which are read to cut the grid rather than through the
# cache of layer reads, each recorded as an outline: by its provenance record, or by its name.
function _specinputs(spec::RasterSpec)
    return _readinputs(spec, role = :region)
end

function _specinputs(spec::ConstructedRasterSpec)
    return reduce(vcat, (_specinputs(l) for l in spec.layers),
                  init = ProvenanceRecord[])
end

function _specinputs(spec::ShapeSpec; role::Symbol = :region)
    path = _localpath(spec.path)
    isfile(path) || return ProvenanceRecord[]
    return [something(_entryrecord(spec.path),
                      ProvenanceRecord(role = role, source = "file",
                                       path = basename(path)))]
end

function _specinputs(spec::NaturalEarthSpec; role::Symbol = :region)
    path = _localpath(_nesource(_checklevel(spec.level)))
    isfile(path) || return ProvenanceRecord[]
    record = something(Provenancer.recordat(path),
                       ProvenanceRecord(role = role, source = "file",
                                        path = basename(path)))
    return [_withcatalogue(record, datasetinfo(NaturalEarthLevel), path)]
end

function _specinputs(spec::ConstructedShapeSpec; role::Symbol = :region)
    return reduce(vcat, (_specinputs(m, role = role) for m in spec.members),
                  init = ProvenanceRecord[])
end

# Anything else a `within` or a combination member may be - `nothing`, a matrix, a box, a circle, a
# synthetic layer - names no file. Deliberately untyped, as the fallback for all of them.
_specinputs(::Any; role::Symbol = :region) = ProvenanceRecord[]

# A shape read as a **layer** rather than as a region: the same files, recorded under the layer's own
# role, so a habitat built from one says what ground it came from.
function _specinputs(spec::ShapeCoverage; role::Symbol = :habitat)
    return _specinputs(spec.shape, role = role)
end

# A shape wrapped in the rule deciding which cells it activates: the same files under the same role,
# since how much of a cell must be covered does not change where the ground came from. Without this
# the untyped fallback above answers "names no file", and the region vanishes from the provenance
# with nothing reporting it.
function _specinputs(spec::ShapeMaskSpec; role::Symbol = :region)
    return _specinputs(spec.shape, role = role)
end
