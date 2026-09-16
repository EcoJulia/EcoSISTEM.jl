# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What an assembled model was built from: one question, asked of anything that holds a study area -
# its report, the area itself, a habitat, an ecosystem - each answering by walking down to the
# report, which holds the records of what was read, and adding what it knows itself - and the answer
# written out as the TOML file a paper's repository commits.

"""
    provenance(report::StudyAreaReport)
    provenance(area::StudyArea)
    provenance(habitat::GridHabitat)
    provenance(eco::AbstractEcosystem)

Return a [`Provenance`](@ref): what this package version was, every input the object was built
from, the grid it sits on and, for an ecosystem, the run. The same answer comes from each object
holding the same study area, and on every MPI rank, since it is taken from the inputs alone.

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
    return Provenance(_softwarerecord(),
                      _uniqueinputs(vcat(_reportinputs(report),
                                         _specinputs(within))),
                      _gridrecord(report), nothing)
end

provenance(area::StudyArea) = provenance(area.report)

provenance(habitat::GridHabitat) = provenance(habitat.area)

function provenance(eco::AbstractEcosystem)
    built = provenance(eco.habitat)
    inputs = _uniqueinputs(vcat(built.inputs, eco.spplist.inputs, eco.inputs))
    return Provenance(built.software, inputs, built.grid, _runrecord(eco))
end

"""
    write_provenance(path::AbstractString, p::Provenance)
    write_provenance(path::AbstractString, obj)

Write what a study area, a habitat or an ecosystem was built from to `path` as TOML, and return the
path: the software, the grid, the run where there is one, and every input with its DOI, licence,
version, citation, file name and checksum - the file a paper's repository commits beside its
figures, so a reader can see what made them.

An input is written as its record holds it, named by its file and never by a directory, so nothing
in the file is particular to the machine that wrote it. Quantities, the coordinate reference system
and the grid's extent are written as text, and fetch times and an epoch as TOML dates.

# Arguments

  - `path`: the file to write; its directory is created if it does not exist.
  - `p`: the [`Provenance`](@ref) to write.
  - `obj`: a study area report, a `StudyArea`, a `GridHabitat` or an ecosystem, whose
    [`provenance`](@ref) is written, or a recorder, whose run's provenance as of its last write is.
"""
function write_provenance(path::AbstractString, p::Provenance)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        return TOML.print(io, _provenancetoml(p), sorted = true)
    end
    return path
end

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

# This package as the software a result came from: its version, the DOI every version is cited by,
# and the Julia it ran on.
function _softwarerecord()
    return (package = "EcoSISTEM", version = string(pkgversion(@__MODULE__)),
            doi = _conceptdoi(), julia = string(VERSION))
end

# The Zenodo DOI the shipped `codemeta.json` declares, or `""` where the file or the DOI is missing.
# Read by pattern, since the file holds one and JSON is not a dependency.
function _conceptdoi()
    path = pkgdir(@__MODULE__, "codemeta.json")
    isfile(path) || return ""
    found = match(r"\"identifier\"\s*:\s*\"(10\.5281/zenodo\.[0-9]+)\"",
                  read(path, String))
    return isnothing(found) ? "" : String(only(found.captures))
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
            elapsed = eco.elapsed)
end

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
                  init = InputRecord[])
end

function _specinputs(spec::ShapeSpec; role::Symbol = :region)
    path = _localpath(spec.path)
    isfile(path) || return InputRecord[]
    return [something(_ourrecordat(path),
                      InputRecord(role = role, dataset = "file",
                                  path = basename(path)))]
end

function _specinputs(spec::NaturalEarthSpec; role::Symbol = :region)
    path = _localpath(_nesource(_checklevel(spec.level)))
    isfile(path) || return InputRecord[]
    record = something(_ourrecordat(path),
                       InputRecord(role = role, dataset = "file",
                                   path = basename(path)))
    return [_withcatalogue(record, datasetinfo(NaturalEarthLevel), path)]
end

function _specinputs(spec::ConstructedShapeSpec; role::Symbol = :region)
    return reduce(vcat, (_specinputs(m, role = role) for m in spec.members),
                  init = InputRecord[])
end

# Anything else a `within` or a combination member may be - `nothing`, a matrix, a box, a circle, a
# synthetic layer - names no file. Deliberately untyped, as the fallback for all of them.
_specinputs(::Any; role::Symbol = :region) = InputRecord[]

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

# A `Provenance` as the tables its TOML file holds: the software, the grid, the run where there is
# one, and a table per input.
function _provenancetoml(p::Provenance)
    toml = Dict{String, Any}("software" => _tomltable(p.software),
                             "grid" => _tomltable(p.grid),
                             "inputs" => Dict{String, Any}[_tomltable(r)
                                               for r in p.inputs])
    isnothing(p.run) || (toml["run"] = _tomltable(p.run))
    return toml
end

# A named tuple's or an input record's fields as a TOML table, leaving out those with nothing to say.
function _tomltable(nt::NamedTuple)
    return Dict{String, Any}(string(k) => _tomlvalue(v)
                             for (k, v) in pairs(nt) if !_tomlblank(v))
end

function _tomltable(r::InputRecord)
    return Dict{String, Any}(string(k) => _tomlvalue(getfield(r, k))
                             for k in fieldnames(InputRecord)
                             if !_tomlblank(getfield(r, k)))
end

# Whether a field has nothing to say: `nothing`, or empty text.
_tomlblank(v) = isnothing(v) || (v isa AbstractString && isempty(v))

# A value in a form TOML holds: text, booleans, floats, dates, integers that fit TOML's own, and
# tables and arrays of those. A seed too large for a signed integer is written as text.
_tomlvalue(v::AbstractString) = String(v)

_tomlvalue(v::Bool) = v

_tomlvalue(v::Union{Float64, Dates.Date, Dates.DateTime}) = v

function _tomlvalue(v::Integer)
    return typemin(Int64) <= v <= typemax(Int64) ? Int64(v) : string(v)
end

function _tomlvalue(v::AbstractDict)
    return Dict{String, Any}(string(k) => _tomlvalue(x) for (k, x) in v)
end

_tomlvalue(v::AbstractVector) = Any[_tomlvalue(x) for x in v]

# An extent as each dimension's lower and upper bound, as text.
function _tomlvalue(v::Extents.Extent)
    return Dict{String, Any}(string(k) => [string(first(b)), string(last(b))]
                             for (k, b) in pairs(Extents.bounds(v)))
end

# An EPSG code as it is written, `EPSG:27700`, rather than the type's own display.
_tomlvalue(v::Rasters.EPSG) = _crsname(v)

# Everything else - a quantity, another form of coordinate reference system, a symbol, another
# calendar's date - as its printed text. Deliberately untyped, as the fallback for all of them.
_tomlvalue(v) = string(v)
