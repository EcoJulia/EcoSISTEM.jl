# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What this package adds to a result's provenance beyond the software and its inputs: the grid the
# result sits on and, for an ecosystem, the run; and the roles an input may enter the model by.

using Unitful

"""
    GridContext(grid::NamedTuple, run::Union{Nothing, NamedTuple})

The context of an EcoSISTEM result's [`Provenance`](@ref): where the cells are and, for an
ecosystem, the run. Built by [`provenance`](@ref) asked of a study area report, a study area, a
habitat or an ecosystem.

# Fields

  - `grid`: a named tuple of the grid's `crs` (`nothing` for a synthetic grid), its `cellsize`, the
    `extent` its cells cover, and how many `cells` it has and how many of them are `active`.
  - `run`: for an ecosystem, a named tuple of its `seed`, its `epoch` (`nothing` without one), the
    run `calendar` its dates were placed by, the time `elapsed`, and the `placement` its individuals
    were placed by, as written (`nothing` for one built from a population function); `nothing` for a
    study area or a habitat.
"""
struct GridContext <: AbstractProvenanceContext
    grid::NamedTuple
    run::Union{Nothing, NamedTuple}
end

function Base.show(io::IO, c::GridContext)
    return print(io, "GridContext(", c.grid.active, " of ", c.grid.cells,
                 " cells", isnothing(c.run) ? "" : ", seed $(c.run.seed)", ")")
end

function Provenancer.contextsummary(c::GridContext)
    return ", $(c.grid.active) of $(c.grid.cells) cells" *
           (isnothing(c.run) ? "" : ", seed $(c.run.seed)")
end

function Provenancer.contextlines(c::GridContext)
    lines = [rpad("grid", 10) *
             "$(c.grid.active) of $(c.grid.cells) cells active, cells of $(c.grid.cellsize)" *
             (isnothing(c.grid.crs) ? ", synthetic" :
              ", crs $(_crsname(c.grid.crs))")]
    isnothing(c.run) ||
        push!(lines,
              rpad("run", 10) * "seed $(c.run.seed), " *
              _elapsedname(c.run.elapsed) * " elapsed" *
              (isnothing(c.run.epoch) ? "" : " from $(c.run.epoch)") *
              _calendarnote(c.run.calendar))
    isnothing(c.run) || isnothing(c.run.placement) ||
        push!(lines, rpad("placement", 10) * c.run.placement)
    return lines
end

# The grid as its own table and the run as another, as the TOML file has always held them. The
# grid's extent is written as each dimension's lower and upper bound as text and its coordinate
# reference system as it is written, `EPSG:27700`, rather than the types' own display; the run's
# fields take TOML's forms where it has them and their printed text otherwise.
function Provenancer.contexttables(c::GridContext)
    tables = Dict{String, Any}("grid" => _tomlgrid(c.grid))
    isnothing(c.run) || (tables["run"] = Provenancer.tomltable(c.run))
    return tables
end

function _tomlgrid(grid::NamedTuple)
    table = Dict{String, Any}(string(k) => Provenancer.tomlvalue(v)
                              for (k, v) in pairs(grid)
                              if !(k in (:crs, :extent)) && !isnothing(v))
    isnothing(grid.crs) || (table["crs"] = _crsname(grid.crs))
    table["extent"] = Dict{String, Any}(string(k) => [string(first(b)),
                                            string(last(b))]
                                        for (k, b) in pairs(Extents.bounds(grid.extent)))
    return table
end

# What was done to the inputs to build on them, for the statement a modified dataset's licence or
# provider asks for: each was read for this grid and resampled or rasterised onto its cells.
function Provenancer.modification(c::GridContext)
    crs = isnothing(c.grid.crs) ? "" : " in $(_crsname(c.grid.crs))"
    return "contains modified data: each input was read for this grid's extent and resampled " *
           "or rasterised onto its cells of $(c.grid.cellsize)$(crs)."
end

# == Functions ======================================================================================

# The parts of the model an input can enter by, which a record given to this package must name.
const _INPUT_ROLES = (:habitat, :region, :species, :phylogeny, :dispersal,
                      :demography, :intervention, :abundance, :state,
                      :software)

# A `provenance` keyword's value as a vector of records, each with a role this package knows: one
# record, or a vector of them. Anything else is refused by name - the untyped fallback exists to
# give a better message than a `MethodError` naming this function.
_recordvector(record::ProvenanceRecord) = [_checkrole(record)]

function _recordvector(records::AbstractVector{ProvenanceRecord})
    return ProvenanceRecord[_checkrole(r) for r in records]
end

function _recordvector(x)
    return throw(ArgumentError("`provenance` takes an `ProvenanceRecord` or a vector of them; got a " *
                               "$(typeof(x))."))
end

# Refuse a record whose role is not one of the model's: a record takes any symbol, and this package
# says which it understands.
function _checkrole(record::ProvenanceRecord)
    record.role in _INPUT_ROLES ||
        throw(ArgumentError("`role` must be one of $(join(_INPUT_ROLES, ", ")); got " *
                            "`$(record.role)`."))
    return record
end

# A run's elapsed time in the unit it is most easily read in: seconds under a day, days under a
# year, and years beyond.
function _elapsedname(t::Unitful.Time)
    t < 1.0u"d" && return string(round(typeof(1.0u"s"), t, digits = 1))
    t < 1.0u"yr" && return string(round(typeof(1.0u"d"), t, digits = 1))
    return string(round(typeof(1.0u"yr"), t, digits = 2))
end
