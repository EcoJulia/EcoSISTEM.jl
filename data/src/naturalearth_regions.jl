# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Regenerate `data/NaturalEarth/regions.csv`, the shipped table of named regions and their bounding
# boxes. Writes that file and nothing else.
#
# ## Why the table is shipped at all
#
# It is a **cache**, not a curated document: every value in it is derived, and the justification for
# committing it is that a bounding box should cost no download. `EcoSISTEM.boundingbox` therefore
# answers offline and instantly, while an actual *shape* still fetches the geometry.
#
# **So it must never be edited by hand.** Correct the generator, or the source, and re-run. A hand
# edit would be silently reverted by the next regeneration - which is the hazard
# `data/src/rate_table.jl` was fixed for, and the reason it is the shape copied here.
#
# ## Why the generator calls the package rather than reimplementing it
#
# Every geometric decision - which features a name covers, how they dissolve into components, how a
# component's area is measured, how the antimeridian is handled - comes from `EcoSISTEM` itself, via
# `_selectfeatures`, `_dissolve` and `_regionbox`. A generator with its own copy of that arithmetic
# would drift from the one a built shape uses, and the box would stop describing the shape - the one
# failure this table exists to prevent.
#
# ## Running it
#
# Downloads about 25 MB of Natural Earth zips on a cold cache, into `EcoSISTEM.assetdir`, and takes a
# few minutes. Not part of the test suite.
#
#     julia --project=data/src data/src/naturalearth_regions.jl

using EcoSISTEM
using EcoSISTEM: NATURALEARTH_LEVELS, _selectfeatures, _levelvalues, _dissolve,
                 _regionbox
using Unitful, Unitful.DefaultSymbols
using CSV

# No country-code column: `ISO_A3` is already a level, so a country is reachable by its code
# through the same lookup as by its name, and a code repeated per row would say nothing new.
#
# How many of the largest components each row describes individually. `LargestLandmass(count = n)`
# is answerable from the table for `n` up to this and needs the geometry beyond it.
const NPARTS = 5

# Coordinates are rounded OUTWARDS - west and south down, east and north up - never to nearest. A
# box rounded to nearest would sometimes no longer contain the region it describes, which is the one
# thing a bounding box has to do.
#
# Three decimal places, which is 111 m at the equator. What sets that is the use: a bounding box is
# passed to a raster reader as a cut, so the step only has to be small against a *cell*, and the
# finest data this package reads is 30 arcsec - about 928 m. An edge can therefore move by an eighth
# of a cell at worst, and no region in the table collapses to zero width. Anything needing the
# region's actual outline uses the geometry rather than the box.
_out(x, mode) = round(x, mode, digits = 3)
_lo(x) = _out(x, RoundDown)
_hi(x) = _out(x, RoundUp)

# Areas order the components and are reported, never used to reconstruct geometry, so they round to
# nearest at six significant figures - a square metre on the smallest region here.
_area(a) = round(ustrip(km^2, a), sigdigits = 6)

# One row per (level, name): the whole selection's box and area, then the largest `NPARTS`
# components individually so that a landmass query needs no download either.
function regionrow(level, name)
    parts = _dissolve(_selectfeatures(level, name))
    isempty(parts) && return nothing
    box = _regionbox(parts)
    row = Dict{Symbol, Any}(:Level => level.name, :Name => name,
                            :Parts => length(parts),
                            :Wraps => box.wraps,
                            :West => _lo(box.west), :South => _lo(box.south),
                            :East => _hi(box.east), :North => _hi(box.north),
                            :AreaKm2 => _area(sum(p -> p.area, parts)))
    for i in 1:NPARTS
        p = i <= length(parts) ? parts[i] : nothing
        row[Symbol("Part$(i)Area")] = isnothing(p) ? "" : _area(p.area)
        # A single component never crosses the date line - Natural Earth splits its polygons there -
        # so a part's box needs no wrap handling and is a plain envelope.
        row[Symbol("Part$(i)West")] = isnothing(p) ? "" : _lo(p.envelope.MinX)
        row[Symbol("Part$(i)South")] = isnothing(p) ? "" : _lo(p.envelope.MinY)
        row[Symbol("Part$(i)East")] = isnothing(p) ? "" : _hi(p.envelope.MaxX)
        row[Symbol("Part$(i)North")] = isnothing(p) ? "" : _hi(p.envelope.MaxY)
    end
    return row
end

const COLUMNS = vcat([:Level, :Name, :Parts, :Wraps, :West, :South, :East,
                         :North,
                         :AreaKm2],
                     [Symbol("Part$(i)$(f)")
                      for i in 1:NPARTS
                      for f in ("Area", "West", "South", "East", "North")])

rows = NamedTuple[]
for level in NATURALEARTH_LEVELS
    names = _levelvalues(level)
    t = @elapsed for name in names
        row = regionrow(level, name)
        isnothing(row) ||
            push!(rows,
                  NamedTuple{Tuple(COLUMNS)}(Tuple(row[c] for c in COLUMNS)))
    end
    println(rpad(level.name, 26), lpad(length(names), 5), " names  ",
            lpad(round(t, digits = 1), 7), " s")
end

const OUT = joinpath(pkgdir(EcoSISTEM), "data", "NaturalEarth", "regions.csv")
mkpath(dirname(OUT))
CSV.write(OUT, rows)
println("\nwrote ", length(rows), " rows to ", OUT)
println("size ", round(filesize(OUT) / 1024, digits = 1), " KB")
