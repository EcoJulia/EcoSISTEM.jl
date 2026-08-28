# SPDX-License-Identifier: LGPL-3.0-or-later

module TestExtractclimate

using EcoSISTEM
# `[C7-VIS]` C: `extract_values` is `public` rather than exported — a debugging aid.
using EcoSISTEM: extract_values
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using DimensionalData
using DimensionalData: X, Y, Ti, dims
using DimensionalData.Lookups: Contains, SelectorError
using Rasters
using RasterDataSources
using Dates
using Test

include("rasterfixtures.jl")

# **THE ORACLE** (`[S3-DO]` item 3). What a query *should* return is derived here from the axes'
# own declarations — half-open `[lo, hi)` cells, 1-based months — and never from what the extractor
# currently does.
#
# **The fixtures must come from what a reader produces, not from the implementation's own
# conventions.** A fixture written the second way — an `AxisArray` with no sampling concept at all,
# and a time axis `(24000:24014) .* month_mean_duration` that is **0-based within the year**, which
# no reader has ever produced — keeps every assertion green and the
# defect alive. The identical defect was already found and fixed once in the plot recipes, whose
# comment records it: *"Only the test's hand-built 0-based axis, which no reader ever produces, kept
# that hidden."* This file was its surviving second instance.
#
# Every fixture below therefore comes from `_testraster`, which builds exactly what `datasetread.jl`
# builds — measured 2026-08-13 against real `WorldClim`/`EarthEnv` reads.

# **Non-square, and with a different step on each axis.** A square grid cannot see a Y/X
# transposition, and equal steps cannot see a single step being reused for both axes — which is a
# bug this very file once guarded against by hand (`_cellselectors`' separate `dlat`/`dlong`).
const LATS = (0.0:20.0:40.0) .* °          # 3 cells: [0,20) [20,40) [40,60)
const LONGS = (0.0:10.0:30.0) .* °         # 4 cells: [0,10) [10,20) [20,30) [30,40)

# Values that encode their own position, so a wrong cell is a wrong *number* rather than a plausible
# one: cell (i, j, k) holds `100i + 10j + k`.
_positional(ny, nx, nz) = [100i + 10j + k for i in 1:ny, j in 1:nx, k in 1:nz]

const MONTHLY = _testraster(WorldClim{Climate},
                            Float64.(_positional(3, 4, 12)) .* K,
                            lat = LATS, long = LONGS,
                            time = (1:12) .* month_mean_duration)
const BANDED = _testraster(WorldClim{BioClim},
                           Float64.(_positional(3, 4, 3)) .* K,
                           lat = LATS, long = LONGS,
                           code = [:bio1, :bio2, :bio3])
# A real **dated** axis, which `_stackdim` preserves from the source rather than rebuilding
# (`eltype(sourcedim.val) <: Dates.TimeType ? sourcedim : _mkstackaxis(...)`). Both `Date` and
# `DateTime` are exercised, and both are needed: a `year =` guard written as
# `eltype(axisvals) <: Unitful.Time` lets a genuinely dated axis fall straight through it.
const DATED = _testraster(WorldClim{Climate},
                          Float64.(_positional(3, 4, 15)) .* K,
                          lat = LATS, long = LONGS,
                          time = Date(2000, 1, 1) .+ Month.(0:14))

# **The fixture is `100i + 10j + k`** — latitude cell `i`, longitude cell `j`, time/band slice
# `k` — so every expected value below states which cell the query is supposed to land in, and a
# selector that quietly picks the neighbour produces a visibly wrong number rather than a
# plausible one. Cells are `Start`-labelled and half-open: a cell labelled `v` with step `s` spans
# `[v, v + s)`.
#
# There is deliberately **no oracle testset** asserting that contract against the axes themselves —
# half-open spans, 1-based months, codes-as-labels — before checking `extract_values` against it.
# Such a testset checks what `DimensionalData` and the fixture declare rather than anything this
# package does, which is more than a debugging helper warrants.

# **`extract_values` checked against the oracle above, never against itself.** Every expected
# value below is the cell the oracle says the query lands in, read out of the fixture by index — so a
# selector that quietly picks the neighbouring cell produces a wrong *number*, not a plausible one.
@testset "extract_values" begin
    @testset "one place, one slice, one value" begin
        # (25°, 15°) is inside lat cell 2 and long cell 2; March is index 3. `100i + 10j + k`.
        @test extract_values(MONTHLY, lat = 25.0°, long = 15.0°,
                             month = March) ==
              223.0K
        # A scalar request drops both spatial dims and the time dim, so the result is a number —
        # no `single` flag decides this, the selectors' own shapes do.
        @test extract_values(MONTHLY, lat = 25.0°, long = 15.0°,
                             month = March) isa Unitful.Quantity
    end

    @testset "the containing cell, including on a boundary" begin
        # Anywhere inside a cell gives that cell…
        @test extract_values(MONTHLY, lat = 39.999°, long = 19.999°,
                             month = January) ==
              extract_values(MONTHLY, lat = 20.0°, long = 10.0°,
                             month = January)
        # …and a coordinate exactly on a cell's *upper* edge belongs to the next cell, which is the
        # half-open contract the oracle pins.
        @test extract_values(MONTHLY, lat = 40.0°, long = 10.0°,
                             month = January) !=
              extract_values(MONTHLY, lat = 39.999°, long = 10.0°,
                             month = January)
        # Off the grid is an error, not the nearest cell. This is the whole of the `Contains`
        # decision (`[S3-F5]`) made visible at the public surface.
        @test_throws SelectorError extract_values(MONTHLY, lat = 90.0°,
                                                  long = 10.0°, month = January)
        @test_throws SelectorError extract_values(MONTHLY, lat = 25.0°,
                                                  long = -1.0°, month = January)
    end

    @testset "lat and long are required, and say so" begin
        @test_throws ErrorException extract_values(MONTHLY, long = 10.0°)
        @test_throws ErrorException extract_values(MONTHLY, lat = 10.0°)
        err = try
            extract_values(MONTHLY, long = 10.0°)
        catch e
            e
        end
        @test occursin("`lat`", err.msg) && occursin("extract_values", err.msg)
    end

    @testset "shape follows the request, and vectors CROSS" begin
        # Two latitudes × three longitudes × one month = 2 × 3.
        m = extract_values(MONTHLY, lat = [5.0°, 25.0°],
                           long = [5.0°, 15.0°, 25.0°], month = March)
        @test size(m) == (2, 3)
        # **Crossed, not zipped**: equal-length vectors give the full grid of combinations, not
        # the pairs. Zipping is an occurrence-record operation and is deliberately not offered.
        @test size(extract_values(MONTHLY, lat = [5.0°, 25.0°],
                                  long = [5.0°, 25.0°], month = March)) ==
              (2, 2)
        # One place, several slices → a vector along time only.
        v = extract_values(MONTHLY, lat = 25.0°, long = 15.0°,
                           month = [March, June])
        @test v isa AbstractVector && length(v) == 2
        # Omitting a keyword takes that dimension whole.
        @test length(extract_values(MONTHLY, lat = 25.0°, long = 15.0°)) == 12
    end

    @testset "a monthly climatology is 1-based" begin
        # **The contract at the public surface.** A month is 1-based: `(m - 1)` arithmetic returns
        # February when asked for March, which is what these assertions exist to catch.
        @test extract_values(MONTHLY, lat = 25.0°, long = 15.0°,
                             month = January) == 221.0K
        @test extract_values(MONTHLY, lat = 25.0°, long = 15.0°,
                             month = December) == 232.0K
        # A month outside 1…12 is refused by the shared 1-based conversion, in month names.
        @test_throws ErrorException extract_values(MONTHLY, lat = 25.0°,
                                                   long = 15.0°, month = 13)
    end

    @testset "a dated axis takes dates, intervals, and every-March" begin
        dates = parent(DimensionalData.lookup(DATED.array, Ti))
        @test extract_values(DATED, lat = 25.0°, long = 15.0°,
                             date = dates[3]) == 223.0K
        # `year = 2000` is gone; a calendar year is an interval, and composes.
        year2000 = extract_values(DATED, lat = 25.0°, long = 15.0°,
                                  date = Date(2000, 1, 1) .. Date(2000, 12, 31))
        @test length(year2000) == 12
        # **The case the old interface could not express at all**: every March, across years.
        @test length(extract_values(DATED, lat = 25.0°, long = 15.0°,
                                    month = March)) == 2
        # …and composed with a span, "March, but only in 2000" — one slice, but returned as a
        # one-element vector. **Deliberate, and the same rule the spatial axes follow**: a
        # dimension is dropped only when every keyword addressing it named exactly one thing, and
        # `date` here named a whole year. The shape follows the *request*, never the answer's size.
        @test extract_values(DATED, lat = 25.0°, long = 15.0°, month = March,
                             date = Date(2000, 1, 1) .. Date(2000, 12, 31)) ==
              [223.0K]
        # Keywords that select nothing in common say so rather than returning an empty answer.
        @test_throws ErrorException extract_values(DATED, lat = 25.0°,
                                                   long = 15.0°, month = March,
                                                   date = Date(2000, 6, 1) ..
                                                          Date(2000, 8, 1))
    end

    @testset "a layer axis is addressed by code" begin
        @test extract_values(BANDED, lat = 5.0°, long = 35.0°, code = :bio3) ==
              143.0K
        @test length(extract_values(BANDED, lat = 5.0°, long = 35.0°,
                                    code = [:bio1, :bio3])) == 2
        # A code is a label: an unknown one is refused, never read as a position.
        @test_throws SelectorError extract_values(BANDED, lat = 5.0°,
                                                  long = 35.0°, code = :bio9)
    end

    @testset "a keyword the dataset cannot answer is refused, by name" begin
        # **The error quality the keyword form buys.** The old `_sliceselector` managed only
        # "`year` only applies to a dataset with a time axis"; each of these says what this dataset
        # *does* have and which keyword reaches it.
        @test_throws ErrorException extract_values(BANDED, lat = 25.0°,
                                                   long = 15.0°, month = March)
        @test_throws ErrorException extract_values(MONTHLY, lat = 25.0°,
                                                   long = 15.0°, code = :bio1)
        @test_throws ErrorException extract_values(MONTHLY, lat = 25.0°,
                                                   long = 15.0°,
                                                   date = Date(2000, 1, 1))
        @test_throws ErrorException extract_values(DATED, lat = 25.0°,
                                                   long = 15.0°,
                                                   offset = 3month_mean_duration)
        code_err = try
            extract_values(BANDED, lat = 25.0°, long = 15.0°, month = March)
        catch e
            e
        end
        @test occursin("`code`", code_err.msg)

        # **`month` on a `Unitful.Time` axis that is not monthly** — the guard `[S3-F13]` asks
        # for, and the case a model output on a ~30-day timestep produces. It names `offset`.
        weekly = _testraster(WorldClim{Climate},
                             Float64.(_positional(3, 4, 5)) .* K,
                             lat = LATS, long = LONGS,
                             time = (1:5) .* EcoSISTEM.Units.week)
        @test_throws ErrorException extract_values(weekly, lat = 25.0°,
                                                   long = 15.0°, month = March)
        weekly_err = try
            extract_values(weekly, lat = 25.0°, long = 15.0°, month = March)
        catch e
            e
        end
        @test occursin("offset", weekly_err.msg)
        # …and `offset` addresses it by its own coordinates.
        @test extract_values(weekly, lat = 25.0°, long = 15.0°,
                             offset = 3EcoSISTEM.Units.week) == 223.0K
    end
end

end
