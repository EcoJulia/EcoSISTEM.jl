# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDiversitySet

using EcoSISTEM
using EcoSISTEM: gettimes
using EcoSISTEM.Units
using DataFrames
using Feather
using Unitful
using Unitful.DefaultSymbols
using Test

# `gettimes` has two branches, and only the first was ever exercised: with nothing recorded it hands
# back every requested time, and with a saved Feather file it hands back only those beyond the latest
# the file holds, so a part-finished run resumes rather than repeating itself. The second branch
# needs a file on disk, which is what the fixture below supplies.
@testset "gettimes without a saved file returns every requested time" begin
    mktempdir() do dir
        ds = EcoSISTEM.DiversitySet(missing, dir, [1.0year, 2.0year, 3.0year])
        @test gettimes(ds) == [1.0year, 2.0year, 3.0year]
    end
end

@testset "gettimes with a saved file resumes past the latest recorded time" begin
    mktempdir() do dir
        # The `time` column is stored bare, as `Feather` wrote it - the unit is reattached on read.
        Feather.write(joinpath(dir, "results.feather"),
                      DataFrame(time = ustrip.(s, [1.0year, 2.0year]),
                                measure = [0.5, 0.6]))
        ds = EcoSISTEM.DiversitySet(missing, dir,
                                    [1.0year, 2.0year, 3.0year, 4.0year])
        # Only the two beyond 2 years remain to be computed...
        @test gettimes(ds) == [3.0year, 4.0year]
        # ...and the loaded frame is cached on the set, in `Unitful.Time`.
        @test !ismissing(ds.data)
        @test eltype(ds.data[!, :time]) <: Unitful.Time
        # The column the reader adds for Diversity's benefit is present and blank.
        @test all(isempty, ds.data[!, :type_name])
        # Asking again uses the cached frame rather than re-reading, and answers the same.
        @test gettimes(ds) == [3.0year, 4.0year]
    end
end

end
