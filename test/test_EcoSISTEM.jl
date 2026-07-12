# SPDX-License-Identifier: LGPL-3.0-or-later

module TestEcoSISTEM

using EcoSISTEM
using Hwloc
using Test

@testset "species_blocksize" begin
    bs = EcoSISTEM.species_blocksize()
    # Used as an inner-loop block width in update! (cld(spp, bs), 1:bs:...), so it
    # must be a positive integer.
    @test bs isa Int
    @test bs >= 1
    # It is the detected CPU cache line divided by the abundance element size
    # (Int64), or the 128-byte-line fallback (16) if hwloc detection failed.
    expected = try
        max(1, Hwloc.cachelinesize() ÷ sizeof(Int64))
    catch
        16
    end
    @test bs == expected
end

end
