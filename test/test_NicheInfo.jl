# SPDX-License-Identifier: LGPL-3.0-or-later

module TestNicheInfo

using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Test

@testset "canonicalunit" begin
    @testset "SolarRadiationAxis" begin
        # WorldClim's srad is kJ·m⁻²·day⁻¹, CHELSA's rsds*/BioClimPlus's rsds_* are MJ·m⁻²·day⁻¹ — the
        # override reconciles both to one scale (see `_tocanon`/Layer.jl).
        @test EcoSISTEM.canonicalunit(SolarRadiation()) == kJ / (m^2 * day)
        @test EcoSISTEM.canonicalunit(SolarRadiationRange()) == kJ / (m^2 * day)
    end
end

@testset "every axis with disagreeing shipped units has a canonicalunit override" begin
    # Catalogue-level audit, generalising the `SolarRadiation` fix: for every concrete niche axis, if the
    # shipped sources (across all tables) disagree on unit, `canonicalunit` must reconcile them — otherwise
    # a regime built from one source and a tolerance defaulted from another can silently diverge in scale
    # (exactly the bug `SolarRadiation` had). Reuses the catalogue helpers that already back
    # `layerinfo`/`layerunit` — no new registry.
    for A in EcoSISTEM.ClimatePref._leafaxes()
        recs = EcoSISTEM.ClimatePref.layersbyaxis(A)
        isempty(recs) && continue
        units = unique(r.unit for r in recs)
        if length(units) > 1
            @test EcoSISTEM.canonicalunit(A()) !== nothing
        end
    end
end

end
