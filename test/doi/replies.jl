# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The saved replies of `doi.org` in this directory, in the form `lookup_doi`'s `fetch` takes:
# `saved(url, type)` gives the reply for a DOI and a metadata type, or `nothing` where the type has
# no file, which stands for the 406 the resolver answers for a form a DOI does not have. Included by
# `test_SourceFacts.jl` and `test_assets.jl`.

const SAVED = @__DIR__

const FIXTURES = Dict("10.5285/ca339c86-3674-4030-b891-35326e71141e" => "eidc-envelopes",
                      "10.15468/dl.n4a2sa" => "gbif-download",
                      "10.5067/MODIS/MOD13C2.061" => "modis-mod13c2",
                      "10.24381/cds.f17050d7" => "era5-monthly",
                      "10.1111/geb.13564" => "geb-paper")

const FORMS = Dict("application/vnd.datacite.datacite+json" => "datacite",
                   "application/vnd.citationstyles.csl+json" => "csl")

function saved(url, type)
    doi = replace(url, "https://doi.org/" => "")
    path = joinpath(SAVED, FIXTURES[doi] * "." * FORMS[type] * ".json")
    return isfile(path) ? read(path, String) : nothing
end
