# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A made-up WFO Plant List release, small enough to write in a test: rows for each branch of the
# rules `WFOPlantList` resolves names by, under a classification from species to kingdom, written as
# the Darwin Core zip a release is, through GDAL's zip file system. Included by
# `test_NameAuthority.jl` and `test_SpeciesPool.jl`.

using EcoSISTEM: ArchGDAL

# The rows: identifier, name, authorship, rank, status, accepted identifier, parent. What each
# name exercises is in `FAKEWFO_EXPECTED` below.
const FAKEWFO_ROWS = [
    ("wfo-0000000400", "Plantae", "", "kingdom", "Accepted", "", ""),
    ("wfo-0000000300", "Fabales", "", "order", "Accepted", "",
     "wfo-0000000400"),
    ("wfo-0000000200", "Fabaceae", "", "family", "Accepted", "",
     "wfo-0000000300"),
    ("wfo-0000000100", "Acacia", "Mill.", "genus", "Accepted", "",
     "wfo-0000000200"),
    ("wfo-0000000101", "Mimosa", "L.", "genus", "Accepted", "",
     "wfo-0000000200"),
    ("wfo-0000000001", "Acacia alpha", "L.", "species", "Accepted", "",
     "wfo-0000000100"),
    ("wfo-0000000002", "Mimosa alpha", "L.", "species", "Synonym",
     "wfo-0000000001", ""),
    ("wfo-0000000003", "Acacia beta", "L.", "species", "Accepted", "",
     "wfo-0000000100"),
    ("wfo-0000000004", "Acacia beta", "Mill.", "species", "Synonym",
     "wfo-0000000005", ""),
    ("wfo-0000000005", "Acacia gamma", "L.", "species", "Accepted", "",
     "wfo-0000000100"),
    ("wfo-0000000006", "Acacia delta", "A", "species", "Accepted", "",
     "wfo-0000000100"),
    ("wfo-0000000007", "Acacia delta", "B", "species", "Accepted", "",
     "wfo-0000000101"),
    ("wfo-0000000008", "Acacia eps", "L.", "species", "Synonym",
     "wfo-0000000001", ""),
    ("wfo-0000000009", "Acacia eps", "Mill.", "species", "Synonym",
     "wfo-0000000005", ""),
    ("wfo-0000000010", "Acacia zeta", "L.", "species", "Accepted", "",
     "wfo-0000000100"),
    ("wfo-0000000011", "Acacia zeta", "Mill.", "species", "Synonym",
     "wfo-0000000010", ""),
    ("wfo-0000000012", "Acacia eta", "L.", "species", "Unchecked", "", ""),
    ("wfo-0000000013", "Acacia theta", "L.", "species", "Accepted", "",
     "wfo-0000000100"),
    ("wfo-0000000014", "Acacia theta", "Mill.", "species", "Unchecked", "",
     ""),
    ("wfo-0000000015", "Acacia iota", "L.", "species", "Unchecked", "", ""),
    ("wfo-0000000016", "Acacia iota", "Mill.", "species", "Unchecked", "", ""),
    ("wfo-0000000017", "Mimosa delta", "L.", "species", "Synonym",
     "wfo-0000000006", ""),
    ("wfo-0000000018", "Mimosa beta", "L.", "species", "Accepted", "",
     "wfo-0000000101"),
    ("wfo-0000000102", "Racosperma", "Mart.", "genus", "Synonym",
     "wfo-0000000100", ""),
    ("wfo-0000000103", "Acacia", "Mill.", "section", "Accepted", "",
     "wfo-0000000100")]

# What each name resolves to, `(name, id, match)`, for each pair of `allow_unchecked` and
# `prefer_accepted`; a name not recognised is `(nothing, nothing, reason)`. Accepted; a lone
# synonym; a homonym accepted once (`prefer_accepted` takes that row); a homonym accepted twice;
# synonyms of two species; rows all reaching one species; Unchecked alone, beside an Accepted row,
# and twice; a synonym of a homonym's species, which takes the homonym's joined key; a name WFO
# does not hold. "Mimosa beta" is a species of the second genus, "Racosperma" a synonym of the genus
# Acacia, and "Acacia" a section as well as a genus, for selecting species by rank.
const FAKEWFO_EXPECTED = let
    res(n, i, m) = (name = n, id = i, match = m)
    miss(m) = (name = nothing, id = nothing, match = m)
    common = Dict("Acacia alpha" =>
                      res("Acacia alpha", "wfo-0000000001", :accepted),
                  "Mimosa alpha" =>
                      res("Acacia alpha", "wfo-0000000001", :synonym),
                  "Acacia delta" => res("Acacia delta",
                      "wfo-0000000006+wfo-0000000007",
                      :homonym),
                  "Acacia eps" => res("Acacia eps",
                      "wfo-0000000008+wfo-0000000009",
                      :homonym),
                  "Acacia zeta" => res("Acacia zeta", "wfo-0000000010",
                      :accepted),
                  "Acacia theta" => res("Acacia theta", "wfo-0000000013",
                      :accepted),
                  "Mimosa delta" => res("Acacia delta",
                      "wfo-0000000006+wfo-0000000007",
                      :synonym),
                  "Algae x" => miss(:not_found))
    beta(pa) = res("Acacia beta",
                   pa ? "wfo-0000000003" : "wfo-0000000003+wfo-0000000004",
                   :homonym)
    Dict((au, pa) => merge(common,
               Dict("Acacia beta" => beta(pa),
                    "Acacia eta" =>
                        au ?
                        res("Acacia eta",
                            "wfo-0000000012",
                            :unchecked) :
                        miss(:unchecked_only),
                    "Acacia iota" =>
                        au ?
                        res("Acacia iota",
                            "wfo-0000000015+wfo-0000000016",
                            :homonym) :
                        miss(:unchecked_only)))
         for au in (false, true), pa in (false, true))
end

# Write the made-up release into `dir`, as the release's own file name, and return its path.
function fakewfo(dir::AbstractString)
    zip = joinpath(dir, "_DwC_backbone_R.zip")
    header = ["taxonID", "scientificName", "scientificNameAuthorship",
        "taxonRank", "taxonomicStatus", "acceptedNameUsageID",
        "parentNameUsageID", "family"]
    lines = [join(header, '\t');
             [join(vcat(["\"$x\"" for x in row], ["\"Fabaceae\""]), '\t')
              for row in FAKEWFO_ROWS]]
    bytes = Vector{UInt8}(join(lines, '\n') * '\n')
    handle = ArchGDAL.GDAL.vsifopenl("/vsizip/" * zip * "/classification.csv",
                                     "wb")
    ArchGDAL.GDAL.vsifwritel(bytes, 1, length(bytes), handle)
    ArchGDAL.GDAL.vsifclosel(handle)
    return zip
end

# Run `f` with the path of a made-up release in a temporary directory, removing the table cached
# from it afterwards, since that lives in the real asset cache.
function withfakewfo(f)
    return mktempdir() do dir
        zip = fakewfo(dir)
        try
            return f(zip)
        finally
            rm(EcoSISTEM._wfocachepath(zip), force = true)
        end
    end
end
