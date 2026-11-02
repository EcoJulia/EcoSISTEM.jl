# Saved replies of doi.org

The metadata `test_SourceFacts.jl` reads the facts of a source from, as `https://doi.org/<doi>`
returns it with an `Accept` header, so that the readers are tested with no network:

| file | DOI | form |
|---|---|---|
| `eidc-envelopes.datacite.json` | `10.5285/ca339c86-3674-4030-b891-35326e71141e` | DataCite, as returned |
| `gbif-download.datacite.json` | `10.15468/dl.n4a2sa` | DataCite, cut to the fields the reader uses and three of its 1,678 related identifiers |
| `modis-mod13c2.datacite.json` | `10.5067/MODIS/MOD13C2.061` | DataCite, as returned |
| `era5-monthly.datacite.json` | `10.24381/cds.f17050d7` | DataCite, as returned |
| `geb-paper.csl.json` | `10.1111/geb.13564` | CSL JSON, the paper's reference list and Crossref's indexing fields cut; DataCite answers 406 for it |

A test names a file by the DOI and the form, and a form with no file stands for a 406.
