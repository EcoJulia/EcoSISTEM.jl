```@meta
CurrentModule = EcoSISTEM
```

# Provenance: what a result was built from

A simulation is built from published data - climate layers, land cover, a region's outline,
occurrence records - and from the software that ran it. [`provenance`](@ref) answers, for a study
area, a habitat or an ecosystem, what that was: the package version and the DOI it is cited by, the
grid, the run, and a record of every input. [`write_provenance`](@ref) writes the same answer as a
TOML file, which is the file to commit beside a paper's figures so that a reader can see what made
them.

## What is recorded

Every file the package fetches gets a record written beside it - where it came from, when, its
size and checksum, whether it is the provider's own file or an extract cut to a request, and from
the catalogue (see [Layers](layers.md)) the file's title, the dataset's DOI, licence, version and
release date, and what its provider asks of anyone using the data: the works to cite and the
notices to carry, quoted as the provider wrote them and gathered on
[Providers' requests](requests.md). A date or year a provider leaves to be filled
in is filled from when the file was fetched wherever a record is shown or written. A listing ends by
saying what was done to the inputs - each was resampled or rasterised onto the study grid - which a
modified dataset's licence often asks to be said. When a study area reads a layer, it keeps a record of each file it read, and
those records stay with the area, and with a habitat or an ecosystem built on it, after the reads
themselves are discarded. The files of a shape or mask that cut the grid are recorded the same way.
A file with no record of its own - one you handed over rather than one the package fetched - is
recorded by its name, with the catalogue's facts about its dataset where it is a layer of one;
[`Provenancer.record_provenance`](@ref EcoSISTEM.Provenancer.record_provenance) gives it a record of its
own, below.

Each input is an [`ProvenanceRecord`](@ref), and a record never holds an absolute path, so the file can
be shared as it stands.

## Recording data the package did not fetch

Occurrence records, trait tables and parameters from the literature reach a model through your own
code, so the package cannot know where they came from. Record such a file once, where it is, with
[`Provenancer.record_provenance`](@ref EcoSISTEM.Provenancer.record_provenance): given the file's path, the
DOI of its source and the role the data play, it looks the DOI's metadata up, writes beside the
file the same `.provenance.toml` a fetch would, and returns the [`ProvenanceRecord`](@ref). A later
call from the same source under the same role finds the record and looks nothing up, so the call
can stay in a script that runs every time; one saying otherwise - another source, another role, or
a file the package fetched from its own catalogue - is refused.

```julia
using EcoSISTEM.Provenancer
envelopes = record_provenance("envelopes.csv", "10.5285/ca339c86-3674-4030-b891-35326e71141e",
                         role = :species)
```

Where a source has no DOI, or you know more than its metadata says, write the facts yourself as a
[`Provenancer.SourceFacts`](@ref EcoSISTEM.Provenancer.SourceFacts). `route` says how the file came
to be, since a download a service cut to your request is an extract rather than the published
dataset:

```@example provenance
using EcoSISTEM
using EcoSISTEM.Provenancer

path = joinpath(mktempdir(), "records_africa.csv")
write(path, "species,decimalLatitude,decimalLongitude\nAcacia tortilis,-3.1,36.7\n")
doi = "10.15468/dl.abc123"
facts = SourceFacts(source = "GBIF occurrence download", doi = doi,
                    publisher = "The Global Biodiversity Information Facility",
                    releasedate = "2026-09-12", licence = "cc-by-nc-4.0", identifies = :data,
                    citation = "GBIF.org (12 September 2026) GBIF Occurrence Download " *
                               "https://doi.org/" * doi)
seeding = record_provenance(path, facts, role = :abundance, route = ExtractedByService(doi))
```

[`provenance`](@ref) reads the record back from the file's path whenever it is wanted again. A
file in a place that cannot be written to takes a `sidecar` keyword naming where its record goes,
which `provenance(path)` takes too. The record says where the file and its sidecar are, so it
stands in for the file's path in a [`RasterFileSpec`](@ref) or a [`ShapeSpec`](@ref). Under the
FAIR data pipeline the call takes the run's `track` as well, and does the run's work for the file
before your code reads it; see [Data pipeline](pipeline.md).

Give the records to what is built from the data with the `provenance` keyword, which takes an
[`ProvenanceRecord`](@ref) or a vector of them: on [`build_species`](@ref) for the species, on
[`build_ecosystem`](@ref) for what belongs to the run as a whole, and on an [`Intervention`](@ref)
for data an intervention is built from - whose records join the ecosystem's once it has acted.

```@example provenance
species = build_species(DefaultEcosystem(), numspecies = 3, verbosity = :silent)
habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
eco = build_ecosystem(species, habitat, seed = 1, provenance = seeding)
provenance(eco)
```

This ecosystem is built on a synthetic grid, so it read no files: its only input is the record
given to it.

## Writing it out

```@example provenance
path = write_provenance(joinpath(mktempdir(), "provenance.toml"), eco)
print(read(path, String))
```

Blank fields are left out. Quantities, the coordinate reference system and the grid's extent are
written as text, and fetch times as TOML dates.

## What a record does not cover

A raster handed around on its own - one read with `read(spec)` and passed to a function - carries
no record; ask the spec or the study area it came from. A recorder - [`RecordAbundance`](@ref),
[`RecordDiversity`](@ref), [`SaveAbundance`](@ref) - keeps the run's provenance as it stood at its
last write, so `write_provenance(path, recorder)` writes it; `SaveAbundance` writes it beside every
file it saves. An array you fill yourself from a callback carries none, so write the ecosystem's
provenance beside it.

## The machinery, for another package

The records, the sidecar, the download that records itself and the DOI lookup are
`EcoSISTEM.Provenancer`'s, a submodule depending on nothing of EcoSISTEM's, so that another
package can record its own inputs the same way: it implements
[`Provenancer.IsFetchable`](@ref EcoSISTEM.Provenancer.IsFetchable) for each kind of thing it
fetches, saying which route the files come by; fetches a file from a URL with
[`Provenancer.fetchfile`](@ref EcoSISTEM.Provenancer.fetchfile), which records it as it lands; and
subtypes [`Provenancer.AbstractProvenanceContext`](@ref EcoSISTEM.Provenancer.AbstractProvenanceContext)
for what it adds to a result's account. Its names are on the [API](api.md) page.
