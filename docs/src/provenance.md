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
size and checksum, and the dataset's DOI, licence, version and citation from the catalogue (see
[Layers](layers.md)). When a study area reads a layer, it keeps a record of each file it read, and
those records stay with the area, and with a habitat or an ecosystem built on it, after the reads
themselves are discarded. The files of a shape or mask that cut the grid are recorded the same way.
A file with no record of its own - one you handed over rather than one the package fetched - is
recorded by its name, with the catalogue's facts about its dataset where it is a layer of one.

Each input is an [`InputRecord`](@ref), and a record never holds an absolute path, so the file can
be shared as it stands.

## Recording data the package did not fetch

Occurrence records, trait tables and parameters from the literature reach a model through your own
code, so the package cannot know where they came from. Say so with the `provenance` keyword, which
takes an [`InputRecord`](@ref) or a vector of them: on [`build_species`](@ref) for the species, on
[`build_ecosystem`](@ref) for what belongs to the run as a whole, and on an [`Intervention`](@ref)
for data an intervention is built from - whose records join the ecosystem's once it has acted.

```@example provenance
using EcoSISTEM

seeding = EcoSISTEM.InputRecord(role = :abundance,
                                dataset = "GBIF occurrence download",
                                doi = "10.15468/dl.abc123",
                                licence = "CC BY-NC 4.0")
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
