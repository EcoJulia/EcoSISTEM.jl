# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Run the data pipeline example (`examples/pipeline/`) and the paper's small figures
# (`examples/paper/paper.yaml`) as `fair` runs against a registry of this test's own, and check what
# the registry records of them. The runs load this package from the checkout under the test
# process's own Julia flags, so they count towards its coverage - the pipeline extension's included.
#
# Only through the test harness, which provides DataPipeline, Git and Sockets:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_pipeline.jl"])'
#
# It runs when named, as above, or with `ECOSISTEM_TEST_PIPELINE=true`, which the CI job that
# requires it sets; otherwise it says so and skips. It needs `fair` - `fair-cli` from PyPI, or
# DataPipeline's own copy, `DataPipeline.install_cli!()` - a POSIX shell and the network.
#
# Its pipeline is kept in EcoSISTEM's scratch space, so nothing of the user's is touched: a home for
# every `fair` command and what it runs, holding a registry installed there the first time (a few
# minutes) and emptied before each run of this file, on a free port rather than 8000; and a copy of
# the repository's tracked files, made afresh each time, to `fair init` in.

module ExtrasPipeline

using Test
using EcoSISTEM
using DataPipeline
using Git
using Scratch: get_scratch!
using Sockets

const D = DataPipeline

# The registry release installed, the one `.github/workflows/docs.yaml` installs.
const REGISTRY_VERSION = "v1.4.0"

const REQUIRED = get(ENV, "ECOSISTEM_TEST_PIPELINE", "") == "true"
const ASKED_FOR = any(a -> occursin("extras_pipeline", a), ARGS)

# A word quoted for a POSIX shell.
shquote(word) = "'" * replace(word, "'" => "'\\''") * "'"

# Write a `julia` for the runs to find first on the path: this process's own binary and flags,
# coverage included, with each run's arguments after them, so that a run's own `--threads` and
# `--project` win. Return the folder it is in.
function juliawrapper(dir)
    bin = mkpath(joinpath(dir, "bin"))
    path = joinpath(bin, "julia")
    write(path,
          "#!/bin/sh\nexec " * join(shquote.(Base.julia_cmd().exec), " ") *
          " \"\$@\"\n")
    chmod(path, 0o755)
    return bin
end

# The repository's tracked files as they stand, the examples' files not yet staged, so that a new
# example runs here before it is, and its `.git`, copied into `repo`.
function copyrepository(repo)
    rm(repo, recursive = true, force = true)
    root = pkgdir(EcoSISTEM)
    files = readlines(`$(Git.git()) -C $root ls-files`)
    append!(files,
            readlines(`$(Git.git()) -C $root ls-files --others --exclude-standard -- examples`))
    for file in files
        source = joinpath(root, file)
        isfile(source) || continue
        mkpath(dirname(joinpath(repo, file)))
        cp(source, joinpath(repo, file))
    end
    cp(joinpath(root, ".git"), joinpath(repo, ".git"))
    return repo
end

# The CLI configuration `fair init --ci` writes, but for the port, which `--ci` fixes at 8000, where
# a registry of the user's own may be running.
function cliconfig(path; home, repo, registry, port)
    remote = try
        readchomp(`$(Git.git()) -C $(pkgdir(EcoSISTEM)) remote get-url origin`)
    catch
        "https://github.com/EcoJulia/EcoSISTEM.jl"
    end
    store(folder) = joinpath(home, folder) * "/"
    token = joinpath(registry, "token")
    write(path, """
          namespaces:
            input: testing
            output: testing
          registries:
            local:
              uri: http://127.0.0.1:$port/api/
              directory: $registry
              data_store: $(store("data_store"))
              token: $token
            origin:
              uri: http://127.0.0.1:$(port + 1)/api/
              data_store: $(store("origin"))
              token: $token
          user:
            email: test@noreply.com
            family_name: Test
            given_names: Interface
            orcid: 000-0000-0000-0000
            uri: https://github.com/FAIRDataPipeline
            uuid: 2ddb2358-84bf-43ff-b2aa-3ac7dc3b49f1
            remote_user: FAIRDataPipeline
          git:
            local_repo: $repo
            remote: origin
            remote_repo: $remote
          """)
    return path
end

# Every code run the registry holds from the copy's commits, in the order they ran.
function coderuns(registry, repo)
    found = D.identify(registry, D.AllCommits(repo = repo, dirty = true))
    return sort!([run for commit in found.commits for run in commit.runs],
                 by = run -> run.run_date)
end

# A data product as the pipeline names one, without its version.
named(product) = product.namespace * ":" * product.name

if !(REQUIRED || ASKED_FOR)
    @info "Skipping the pipeline runs: name `extras_pipeline` in `test_args`, or set " *
          "`ECOSISTEM_TEST_PIPELINE=true`, to run them."
elseif Sys.iswindows()
    error("the pipeline runs need a POSIX shell for their `julia`, and this is Windows.")
else
    @testset "Pipeline runs" begin
        space = get_scratch!(EcoSISTEM, "pipeline-tests")
        # The copy is within the home, where the CLI's search for a project above it stops.
        home = joinpath(space, "home")
        repo = joinpath(home, "repo")
        registry = joinpath(home, ".fair", "registry")
        pristine = joinpath(space, "db.sqlite3.pristine")
        bin = juliawrapper(space)
        # Every `fair` command, and what it runs, has the test's home, this `julia` first on the
        # path, its depot, a load path without the test sandbox, so that a run uses its own
        # environment, the small paper figures, all of them, and no display.
        environment = ["HOME" => home, "USERPROFILE" => home,
            "PATH" => string(bin, ":", Sys.BINDIR, ":", get(ENV, "PATH", "")),
            "JULIA_DEPOT_PATH" => join(DEPOT_PATH, ":"),
            "JULIA_LOAD_PATH" => "@:@stdlib", "JULIA_PROJECT" => nothing,
            "ECOSISTEM_SCALE" => "small", "ECOSISTEM_PAPER_FIGURES" => nothing,
            "GKSwstype" => "100",
            "GIT_AUTHOR_NAME" => "Interface Test",
            "GIT_AUTHOR_EMAIL" => "test@noreply.com",
            "GIT_COMMITTER_NAME" => "Interface Test",
            "GIT_COMMITTER_EMAIL" => "test@noreply.com"]
        withenv(environment...) do
            fair"--version"
            if !isfile(joinpath(registry, "manage.py")) || !isfile(pristine)
                @info "Installing the test's own registry at $registry, once: a few minutes"
                rm(registry, recursive = true, force = true)
                mkpath(home)
                cd(space) do
                    return fair"registry install --version $REGISTRY_VERSION"
                end
                cp(joinpath(registry, "db.sqlite3"), pristine, force = true)
            end
            # An empty registry and store, so that every run of this file takes the same path. A
            # registry an interrupted run left running is stopped first, by its install's own script,
            # which stops nothing else: `fair registry stop` before `fair init` would ask whatever
            # server holds port 8000, a registry of the user's own perhaps.
            run(`$(joinpath(registry, "scripts", "stop_fair_registry"))`)
            cp(pristine, joinpath(registry, "db.sqlite3"), force = true)
            rm(joinpath(home, "data_store"), recursive = true, force = true)

            copyrepository(repo)
            # The examples' environments take this package from the checkout, so that what the
            # runs execute is what this process covers. That changes tracked files of the copy,
            # hence `--dirty` on every run.
            julia = joinpath(bin, "julia")
            for example in ("pipeline", "paper")
                project = joinpath(repo, "examples", example)
                run(`$julia --project=$project -e "using Pkg; Pkg.develop(path = $(repr(pkgdir(EcoSISTEM)))); Pkg.instantiate()"`)
            end

            # Searched for from 8030, clear of DataPipeline's own test harness, which starts at 8010.
            port, server = listenany(ip"127.0.0.1", 8030)
            close(server)
            config = cliconfig(joinpath(space, "cli-config.yaml"); home, repo,
                               registry,
                               port)
            cd(repo) do
                fair"init --using $config --local"
                # Started once and held, so the registry can be asked between the runs; starting
                # writes its token.
                fair"registry start"
                try
                    token = strip(read(joinpath(registry, "token"), String))
                    endpoint = D.RegistryEndpoint("http://127.0.0.1:$port/api/",
                                                  token = token)
                    @testset "Fetch runs" begin
                        fair"run --local --dirty examples/pipeline/fetch-worldclim.yaml"
                        fair"run --local --dirty examples/pipeline/fetch-naturalearth.yaml"
                        fair"run --local --dirty examples/pipeline/fetch-box.yaml"
                        runs = coderuns(endpoint, repo)
                        @test length(runs) == 3
                        # Each source's files and their records, in that source's namespace.
                        @test Set(named.(runs[1].outputs)) ==
                              Set(["WorldClim:worldclim/data/wc2.1_10m_bio_1",
                                      "WorldClim:worldclim/data/wc2.1_10m_bio_12",
                                      "WorldClim:worldclim/provenance/wc2.1_10m_bio_1",
                                      "WorldClim:worldclim/provenance/wc2.1_10m_bio_12"
                                  ])
                        @test Set(named.(runs[2].outputs)) ==
                              Set([
                                      "NaturalEarth:naturalearth/data/ne_10m_admin_0_map_units",
                                      "NaturalEarth:naturalearth/provenance/ne_10m_admin_0_map_units"
                                  ])
                        # The file of your own, with the source it was cut from.
                        @test Set(named.(runs[3].outputs)) ==
                              Set(["NaturalEarth:box/data/scotland_box",
                                      "NaturalEarth:box/provenance/scotland_box"
                                  ])
                        # With everything registered, a fetch run stops in its rehearsal.
                        fair"run --local --dirty examples/pipeline/fetch-worldclim.yaml"
                        fair"run --local --dirty examples/pipeline/fetch-box.yaml"
                        @test length(coderuns(endpoint, repo)) == 3
                    end

                    @testset "A download of your own, rehearsed" begin
                        # A table fetched from a URL through a run's track: a run writing it and its
                        # record fetches and registers it, and a run doing neither is refused.
                        # Rehearsals register nothing, so the runs on either side are untouched.
                        served = joinpath(home, "served")
                        mkpath(served)
                        cp(joinpath(pkgdir(EcoSISTEM), "data", "NaturalEarth",
                                    "regions.csv"),
                           joinpath(served, "regions.csv"), force = true)
                        downloads = joinpath(home, "downloads")
                        rm(downloads, recursive = true, force = true)
                        facts = EcoSISTEM.SourceFacts(EcoSISTEM.datasetinfo(EcoSISTEM.NaturalEarthLevel),
                                                      title = "Natural Earth's regions (regions)")
                        table = EcoSISTEM.Provenancer.URLFile("file://" *
                                                              joinpath(served,
                                                                       "regions.csv"),
                                                              joinpath(downloads,
                                                                       "regions.csv"),
                                                              facts = facts)
                        metadata = """
                        run_metadata:
                          description: A table of your own downloaded from a URL
                          default_input_namespace: NaturalEarth
                          default_output_namespace: NaturalEarth
                          script: |
                            true
                        """
                        writes = joinpath(repo, "fetch-table.yaml")
                        write(writes,
                              metadata * """
                              write:
                              - data_product: table/data/*
                                description: A table of your own
                                file_type: csv
                              - data_product: table/provenance/*
                                description: EcoSISTEM's provenance record of a file under table
                                file_type: toml
                              """)
                        neither = joinpath(repo, "neither-table.yaml")
                        write(neither, metadata)
                        D.rehearse(writes, local_only = true) do handle
                            track = EcoSISTEM.Track(handle, "table")
                            @test only(EcoSISTEM.missingfiles(table,
                                                              track = track)).reason ==
                                  "not held: the fetch run fetches and registers it"
                            r = EcoSISTEM.Provenancer.fetchfile(table,
                                                                role = :region,
                                                                track = track)
                            @test r.location == table.path &&
                                  isfile(table.path) &&
                                  !isnothing(r.fetched) && r.primary === true
                            return nothing
                        end
                        D.rehearse(neither, local_only = true) do handle
                            @test_throws "neither writes nor reads" EcoSISTEM.Provenancer.fetchfile(table,
                                                                                                    role = :region,
                                                                                                    track = EcoSISTEM.Track(handle,
                                                                                                                            "table"))
                            return nothing
                        end
                    end

                    @testset "Model run" begin
                        fair"pull --local examples/pipeline/model.yaml"
                        fair"run --local --dirty examples/pipeline/model.yaml"
                        runs = coderuns(endpoint, repo)
                        @test length(runs) == 4
                        @test named.(runs[4].outputs) ==
                              ["EcoSISTEM:Scotland-plot"]
                        trace = D.provenance(endpoint,
                                             "EcoSISTEM:Scotland-plot")
                        figure = trace.products[only(trace.roots)]
                        model = trace.runs[figure.made_by]
                        inputs = [trace.products[url] for url in model.inputs]
                        @test Set(named.(inputs)) ==
                              Set(vcat(named.(runs[1].outputs),
                                       named.(runs[2].outputs),
                                       named.(runs[3].outputs)))
                        # Each derived from its own source's deposit alone, which has no file.
                        @testset "$(named(input))" for input in inputs
                            fetch = trace.runs[input.made_by]
                            deposit = trace.products[only(fetch.inputs)]
                            root = first(split(input.name, '/'))
                            @test named(deposit) ==
                                  "$(input.namespace):$root/deposit"
                            @test isnothing(deposit.hash)
                            @test !isnothing(deposit.source)
                        end
                    end

                    @testset "Paper figures" begin
                        fair"pull --local examples/paper/paper.yaml"
                        fair"run --local --dirty examples/paper/paper.yaml"
                        runs = coderuns(endpoint, repo)
                        @test length(runs) == 5
                        outputs = ["EcoSISTEM:Examples/Figures/" .*
                                   ["NicheWidth", "IslandTesting",
                                       "DispersalSD",
                                       "OptVarPanel", "InvasionCircle",
                                       "Invasion"];
                                   "EcoSISTEM:Examples/Tables/" .*
                                   ["IslandParameters", "ContinentParameters"];
                                   "EcoSISTEM:Examples/Config"]
                        @test Set(named.(runs[5].outputs)) == Set(outputs)
                        # Each derived from the parameters alone.
                        @testset "$output" for output in outputs
                            trace = D.provenance(endpoint, output)
                            run = trace.runs[trace.products[only(trace.roots)].made_by]
                            @test [named(trace.products[url])
                                   for url in run.inputs] ==
                                  ["EcoSISTEM:Examples/Parameters"]
                        end
                    end
                    @testset "A file given by hand, and a deposit cited apart" begin
                        # Written into the copy: a script recording a table given by hand, and one
                        # made for a paper, each in a fetch run of its own, and reading the first
                        # back in a model run; the configs; and the citations their facts name by
                        # key, so nothing is looked up.
                        given = joinpath(repo, "given")
                        mkpath(given)
                        write(joinpath(given, "citations.bib"),
                              """
                              @misc{giventable,
                                author = {Test, Interface},
                                title = {A table given by hand},
                                publisher = {EcoSISTEM tests},
                                date = {2026-01-01},
                                url = {https://example.org/given-table}
                              }
                              @misc{citedtable,
                                author = {Test, Interface},
                                title = {A table made for a paper},
                                publisher = {EcoSISTEM tests},
                                date = {2026-01-01},
                                url = {https://example.org/cited-table}
                              }
                              @misc{citedpaper,
                                author = {Test, Interface},
                                title = {The paper the table was made for},
                                publisher = {EcoSISTEM tests},
                                date = {2026-01-01},
                                url = {https://example.org/cited-paper}
                              }
                              """)
                        write(joinpath(given, "given.jl"),
                              raw"""
                              using DataPipeline
                              using EcoSISTEM
                              using EcoSISTEM.Provenancer: SourceFacts, DepositIsFile,
                                                          DepositCitedApart, record_provenance

                              handle = DataPipeline.initialise()
                              DataPipeline.register_citations!(handle,
                                                               joinpath(@__DIR__, "citations.bib"))
                              kind = only(ARGS)
                              given = SourceFacts(source = "A table given by hand", key = "giventable",
                                                  title = "A table given by hand", identifies = :data,
                                                  releasedate = "2026-01-01", namespace = "Given",
                                                  deposit = DepositIsFile())
                              if kind == "cited"
                                  path = joinpath(@__DIR__, "made.csv")
                                  write(path, "a,b\n3,4\n")
                                  paper = SourceFacts(source = "The paper", key = "citedpaper",
                                                      identifies = :paper, releasedate = "2026-01-01")
                                  facts = SourceFacts(source = "A table made for a paper",
                                                      key = "citedtable",
                                                      title = "A table made for a paper",
                                                      identifies = :data, releasedate = "2026-01-01",
                                                      namespace = "Given",
                                                      deposit = DepositCitedApart(paper))
                                  record_provenance(path, facts, role = :species,
                                                    track = EcoSISTEM.Track(handle, "cited"))
                              else
                                  path = joinpath(@__DIR__, "table.csv")
                                  isfile(path) || write(path, "a,b\n1,2\n")
                                  r = record_provenance(path, given, role = :species,
                                                        track = EcoSISTEM.Track(handle, "given"))
                                  kind == "model" &&
                                      write(DataPipeline.link_write!(handle, "given-read"),
                                            string(countlines(r.location)))
                              end
                              DataPipeline.finalise(handle)
                              """)
                        metadata(description, argument) = """
                        run_metadata:
                          description: $description
                          default_input_namespace: Given
                          default_output_namespace: Given
                          script: |
                            julia --project=examples/pipeline given/given.jl $argument
                        register:
                        - namespace: Given
                          full_name: Given by hand
                          website: https://example.org/
                        """
                        record(root) = """
                        - data_product: $root/provenance/*
                          description: EcoSISTEM's provenance record of a file under $root
                          file_type: toml
                        """
                        write(joinpath(given, "fetch-given.yaml"),
                              metadata("A table given by hand, registered as its source itself",
                                       "given") * "write:\n" * record("given"))
                        write(joinpath(given, "fetch-cited.yaml"),
                              metadata("A table made for a paper, its deposit the paper",
                                       "cited") *
                              """
                              write:
                              - data_product: cited/data/*
                                description: A table made for a paper
                                file_type: csv
                              """ * record("cited"))
                        write(joinpath(given, "model-given.yaml"),
                              metadata("The table given by hand, read back",
                                       "model") *
                              """
                              read:
                              - data_product: given/*/*
                                use:
                                  namespace: Given
                              write:
                              - data_product: given-read
                                description: The number of lines in the table given by hand
                                file_type: txt
                              """)
                        fair"run --local --dirty given/fetch-given.yaml"
                        fair"run --local --dirty given/fetch-cited.yaml"
                        fair"pull --local given/model-given.yaml"
                        fair"run --local --dirty given/model-given.yaml"
                        runs = coderuns(endpoint, repo)
                        @test length(runs) == 8
                        # The file given by hand is its run's one input, registered as its source
                        # and primary, with no deposit; its record is the one output, derived from
                        # it.
                        @test named.(runs[6].outputs) ==
                              ["Given:given/provenance/table"]
                        trace = D.provenance(endpoint,
                                             "Given:given/provenance/table")
                        fetch = trace.runs[trace.products[only(trace.roots)].made_by]
                        table = trace.products[only(fetch.inputs)]
                        @test named(table) == "Given:given/data/table" &&
                              isnothing(table.made_by) && !isnothing(table.hash)
                        @test occursin("given-table",
                                       table.source.identifier) &&
                              table.source.primary
                        # A table made for a paper is derived from a deposit with no file, the
                        # paper, and keeps its own source.
                        @test Set(named.(runs[7].outputs)) ==
                              Set([
                                      "Given:cited/data/made",
                                      "Given:cited/provenance/made"
                                  ])
                        trace = D.provenance(endpoint, "Given:cited/data/made")
                        made = trace.products[only(trace.roots)]
                        deposit = trace.products[only(trace.runs[made.made_by].inputs)]
                        @test named(deposit) == "Given:cited/deposit" &&
                              isnothing(deposit.hash) &&
                              occursin("cited-paper",
                                       deposit.source.identifier) &&
                              !deposit.source.primary
                        @test occursin("cited-table", made.source.identifier) &&
                              made.source.primary
                        # A model run reads the file registered as its source, and its record.
                        @test named.(runs[8].outputs) == ["Given:given-read"]
                        trace = D.provenance(endpoint, "Given:given-read")
                        read = trace.runs[trace.products[only(trace.roots)].made_by]
                        @test Set(named(trace.products[url])
                                  for url in read.inputs) ==
                              Set(["Given:given/data/table",
                                      "Given:given/provenance/table"])
                    end
                finally
                    fair"registry stop"
                end
            end
        end
    end
end

end
