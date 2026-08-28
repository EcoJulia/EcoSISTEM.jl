# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Run every Pluto notebook in `notebooks/` as an integration test, in the `notebooks/` environment.
#
# On its own, either through the test harness (which provisions the test environment):
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_notebooks.jl"])'
#
# or as a bare script, which works here because this file activates its own environment:
#
#     cd test && julia --project=.. extras_notebooks.jl
#
# `runtests.jl` also picks it up automatically, along with every other `test/extras_*.jl` — but only
# after the unit tests have passed, since a failing testset throws before the loop is reached.
#
# Skips itself on Windows — the notebooks download raster data, as `test_datasetread.jl` does. The
# skip says so rather than passing silently.

module ExtrasNotebooks

# Headless plots — see the fuller note in `runtests.jl`. Repeated here because this file is
# documented as a standalone script, and running it *is* running tests. `get!` is idempotent, so
# setting it in both places is harmless.
get!(ENV, "GKSwstype", "100")

using Test
using EcoSISTEM
using Pkg

if Sys.iswindows()
    @info "Skipping the notebooks: they download raster data, which Windows runners cannot do."
else
    # Restore whatever environment we were called with — see `extras_examples.jl` for why this
    # matters (the hygiene checks' packages live in the *test* environment, not this one).
    original_project = dirname(Base.active_project())
    try
        @testset "Notebooks folder" begin
            println()
            @info "Running from notebooks folder ..."
            notebookdir = pkgdir(EcoSISTEM, "notebooks")
            notebooks = filter(str -> occursin(r"\.jl$", str),
                               readdir(notebookdir))
            # Snapshot the notebook files. Running a Pluto notebook can rewrite it in place with an
            # embedded Project/Manifest, but the formatter / RSMD hygiene tests require the tracked
            # files to stay byte-identical — so restore them from this snapshot no matter what
            # happens below.
            originals = Dict(nb => read(joinpath(notebookdir, nb), String)
                             for nb in notebooks)
            try
                cd(notebookdir) do
                    Pkg.activate(notebookdir)
                    Pkg.resolve()
                    for nb in notebooks
                        println("    * Running $nb ...")
                        # These are Pluto notebooks: interactive `@bind` widgets evaluate to
                        # `missing` outside Pluto. Substitute each slider's HTML default value so
                        # the notebook runs non-interactively, from a temp copy so the tracked file
                        # is never touched.
                        runnable = replace(originals[nb],
                                           r"@bind\s+(\w+)\s+html\"[^\"]*value='([^']*)'[^\"]*\"" =>
                                               s"\1 = \2")
                        tmp = joinpath(tempdir(), "ecosistem_nb_" * nb)
                        write(tmp, runnable)
                        try
                            # `@test`, not `@test_nowarn`: notebooks legitimately log (data
                            # downloads, cell-size info) — we only require they run without error.
                            @test (include(tmp); true)
                        finally
                            rm(tmp, force = true)
                        end
                    end
                end
            finally
                # Restore the tracked notebook files (dropping any Pluto-generated
                # Project/Manifest) so the hygiene tests still see a clean repo.
                for nb in notebooks
                    write(joinpath(notebookdir, nb), originals[nb])
                end
            end
        end
    finally
        Pkg.activate(original_project)
    end
end

end
