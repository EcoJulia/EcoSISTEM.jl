# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The data pipeline for the paper's HPC runs, set up once by `setup.bash`: DataPipeline's own copy
# of `fair`, a local registry, and this clone initialised with a CLI configuration written here, all
# under `ECOSISTEM_HPC_ROOT`. Run from the root of the clone, in `examples/paper`'s environment,
# after `env.sh`. `fair` and a new registry are each the newest release. Run again, it brings the
# copy of `fair` up to date and keeps the registry already installed, with its records.

using DataPipeline

const ROOT = ENV["ECOSISTEM_HPC_ROOT"]
const PORT = ENV["ECOSISTEM_HPC_PORT"]
const REGISTRY = ENV["FAIR_REGISTRY_DIR"]
const USERFILE = joinpath(ROOT, "fair", "user.yaml")

isfile(joinpath("examples", "paper", "paper.yaml")) ||
    error("run this from the root of the EcoSISTEM clone, not ", pwd())
isfile(USERFILE) ||
    error("no $USERFILE: copy examples/paper/hpc/fair-user.yaml there and fill it in")
any(line -> occursin(r"^[^#]*<", line), eachline(USERFILE)) &&
    error("$USERFILE still has a `<...>` to replace")

DataPipeline.install_cli()

# Without `--version`, `fair` installs the registry's newest release.
if !isfile(joinpath(REGISTRY, "manage.py"))
    mkpath(dirname(REGISTRY))
    fair"registry install --directory $REGISTRY"
end

# The CLI configuration: the local registry on its port and its store beside it, the public registry
# as the remote, this clone, and the user from `USERFILE`.
remote = try
    readchomp(`git remote get-url origin`)
catch
    "https://github.com/EcoJulia/EcoSISTEM.jl"
end
user = read(USERFILE, String)
config = joinpath(ROOT, "fair", "cli-config.yaml")
write(config, """
      namespaces:
        input: EcoSISTEM
        output: EcoSISTEM
      registries:
        local:
          uri: http://127.0.0.1:$PORT/api/
          directory: $REGISTRY
          data_store: $(joinpath(ROOT, "fair", "data_store"))/
          token: $(joinpath(REGISTRY, "token"))
        origin:
          uri: https://data.fairdatapipeline.org/api/
          data_store: https://data.fairdatapipeline.org/data/
          token: $(joinpath(ROOT, "fair", "remote-token"))
      git:
        local_repo: $(pwd())
        remote: origin
        remote_repo: $remote
      $user""")

# Asks nothing, given the configuration; for a clone already initialised, registers its user again.
fair"init --using $config --local"
println("The pipeline is set up: registry ",
        readchomp(`git -C $REGISTRY describe --tags`), " at ", REGISTRY,
        ", data store at ", joinpath(ROOT, "fair", "data_store"))
