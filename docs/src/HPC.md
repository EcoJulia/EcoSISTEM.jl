# Running on DiRAC

In order to run EcoSISTEM.jl on DiRAC, it is first necessary to [request an account](https://safe.epcc.ed.ac.uk/dirac/)
and then to request access to the `dc003` project. Once that is granted, it is possible to log
into the server via ssh at the address
```
login.hpc.cam.ac.uk
```

Julia is already installed on the cluster, so the only things that are needed are a Project.toml
file and a Julia script defining and setting up the environment for the experiment, a Julia
script with the experiment, and a submission script. Here we give examples using the
[Scottish experiment](https://github.com/ScottishCovidResponse/EcoSISTEM.jl/blob/dev/examples/Epidemiology/Scotland_run.jl).
The other files needed for this experiment can be found [here](https://github.com/ScottishCovidResponse/EcoSISTEM.jl/blob/dev/examples/Epidemiology/HPC/)

## Project.toml

This file should be enough to setup the environment for the experiment. Thus, it should
change depending on the packages required by the experiment script. E.g.:

```julia
name = "RAMP"

[deps]
AxisArrays = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
julia = "1.1"

[extras]

[targets]
```

**This file should not include EcoSISTEM.jl, as it is unregistered.** This will be taken care
of by the setup script.

## Setup script

This file is responsible for cloning EcoSISTEM.jl, as it is unregistered, and also for
instantiating the environment.

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/ScottishCovidResponse/EcoSISTEM.jl.git"))
Pkg.instantiate()
```

## Experiment script

This file should contain the experiment to be run. There are [examples](https://github.com/ScottishCovidResponse/EcoSISTEM.jl/blob/dev/examples/Epidemiology/) available.

## Submission script

This file parses options, sets up the local directories and submits the job to the queue.

```sh
#!/bin/bash

# Usage: ./submit_scot --input <scriptfilename> --nprocs <numberofprocesses> --walltime <walltime> --dir <pathtorundir>

set -e
set -u

while [ $# -gt 0 ]; do
    if [ $1 = "--input" ]; then
        inputfile=$2
        shift 2
    elif [ $1 == "--samples" ]; then
        nsamples=$2
        shift 2
    elif [ $1 == "--dir" ]; then
        dir=$2
        shift 2
    elif [ $1 == "--nprocs" ]; then
        nprocs=$2
        shift 2
    elif [ $1 == "--walltime" ]; then
        walltime=$2
        shift 2
    else
        echo "Unrecognised arguments: $*" >&2
        exit 1
    fi
done

jobname=$(basename $dir)

mkdir $dir
mkdir $dir/Output
cp $inputfile $dir
cp Project.toml $dir
cp setup.jl $dir
inputbasename=$(basename $inputfile)

cat >$dir/submit.sh <<EOF
#!/bin/bash
#SBATCH --job-name $jobname
#SBATCH --account DIRAC-DC003-CPU
#SBATCH --ntasks $nprocs
#SBATCH --time $walltime
#SBATCH --mail-type ALL
#SBATCH --no-requeue
#SBATCH --partition skylake
#SBATCH --output log.txt
. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4
module load julia/1.4
export JULIA_NUM_THREADS=$nprocs
stdbuf -oL -eL julia --project=@. setup.jl
stdbuf -oL -eL julia --project=@. $inputbasename
EOF

cd $dir
chmod u+x submit.sh
sbatch submit.sh
```

Usage example:
```sh
./submit_scot --input Scotland_run.jl --nprocs 6 --walltime 1:00:00 --dir /home/username/run1
```
The command above must be run from the directory where `submit_scot`, `Scotland_run.jl` (the
experiment script in this case), `setup.jl` and `Project.toml` are. These will be copied to
a new folder at the provided path, `/home/username/run1`. At that path, a `log.txt` file
will contain the logs of the job. The `mail-type ALL` option makes it such that an email is
sent informing of the conclusion of the job, and it can be changed or removed.

## Further resources

More detailed information on the cluster, including further options and commands, can be found [here](https://docs.hpc.cam.ac.uk/hpc/user-guide/quickstart.html).
