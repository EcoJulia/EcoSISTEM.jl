# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEMDataPipelineExt

# `import DataPipeline` is load-bearing and easy to lose: a trigger package is *loaded* when the
# extension activates, but its name is only in scope if the module asks for it. Without this line
# `assetdir(DataPipeline)` below is an `UndefVarError` — and the module still precompiles, because a
# function body resolves at call time, so the failure appears only when `unziptemp` is finally used.
import DataPipeline
import p7zip_jll

using EcoSISTEM

# Unzip `path` into a fresh directory under EcoSISTEM's DataPipeline scratch area, and return it.
#
# The destination is a `mktempdir` under `assetdir(DataPipeline)` rather than the system temporary
# directory, so an unpacked archive lands with the rest of the pipeline's cached data and is cleaned
# up on the same terms.
function EcoSISTEM.unziptemp(path::String)
    newpath = mktempdir(EcoSISTEM.assetdir(DataPipeline))
    run(`$(p7zip_jll.p7zip()) x -tzip -y -o$(newpath) $(path)`)
    @debug "Unzipped to $newpath"
    return newpath
end

end
