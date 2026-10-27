# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEMDataPipelineExt

import DataPipeline
using DataPipeline: DataRegistryHandle, ConfigDataProduct, ConfigFileException,
                    ReadWriteException, link_read!, link_register!, link_write!,
                    register_authors!, register_citations!

using EcoSISTEM
using EcoSISTEM: Track, NaturalEarthSpec, RasterSpec

using Dates: Dates
using TOML: TOML

include("tracked.jl")

end
