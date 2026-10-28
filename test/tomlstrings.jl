# SPDX-License-Identifier: LGPL-3.0-or-later

# Every string a TOML file holds, at any depth, for the tests checking that a written record names
# nothing of the machine that wrote it. They search the parsed values, not the file's text: the
# text holds a Windows path's backslashes escaped, so a path would never be found in it there.
# `include`d into each test module that needs it. Not named `test_*.jl`, so `runtests.jl` neither
# runs it directly nor expects a matching `src/` file.

tomlstrings(path::AbstractString) = _tomlstrings(EcoSISTEM.TOML.parsefile(path))

_tomlstrings(x::AbstractString) = [String(x)]

function _tomlstrings(x::AbstractDict)
    return reduce(vcat, map(_tomlstrings, collect(values(x))), init = String[])
end

function _tomlstrings(x::AbstractVector)
    return reduce(vcat, map(_tomlstrings, x), init = String[])
end

_tomlstrings(::Any) = String[]
