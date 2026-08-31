# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Runs the five landscape configurations in `examples/landscapes/`, recreated from
# the former `test/examples/{Island,Patch,Peaked,RegionEven,RegionGrad}.jl` (since removed).
#
# The point of the set is the **comparison**: the same ecology under different environments and
# different boundary conditions. So this does not merely check each builds - it checks they behave
# differently, which is the property that would silently rot if only construction were tested.
#
# **Scale is decided by `ECOSISTEM_SCALE`** (see `landscapes/configurations.jl`): the test suite sets
# it to `small`, and running this file directly gets the originals' full size - 10,000 species on a
# 50 × 50 grid for `peaked`. That is why everything below is derived from **one** simulation per
# landscape: at full scale a second run of `peaked` for the sake of an assertion is not cheap.

# **A module, deliberately.** `test/extras_examples.jl` includes every top-level example into
# **one** module, and more than one of them defines `runscale()` and `CONFIGURATIONS`. Julia 1.12
# lets a later `const` silently overwrite an earlier one - no warning - so without the wrapper these
# examples would quietly reconfigure each other depending on the order they happened to run in.

module Landscapes

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

include(joinpath(@__DIR__, "landscapes", "configurations.jl"))

const DURATION = 2.0year
const TIMESTEP = 1.0month_mean_duration

# `println`, not `@info`: logging goes to **stderr**, and `test/extras_examples.jl` wraps each
# example in `@test_nowarn`, which fails on *any* stderr output rather than only on warnings.
println("Running the landscapes at $(runscale()) scale.")

# --- one pass: build and simulate each landscape exactly once ----------------
results = map(LANDSCAPES) do name
    configuration = _configuration(name)
    eco = landscape(name)
    simulate!(eco, DURATION, TIMESTEP)
    return (name = name, eco = eco, grid = configuration.grid,
            total = sum(eco.abundances.matrix),
            occupied = count(>(0), sum(eco.abundances.matrix, dims = 1)))
end
results = NamedTuple{LANDSCAPES}(results)

for r in results
    @assert r.total > 0 "$(r.name) went extinct - nothing left to compare"
    @assert r.occupied > 0
end

# --- an island is not a patch ------------------------------------------------
# The whole point of having both. They are the same habitat, the same size and the same species;
# the island just loses its dispersers at the edge instead of reflecting them inland
# (`Island(false)`), and that alone costs it about a third of its abundance.
#
# This assertion is only meaningful *because* of the `false`. Under a plain `Island()` nothing
# is lost at the edge, the two totals land ~1% apart, and the sign flips with the seed - so asserting
# a difference there would have been asserting on RNG noise, and would have passed whether or not the
# boundary did anything. Measured across five seeds at small scale, the ratio sits at 0.68 every time.
let ratio = results.island.total / results.patch.total
    @assert 0.4 < ratio < 0.9 "island/patch abundance ratio was $ratio"
end

# --- the environments differ too, not just the boundaries --------------------
# Identical species and boundary; only the temperature field differs, so the outcome must too.
@assert results.region_even.total != results.region_gradient.total

# --- the boundary geometry, asserted exactly ---------------------------------
# Separately from how many individuals survive, the boundary decides *where* they may go, and that
# is exact rather than statistical: `calc_lookup_moves!` zeroes the probability of every offset that
# would leave the grid. From a corner cell that blocks the offsets pointing outward - all of them for
# `Island`, the north-south ones only for `Cylinder`, none at all for `Torus`.
#
# Asserted as an **ordering**, not as fixed fractions. The blocked share depends on how far the
# dispersal kernel reaches relative to the grid, and the full-scale configurations have much larger
# cells (`Peaked` is 9 km against `Island`'s 1 km), so a band that holds at one scale need not hold
# at the other. The ordering is a property of the boundary itself and holds at both.
let blocked = map(LANDSCAPES) do name
        eco = results[name].eco
        bound = eco.habitat.topology
        EcoSISTEM.calc_lookup_moves!(bound, 1, 1, 1, eco, 0)   # from the corner
        lookup = EcoSISTEM.getlookup(eco, 1)
        return count(iszero, lookup.pnew) / length(lookup.pnew)
    end
    blocked = NamedTuple{LANDSCAPES}(blocked)

    # A torus wraps both ways, so nothing is ever out of reach.
    @assert iszero(blocked.patch)
    @assert iszero(blocked.peaked)
    # A cylinder blocks only north-south, `Island` blocks every direction - so neither can block
    # less than the torus, and the cylinder cannot block more than the unbounded case.
    @assert blocked.island >= blocked.region_even >= blocked.patch
    @assert blocked.region_gradient == blocked.region_even
end

# --- a peak is not a gradient ------------------------------------------------
# Worth asserting rather than assuming: a gradient has one warm end, a peak has a warm centre with
# cool edges on *both* sides, so the same warm-adapted species is squeezed rather than pushed.
#
# Both specs vary along the **row** (`Y`) axis and are constant along the columns, and the regime
# matrix is indexed `(y, x)` as everywhere else in the package - so these read down a column, not
# across a row.
let peakcol = parent(results.peaked.eco.habitat.regime.matrix)[:, 1],
    gradcol = parent(results.region_gradient.eco.habitat.regime.matrix)[:, 1]

    mid = length(peakcol) ÷ 2 + 1
    @assert peakcol[mid] > peakcol[1]     # warm centre ...
    @assert peakcol[mid] > peakcol[end]   # ... cool at *both* ends

    # Whereas the gradient rises monotonically from one end to the other - one warm edge, not two
    # cool ones. This is the distinction the two spec types exist to draw.
    @assert issorted(gradcol)
    @assert gradcol[end] > gradcol[1]
end

end
