using Simulation
using Diversity
using Statistics
using Compat.Test

eco = TestCache()
times = 1month:1month:10years
cache = CachedEcosystem(eco, "/Users/claireh/Documents/PhD/GIT/Simulation/test/Cache", times)
@test_nowarn cache = CachedEcosystem(eco, "test/Cache", times)

divset = DiversitySet(cache, collect(times))
@test_nowarn divset = DiversitySet(cache, collect(times))

function cachefun(cache::CachedEcosystem)
    for tm in gettimes(divset)
        updatesimulation!(cache, tm)
        newdiv = norm_sub_alpha(cache, 1.0)
        if ismissing(divset.data)
            divset.data = newdiv
        else
            append!(divset, newdiv)
        end
    end
end

@test_nowarn cachefun(cache)
