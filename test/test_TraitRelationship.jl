using EcoSISTEM
using Test
using Distributions
using Unitful
using Unitful.DefaultSymbols

@test_nowarn Gauss{Unitful.Temperature}()
@test_nowarn Match{Int64}()
@test_nowarn Trapeze{Int64}()
@test_nowarn Unif{Int64}()
@test_nowarn NoRelContinuous{Int64}()
@test_nowarn NoRelDiscrete{Int64}()
