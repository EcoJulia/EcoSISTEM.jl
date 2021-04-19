# Datasets currently in use by EcoSISTEM.jl

Preliminary list of parameters/datasets.

| Name of parameter / dataset | Description | Value | Source | Other info (e.g. stability) |
|----------------------------|-------------|-------|--------|-----------------------------|
|  p_s                          |     Probability of developing symptoms        |      0.96 |  [From Thibaud's original model](http://gabgoh.github.io/COVID/index.html)      |                              |
|          p\_h                  |       Probability of hospitalisation      |   0.2    |    Guess    |                             |
|           cfr\_home                 |       Case fatality ratio (at home)      |      0.1 |     Guess   |                             |
|             cfr\_hospital               |      Case fatality ratio (at hospital)       |   0.1    |    Guess    |                             |
|              T\_lat              |      Latent period       |    5 days   |      [From Thibaud's original model](http://gabgoh.github.io/COVID/index.html)  |                             |
|               T\_asym             |      Asymptomatic period       |   3 days    |    [From Thibaud's original model](http://gabgoh.github.io/COVID/index.html)   |                             |
|              T_sym              |       Symptomatic period      |    5 days   |     [From Thibaud's original model](http://gabgoh.github.io/COVID/index.html)   |                             |
|              T\_hosp              |      Hospitalisation period       |   5 days    |   [From Thibaud's original model](https://www.icnarc.org/Our-Audit/Audits/Cmp/Reports) |                             |
|              T\_rec              |     Recovery period        |   11 days    |     [From Thibaud's original model](http://gabgoh.github.io/COVID/index.html)    |                             |
|              mu\_1              |      Probability of becoming Asymptomatic       |    1/T\_lat   |        |                             |
|              mu\_2              |      Probability of becoming Symptomatic       |    p\_s \* 1/T\_asym   |        |                             |
|              hospitalisation              |       Probability of becoming Hospitalised      |   p\_h \* 1/T\_sym    |        |                             |
|           sigma\_1                 |       Probability of Recovery from Asymptomatic      |    (1 - p\_s) * 1/T\_asym   |        |                             |
|           sigma\_2                 |     Probability of Recovery from Symptomatic        |    (1 - p\_h) \* (1 - cfr\_home) * 1/T\_rec   |        |                             |
|           sigma_hospital                 |    Probability of Recovery from Hospital         |   (1 - cfr_hosp) \* 1/T\_hosp    |        |                             |
|            death_home                |    Probability of Death at home         |    cfr\_home \* 2/T\_hosp   |        |                             |
|          death_hospital                  |    Probability of Death at hospital         |   cfr\_hosp \* 1/T\_hosp    |        |                             |
|             ScotlandDensity2011               |     Scottish population density at 1km grid        |       |    UK census 2011 - A Reeves 'Covid19-ScottishCensusData' repo    |                            |
|            dispersal\_dist                |       Average dispersal distance of virus per disease category      |    2.0km per infectious disease category   |       Guess |        Varies depending on grid size                     |
|              mean\_pref              |      Mean temperature preference of virus       |   298K    |   Guess     |         Currently tuned to fit environment perfectly                    |
|            var\_pref                |     Temperature niche width of virus        |   0.1K    |   Guess     |      Currently tuned to fit environment perfectly                       |
|            birth                |     Probability of giving birth per individual       |   1.3e-4/day (20-40 year olds), 0 otherwise    |     Guess   |                             |
|             death               |     Probability of giving natural mortality per individual        |    2.7e-5/day   |  Guess     |                             |
|           virus\_growth\_asymp                 |     Rate of generating virus per asymptomatic individual        |   0.1/day    |       Guess |                             |
|           virus\_growth\_symp                 |    Rate of generating virus per symptomatic individual         |   0.1/day    |    Guess    |                             |
|           beta\_force                 |    Force of infection         |  10.0/day     |   Guess     |                             |
|           beta\_env                 |    Environmental transmission         |  10.0/day     |   Guess     |                             |
