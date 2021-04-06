**Table of datasets currently in use by EcoSISTEM.jl**

Preliminary list of parameters/datasets.

| Name of parameter/ dataset | Description | Value | Source | Other info (e.g. stability) |
|----------------------------|-------------|-------|--------|-----------------------------|
|  p_s                          |     Probability of developing symptoms        |      0.96 |  http://gabgoh.github.io/COVID/index.html (From Thibaud's original model)      |                              |
|          p_h                  |       Probability of hospitalisation      |   0.2    |    Guess    |                             |
|           cfr_home                 |       Case fatality ratio (at home)      |      0.1 |     Guess   |                             |
|             cfr_hospital               |      Case fatality ratio (at hospital)       |   0.1    |    Guess    |                             |
|              T_lat              |      Latent period       |    5 days   |      http://gabgoh.github.io/COVID/index.html (From Thibaud's original model)  |                             |
|               T_asym             |      Asymptomatic period       |   3 days    |    http://gabgoh.github.io/COVID/index.html (From Thibaud's original model)   |                             |
|              T_sym              |       Symptomatic period      |    5 days   |     http://gabgoh.github.io/COVID/index.html (From Thibaud's original model)   |                             |
|              T_hosp              |      Hospitalisation period       |   5 days    |   https://www.icnarc.org/Our-Audit/Audits/Cmp/Reports (From Thibaud's original model) |                             |
|              T_rec              |     Recovery period        |   11 days    |    http://gabgoh.github.io/COVID/index.html (From Thibaud's original model)    |                             |
|              mu_1              |      Probability of becoming Asymptomatic       |    1/T_lat   |        |                             |
|              mu_2              |      Probability of becoming Symptomatic       |    p_s * 1/T_asym   |        |                             |
|              hospitalisation              |       Probability of becoming Hospitalised      |   p_h * 1/T_sym    |        |                             |
|           sigma_1                 |       Probability of Recovery from Asymptomatic      |    (1 - p_s) * 1/T_asym   |        |                             |
|           sigma_2                 |     Probability of Recovery from Symptomatic        |    (1 - p_h) * (1 - cfr_home) * 1/T_rec   |        |                             |
|           sigma_hospital                 |    Probability of Recovery from Hospital         |   (1 - cfr_hosp) * 1/T_hosp    |        |                             |
|            death_home                |    Probability of Death at home         |    cfr_home * 2/T_hosp   |        |                             |
|          death_hospital                  |    Probability of Death at hospital         |   cfr_hosp * 1/T_hosp    |        |                             |
|             ScotlandDensity2011               |     Scottish population density at 1km grid        |       |    UK census 2011 - A Reeves 'Covid19-ScottishCensusData' repo    |                            |
|            dispersal_dist                |       Average dispersal distance of virus per disease category      |    2.0km per infectious disease category   |       Guess |        Varies depending on grid size                     |
|              mean_pref              |      Mean temperature preference of virus       |   298K    |   Guess     |         Currently tuned to fit environment perfectly                    |
|            var_pref                |     Temperature niche width of virus        |   0.1K    |   Guess     |      Currently tuned to fit environment perfectly                       |
|            birth                |     Probability of giving birth per individual       |   1.3e-4/day (20-40 year olds), 0 otherwise    |     Guess   |                             |
|             death               |     Probability of giving natural mortality per individual        |    2.7e-5/day   |  Guess     |                             |
|           virus_growth_asymp                 |     Rate of generating virus per asymptomatic individual        |   0.1/day    |       Guess |                             |
|           virus_growth_symp                 |    Rate of generating virus per symptomatic individual         |   0.1/day    |    Guess    |                             |
|           beta_force                 |    Force of infection         |  10.0/day     |   Guess     |                             |
|           beta_env                 |    Environmental transmission         |  10.0/day     |   Guess     |                             |
