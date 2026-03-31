# BayesGreySealPupModel
Bayesian state-space model to estimate grey seal pup production from serial aerial surveys, with survival and carcass persistence explicitly modelled. Code and data accompanying Goldman et al. (submitted to RSOS). 

R script and model code descriptions: 

## Stan models and fit code
* `greysealpups_3stage.stan` - Mirrors the Russell et al. (2019) and Jacobson et al. (2025) maximum likelihood models. Dead pups are added to whitecoat counts (Objective 1)
* `greysealpups_Deads.stan` - Extends the standard model with explicit mortality (age-dependent survival on the logit scale) and carcass persistence. Dead pups are modelled as a separate observed state. TTM and TTL can be either fixed or estimated depending on the objective (2-5)
* `Stan/greysealpups_Deads_estTTL.stan` — Full model with estimated TTL, fixed TTM (Objective 5)
* `Stan/greysealpups_Deads_estTTMTTL.stan` — Full model with estimated TTM and TTL (Objective 5)
 
## Simulation code
* `HgSim_Intervals.R` – Contains the generate_obs() function used to simulate grey seal pup survey data
* `run_mortality_simulations.R` - runs 100 simulations for Objectives 1 and 2, fitting both the standard and full models

## Empirical data fit code 
* `fit_empirical.R` - fits the full model to Brownsman and Staple Island data (2023) using UAV and fixed-wing surveys with updated detection rates, and compares fixed vs estimated TTL and TTM (Objectives 4 and 5)

## Empirical data - available upon request

## Results
* `summaryObj1.csv` - Summarized posterior estimates from 100 simulations (mortality is simulated but omitted from the model), including mean, HDI, relative bias, and coverage for key parameters
* `summaryObj2.csv` – Summarized posterior estimates from 100 simulations (mortality is simulated and included in the model)


