# =============================================================================
# run_mortality_simulations.R
#
# Objectives 1 & 2: Evaluate bias from omitting vs explicitly modelling mortality.
#
# Simulates 100 datasets with plausible pup mortality and fits each with:
#   - greysealpups_3stage.stan  (Obj 1: standard model, deads added to whitecoats)
#   - greysealpups_Deads.stan   (Obj 2: full model, mortality explicitly modelled)
#
# Two run modes:
#   "rerun"  - run all 100 simulations from scratch 
#   "reload" - load previously saved .RData files and summarise
# =============================================================================

library(rstan)
library(furrr)
library(future)
library(dplyr)

source("R/HgSim_Intervals.R")

rstan_options(auto_write = TRUE)

# =============================================================================
# SETTINGS
# =============================================================================

run_mode        <- "reload"   # "reload" or "rerun"
num_simulations <- 100
output_dir_obj1 <- "output/mortality_FitNoDeadsSimDeads"
output_dir_obj2 <- "output/mortality_FitDeads"

# =============================================================================
# RERUN: Simulate and fit both models
# =============================================================================

if (run_mode == "rerun") {
  
  message("Running full simulations (this may take hours)...")
  
  # Compile models
  mod_3s    <- stan_model("stan/greysealpups_3stage.stan")
  mod_deads <- stan_model("stan/greysealpups_Deads.stan")
  
  plan(multisession, workers = parallel::detectCores() - 1)
  
  run_simulation <- function(sim_num, model_type = "Obj1") {
    
    sim_data <- generate_obs(
      seed       = sim_num,
      birth.type = "skewnorm",
      obs.type   = "poisson",
      obsdays    = seq(30, 80, by = 12),
      nborn      = 500,
      mu.bday    = 34,
      sd.bday    = 12,
      skew.bday  = 4,
      mu.leave   = 31.5,
      sd.leave   = 7,
      mu.moult   = 23,
      sd.moult   = 5,
      pow        = 0.95,
      pom        = 0.95,
      pod        = 0.95,
      pcw        = 1,
      pcm        = 0.91,
      use_pcw    = TRUE,
      add_dead   = TRUE
    )
    
    # Obj 1: deads added to whitecoats (standard model)
    # Obj 2: deads as separate observed state (full model)
    if (model_type == "Obj1") {
      obs_vector <- cbind(
        sim_data$obs_vector[, 1] + sim_data$obs_vector[, 3],
        sim_data$obs_vector[, 2]
      )
      datalist <- list(
        nstates     = 3,
        maxDayBirth = 100,
        ndays       = 123,
        nobsstates  = 2,
        nsurveys    = length(sim_data$obsdays),
        obs_vector  = obs_vector,
        obsdays     = sim_data$obsdays
      )
      mod <- mod_3s
    } else {
      obs_vector <- cbind(
        sim_data$obs_vector[, 1],
        sim_data$obs_vector[, 2],
        sim_data$obs_vector[, 3]
      )
      datalist <- list(
        nstates     = 5,
        maxDayBirth = 100,
        ndays       = 123,
        nobsstates  = 3,
        nsurveys    = length(sim_data$obsdays),
        obs_vector  = obs_vector,
        obsdays     = sim_data$obsdays
      )
      mod <- mod_deads
    }
    
    fit <- sampling(mod, datalist, iter = 2000, chains = 4)
    return(fit)
  }
  
  # Run Obj 1: standard model (no mortality modelled)
  all_fits_Obj1 <- future_map(
    1:num_simulations,
    ~ run_simulation(.x, "Obj1"),
    .options = furrr_options(seed = TRUE)
  )
  
  # Run Obj 2: full model (mortality modelled)
  all_fits_Obj2 <- future_map(
    1:num_simulations,
    ~ run_simulation(.x, "Obj2"),
    .options = furrr_options(seed = TRUE)
  )
  
  # Save
  dir.create(output_dir_obj1, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_dir_obj2, showWarnings = FALSE, recursive = TRUE)
  fit_list <- all_fits_Obj1
  save(fit_list, file = file.path(output_dir_obj1, "all_fits_Obj1.RData"))
  fit_list <- all_fits_Obj2
  save(fit_list, file = file.path(output_dir_obj2, "all_fits_Obj2.RData"))
  
  message("Simulations complete and saved.")
}

# =============================================================================
# RELOAD: Load previously saved fits
# =============================================================================

if (run_mode == "reload") {
  
  message("Loading saved fits...")
  
  load(file.path(output_dir_obj1, "all_fits_Obj1.RData"))
  all_fits_Obj1 <- fit_list
  rm(fit_list)
  
  load(file.path(output_dir_obj2, "all_fits_Obj2.RData"))
  all_fits_Obj2 <- fit_list
  rm(fit_list)
}

# =============================================================================
# SUMMARISE: Compute bias, coverage, CRIs across 100 simulations
# =============================================================================

true_values <- c(nborn = 500, muBday = 34, sdBday = 12, alphaBday = 4)

summarize_all_params <- function(fit_list, true_values) {
  
  param_names <- names(true_values)
  
  do.call(rbind, lapply(param_names, function(p) {
    
    all_samples <- lapply(fit_list, function(f) {
      arr <- as.array(f)
      if (!p %in% dimnames(arr)[[3]]) stop(paste("Parameter", p, "not found"))
      as.vector(arr[, , p])
    })
    
    means    <- sapply(all_samples, mean)
    ci_lower <- sapply(all_samples, function(x) quantile(x, 0.025))
    ci_upper <- sapply(all_samples, function(x) quantile(x, 0.975))
    
    mean_estimate <- mean(means)
    sd_estimate   <- sd(means)
    
    data.frame(
      Parameter      = p,
      True_Value     = true_values[p],
      Mean_Estimate  = round(mean_estimate, 2),
      SD             = round(sd_estimate, 2),
      Lower_95_CRI   = round(mean(ci_lower), 2),
      Upper_95_CRI   = round(mean(ci_upper), 2),
      Relative_Pct_Bias = round((mean_estimate - true_values[p]) / true_values[p] * 100, 2),
      Coverage       = round(mean(ci_lower <= true_values[p] & ci_upper >= true_values[p]), 2)
    )
  }))
}

# Manuscript Table 3: Standard model (Obj 1, without mortality)
summary_Obj1 <- summarize_all_params(all_fits_Obj1, true_values)
print(summary_Obj1)

# Manuscript Table 4: Full model (Obj 2, with mortality)
summary_Obj2 <- summarize_all_params(all_fits_Obj2, true_values)
print(summary_Obj2)

# Save summary CSVs
write.csv(summary_Obj1, file.path(output_dir_obj1, "summary_Obj1.csv"), row.names = FALSE)
write.csv(summary_Obj2, file.path(output_dir_obj2, "summary_Obj2.csv"), row.names = FALSE)
