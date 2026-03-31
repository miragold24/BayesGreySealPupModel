# =============================================================================
# fit_empirical.R
#
# Fit models to Brownsman & Staple Island data (2023).
# Objective 4: Full model with fixed TTM/TTL, previous and updated rates.
# Objective 5: Full model with estimated TTL; estimated TTM+TTL (force moult).
#
# =============================================================================

library(rstan)
library(tidyverse)
library(HDInterval)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

# =============================================================================
# Stan models
# =============================================================================

mod_deads           <- stan_model("Stan/greysealpups_Deads.stan")
mod_deads_estTTL    <- stan_model("Stan/greysealpups_Deads_estTTL.stan")
mod_deads_estTTMTTL <- stan_model("Stan/greysealpups_Deads_estTTMTTL.stan")

# =============================================================================
# Detection rates
# =============================================================================

# Previous rates (Russell et al. 2019)
det_prev <- list(pobs_w = 0.95, pobs_m = 0.95, pobs_d = 0.95, pcm = 0.91, pcw = 1.0)

# Updated fixed-wing rates (Russell et al. 2025, SCOS 25/06)
det_PAS  <- list(pobs_w = 0.98, pobs_m = 0.90, pobs_d = 0.85, pcm = 0.80, pcw = 1.0)

# Updated UAV rates (Russell et al. 2025, SCOS 25/06; Appendix S2)
det_UAV  <- list(pobs_w = 0.98, pobs_m = 0.95, pobs_d = 0.93, pcm = 0.78, pcw = 1.0)

# =============================================================================
# Functions
# =============================================================================

prepare_data <- function(df, add_dead = TRUE) {
  df <- df %>%
    arrange(Year, Date) %>%
    group_by(Year) %>%
    mutate(
      start_date = as.Date(paste(Year, 9, 1, sep = "-")),
      nday = as.numeric(difftime(Date, start_date, units = "days")) + 1
    ) %>%
    ungroup()
  
  obs_matrix <- if (add_dead) {
    with(df, cbind(Whitecoats, Moulteds, DeadPups))
  } else {
    with(df, cbind(Whitecoats + DeadPups, Moulteds))
  }
  
  list(
    raw = df,
    datalist = list(
      nstates    = ifelse(add_dead, 5, 3),
      maxDayBirth = 100,
      ndays      = max(df$nday),
      nobsstates = ncol(obs_matrix),
      nsurveys   = nrow(df),
      obs_vector = obs_matrix,
      obsdays    = df$nday
    )
  )
}

add_det_rates <- function(datalist, det) {
  c(datalist, det)
}

fit_model <- function(model, datalist, det_rates, pars, out_path = NULL) {
  dl <- add_det_rates(datalist, det_rates)
  fit <- sampling(model, data = dl, iter = 2000, chains = 4,
                  control = list(adapt_delta = 0.95), pars = pars)
  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(fit, out_path)
  }
  fit
}

# =============================================================================
# Load data (Brownsman & Staple, 2023)
# =============================================================================

# UAV surveys 
uav_raw <- read.csv("Data/Count_NEEbs23_UAV-all_summary.csv")
uav_raw$Date <- as.Date(uav_raw$Date)
names(uav_raw)[names(uav_raw) == "Whites"]    <- "Whitecoats"
names(uav_raw)[names(uav_raw) == "DeadWhite"] <- "DeadPups"

# Fixed-wing surveys 
pas_raw <- read.csv("Data/HgPupCountsForProdModel_2023_Farnes-ALL-PAS_20250429-2djfr_H4DreplacePAS.csv")
pas_raw$Date <- as.Date(pas_raw$Date, format = "%d/%m/%Y")
pas_raw <- subset(pas_raw, Colony %in% c("Brownsman & Staple"))[-5, ]

uav_dat <- prepare_data(uav_raw, add_dead = TRUE)
pas_dat <- prepare_data(pas_raw, add_dead = TRUE)

# =============================================================================
# Objective 4: Fixed TTM/TTL, compare detection rates 
# =============================================================================

pars_fixed <- c("nborn", "muBday", "sdBday", "alphaBday")

# Fixed-wing, previous rates
fit_pas_prev <- fit_model(mod_deads, pas_dat$datalist, det_prev, pars_fixed,
                          "output/empirical/fit_PAS_Deads_prevRates.rds")

# Fixed-wing, updated rates
fit_pas_upd <- fit_model(mod_deads, pas_dat$datalist, det_PAS, pars_fixed,
                         "output/empirical/fit_PAS_Deads_updRates.rds")

# UAV, previous rates
fit_uav_prev <- fit_model(mod_deads, uav_dat$datalist, det_prev, pars_fixed,
                          "output/empirical/fit_UAV_Deads_prevRates.rds")

# UAV, updated rates
fit_uav_upd <- fit_model(mod_deads, uav_dat$datalist, det_UAV, pars_fixed,
                         "output/empirical/fit_UAV_Deads_updRates.rds")

# =============================================================================
# Objective 5: Estimated TTL 
# =============================================================================

pars_estTTL <- c("nborn", "muBday", "sdBday", "alphaBday", "muLeave", "sdLeave")

# UAV, updated rates
fit_uav_estTTL <- fit_model(mod_deads_estTTL, uav_dat$datalist, det_UAV, pars_estTTL,
                            "output/empirical/fit_UAV_Deads_estTTL.rds")

# Fixed-wing, updated rates
fit_pas_estTTL <- fit_model(mod_deads_estTTL, pas_dat$datalist, det_PAS, pars_estTTL,
                            "output/empirical/fit_PAS_Deads_estTTL.rds")

# =============================================================================
# Objective 5: Estimated TTM + TTL with force moult 
# =============================================================================

pars_estTTMTTL <- c("nborn", "muBday", "sdBday", "alphaBday",
                    "muLeave", "sdLeave", "muMoult", "sdMoult")

# Fixed-wing, updated rates
fit_pas_TTMTTL <- fit_model(mod_deads_estTTMTTL, pas_dat$datalist, det_PAS, pars_estTTMTTL,
                            "output/empirical/fit_PAS_Deads_estTTMTTL.rds")