start_time           <- Sys.time()
#=======================================================================================================================#
#                                               Control Variables:
#=======================================================================================================================#
# Validation Variable
check_Initial_Values <- FALSE
check_season         <- FALSE
check_births         <- FALSE
check_iGraph         <- FALSE
check_inputDistrib   <- TRUE
check_clustering     <- TRUE

# Save Data variables
write_output         <- TRUE
save_plots           <- TRUE
save_data            <- TRUE

# Experimental Variables
seasons_calc         <- TRUE
iterate              <- TRUE
cluster              <- TRUE
num_samples          <- 1700
k_min                <- 10
k_max                <- 30
pred_Introduce       <- FALSE
albs_Introduce       <- FALSE
midseason_pred       <- FALSE
midseason_albs       <- TRUE

#=======================================================================================================================#
#                                                   Libraries:
#=======================================================================================================================#
library(ggplot2)
library(deSolve)
library(lhs)
library(dplyr)
library(tidyr)
library(tibble)
library(cowplot)
library(progress)
library(ggh4x)
library(here)

#=======================================================================================================================#
#                                                Generate Variables:
#=======================================================================================================================#
generate_popvar          <- function(){
  popvar                       <- data.frame(total_gentoo = 2000)                         # Total Gentoo population [Source: Dr. Gamble]
  popvar$ResPer_gentoo         <- 1                                                       # Percent of Gentoo that are resident to the island [Source: 1]
  popvar$Res_gentoo            <- popvar$total_gentoo * popvar$ResPer_gentoo              # Calculate the number of resident Gentoo
  popvar$gentoo_migrate_in     <- popvar$total_gentoo * (1 - popvar$ResPer_gentoo)        # Calculate the number of migratory Gentoo
  popvar$gentoo_chicksurvival  <- 0.51                                                    # Breeding success rate [Source: 7]
  
  popvar$total_preds           <- 750                                                     # Total predator population [Source: Dr. Gamble]
  popvar$ResPer_preds          <- 0.1                                                     # Percent of predator that are resident to the island [Source: 1]
  popvar$Res_preds             <- popvar$total_preds * popvar$ResPer_preds                # Calculate the number of resident predators 
  popvar$preds_migrate_in      <- popvar$total_preds * (1 - popvar$ResPer_preds)          # Calculate the number of migratory predators
  popvar$preds_chicksurvival   <- 0.34                                                    # Breeding success rate [Source: 10]                                                 
  
  popvar$total_albatross       <- 50000                                                   # Total Albatross population [Source: Dr. Gamble]
  popvar$ResPer_albatross      <- 0                                                       # Percent of Albatross that are resident to the island [Source: 1]
  popvar$Res_albatross         <- popvar$total_albatross * popvar$ResPer_albatross        # Resident Albatross 
  popvar$albatross_migrate_in  <- popvar$total_albatross * (1 - popvar$ResPer_albatross)  # 100% of albatross are migratory
  popvar$albs_chicksurvival    <- 0.47                                                    # Breeding success rate [Source: 9]
  
  return(popvar)
}
generate_Intitial_values <- function(popvar){
  Initial_values <- c(
    S_Gentoo    = popvar$Res_gentoo,      # Susceptible
    I_Gentoo    = 0,                      # Infected
    R_Gentoo    = 0,                      # Recovered
    D_Gentoo    = 0,                      # Dead
    
    S_Predators = popvar$Res_preds,       # Susceptible
    I_Predators = 0,                      # Infected
    R_Predators = 0,                      # Recovered
    D_Predators = 0,                      # Dead
    
    S_Albatross = popvar$Res_albatross,   # Susceptible
    I_Albatross = 0,                      # Infected
    R_Albatross = 0,                      # Recovered
    D_Albatross = 0,                      # Dead
    
    timer       = 0
  )
  
  return(Initial_values)
}
generate_timevar         <- function(){
  timevar                   <- data.frame(t_arrival = 81) # Day at which migratory birds begin arriving                  [Source: 9]
  timevar$mig_day_arrival   <- 10                         # Number of days between arrival and 24 days before egg-laying [Source: 9]
  timevar$breed_start       <- 115                        # When egg-laying/incubation starts [Source: 7/8]
  timevar$hatch_start       <- 153                        # When hatching starts              [Source: 7/8]
  timevar$hatch_duration    <- 7                          # Hatching duration                 [Source: 8]
  timevar$t_departure       <- 300                        # 30 days post-fledging in April    [Source: 8]
  timevar$mig_day_departure <- 31
  timevar$horizon           <- 365
  
  return(timevar)
}
generate_muH             <- function(){
  min_Life_Gentoo       <- 10                                         # [Source: 2]
  max_Life_Gentoo       <- 15                                         # [Source: 2]
  avg_Life_Gentoo       <- (max_Life_Gentoo + min_Life_Gentoo) / 2    # Calculate average lifespan
  avg_Life_Gentoo_Days  <- avg_Life_Gentoo * 365                      # Calculate average lifespan in days
  muH_Gentoo            <- 1 / avg_Life_Gentoo_Days                   # Calculate MuH for the Gentoo
  
  min_Life_Preds        <- 0                                          # [Source: NA]
  max_Life_Preds        <- 28.8                                       # [Source: 4]
  avg_Life_Preds        <- (max_Life_Preds + min_Life_Preds) / 2      # Calculate average lifespan
  avg_Life_Preds_Days   <- avg_Life_Preds * 365                       # Calculate average lifespan in days
  muH_Preds             <- 1 / avg_Life_Preds_Days                    # Calculate MuH for the Predators
  
  min_Life_Albs         <- 0                                          # [Source: NA]
  max_Life_Albs         <- 47                                         # [Source: 3]
  avg_Life_Albs         <- (max_Life_Albs + min_Life_Albs) / 2        # Calculate average lifespan
  avg_Life_Albs_Days    <- avg_Life_Albs * 365                        # Calculate average lifespan in days
  muH_Albs              <- 1 / avg_Life_Albs_Days                     # Calculate MuH for the Albatross
  
  return(list(
    muH_Gentoo = muH_Gentoo, 
    muH_Preds  = muH_Preds, 
    muH_Albs   = muH_Albs))
}
generate_gamma_sigma     <- function(U){
  gentoo_inf_dur    <- qunif(U[,  5], min = 1,    max = 15)
  gentoo_inf_mort   <- qunif(U[,  6], min = 0.50, max = 0.95)
  muIs_gentoo       <- gentoo_inf_mort / gentoo_inf_dur
  gammas_gentoo     <- 1 / gentoo_inf_dur
  gentoo_imm_dur    <- 182.5                                      # Estimate without a basis
  #sigma_gentoo      <- 1 / gentoo_imm_dur                         # Rate of Gentoo immunity loss, >0 makes it an SIRS
  sigma_gentoo      <- 0                                          # Rate of Gentoo immunity loss, >0 makes it an SIRS
  
  preds_inf_dur     <- qunif(U[,  7], min = 1,    max = 15)
  preds_inf_mort    <- qunif(U[,  8], min = 0.01, max = 0.5)
  muIs_Preds        <- preds_inf_mort / preds_inf_dur
  gammas_preds      <- 1 / preds_inf_dur
  preds_imm_dur     <- 182.5                                      # Estimate without a basis
  #sigma_preds       <- 1 / preds_imm_dur                          # Rate of predator immunity loss, >0 makes it an SIRS
  sigma_preds       <- 0                                          # Rate of predator immunity loss, >0 makes it an SIRS
  
  Albs_inf_dur      <- qunif(U[,  9], min = 1,    max = 15)
  Albs_inf_mort     <- qunif(U[, 10], min = 0.50, max = 0.95)
  muIs_Albs         <- Albs_inf_mort / Albs_inf_dur
  gammas_albatross  <- 1 / Albs_inf_dur
  albatross_imm_dur <- 182.5                                      # Estimate without a basis
  #sigma_albatross   <- 1 / albatross_imm_dur                      # Rate of Albatross immunity loss, >0 makes it an SIRS
  sigma_albatross   <- 0                                          # Rate of Albatross immunity loss, >0 makes it an SIRS
  
  return(list(gammas_gentoo    = gammas_gentoo, 
              sigma_gentoo     = sigma_gentoo, 
              muIs_gentoo      = muIs_gentoo,
              gammas_preds     = gammas_preds,
              sigma_preds      = sigma_preds, 
              muIs_Preds       = muIs_Preds,
              gammas_albatross = gammas_albatross,
              sigma_albatross  = sigma_albatross,
              muIs_Albs        = muIs_Albs
  )
  )
}

#=======================================================================================================================#
#                                                Event Functions:
#=======================================================================================================================#
create_migin_eventdata  <- function(timevar, popvar){
  # Duration (days) for each species' arrival window
  L <- timevar$mig_day_arrival
  if (L <= 0) stop("timevar$mig_day_arrival must be > 0")
  
  # --- Arrival windows (sequential) ---
  # Albatross arrive first
  t_alb_start <- timevar$t_arrival
  t_alb_end   <- t_alb_start + L - 1
  alb_times   <- seq(from = t_alb_start, to = t_alb_end, length.out = L)
  #print(alb_times)
  
  # Predators arrive after albatrosses are done, for the same duration
  t_pred_start <- t_alb_end + 1
  t_pred_end   <- t_pred_start + L - 1
  pred_times   <- seq(from = t_pred_start, to = t_pred_end, length.out = L)
  
  # --- Per-day arrivals (constant per day over each window) ---
  pred_S_arrival <- (1 - popvar$prev_infection_preds) * popvar$preds_migrate_in     / L
  pred_I_arrival <- (    popvar$prev_infection_preds) * popvar$preds_migrate_in     / L
  alb_S_arrival  <- (1 - popvar$prev_infection_albs)  * popvar$albatross_migrate_in / L
  alb_I_arrival  <- (    popvar$prev_infection_albs)  * popvar$albatross_migrate_in / L
  
  # --- Build eventdata1 (sequential blocks) ---
  df_alb <- data.frame(
    var    = rep(c("S_Albatross", "I_Albatross"), each = length(alb_times)),
    time   = rep(alb_times, times = 2),
    value  = rep(c(alb_S_arrival, alb_I_arrival), each = length(alb_times)),
    method = "add"
  )
  
  df_pred <- data.frame(
    var    = rep(c("S_Predators", "I_Albatross"), each = length(pred_times)),
    time   = rep(pred_times, times = 2),
    value  = rep(c(pred_S_arrival, pred_I_arrival), each = length(pred_times)),
    method = "add"
  )
  
  df_newI <- data.frame(
    var    = c("I_Predators",      "I_Albatross",      "S_Predators",      "S_Albatross"),
    time   = c(popvar$ms_pred_day, popvar$ms_albs_day, popvar$ms_pred_day, popvar$ms_albs_day),
    value  = c(5,                  5,                  -5,                 -5),
    method = "add"
  )
  migineventdata      <- rbind(df_alb, df_pred, df_newI)
  migineventdata$time <- round(migineventdata$time, 3)
  migineventdata      <- migineventdata[order(migineventdata$time), ]
  
  return(migineventdata)
}
create_hatch_eventdata  <- function(timevar, popvar, N_Gent, N_Pred, N_Albs){
  # Duration (days) for each species' arrival window
  L <- timevar$hatch_duration
  if (L <= 0) stop("timevar$hatch_duration must be > 0")
  
  # --- Hatching window ---
  t_hatch_start <- timevar$hatch_start
  t_hatch_end   <- t_hatch_start + L - 1
  hatch_times   <- seq(from = t_hatch_start, to = t_hatch_end, length.out = L)
  
  # --- Per-day hatches (constant per day over each window) ---
  gent_S_hatch <- (N_Gent / 2) * popvar$gentoo_chicksurvival / L
  pred_S_hatch <- (N_Pred / 2) * popvar$preds_chicksurvival  / L
  albs_S_hatch <- (N_Albs / 2) * popvar$albs_chicksurvival   / L
  
  # --- Build hatcheventdata ---
  hatcheventdata <- data.frame(
    var    = rep(c("S_Gentoo","S_Predators", "S_Albatross"),  each = length(hatch_times)),
    time   = rep(hatch_times, times = 3),
    value  = rep(c(gent_S_hatch, pred_S_hatch, albs_S_hatch), each = length(hatch_times)),
    method = "add"
  )
  
  df_newI <- data.frame(
    var    = c("I_Predators",      "I_Albatross",      "S_Predators",      "S_Albatross"),
    time   = c(popvar$ms_pred_day, popvar$ms_albs_day, popvar$ms_pred_day, popvar$ms_albs_day),
    value  = c(5,                  5,                  -5,                 -5),
    method = "add"
  )
  
  hatcheventdata      <- rbind(hatcheventdata, df_newI)
  hatcheventdata$time <- round(hatcheventdata$time, 3)
  hatcheventdata      <- hatcheventdata[order(hatcheventdata$time), ]
  
  return(hatcheventdata)
}
create_migout_eventdata <- function(timevar, popvar, statedeparturevector){
  # Duration (days) for each species' arrival window
  L  <- timevar$mig_day_departure
  if (L <= 0) stop("timevar$t_departure must be > 0")
  
  # --- Hatching window ---
  t_migout_start <- timevar$t_departure
  t_migout_end   <- t_migout_start + L - 1
  migout_times   <- seq(from = t_migout_start, to = t_migout_end, length.out = L)
  
  S_Gentoo    <- statedeparturevector["S_Gentoo"]
  S_Predators <- statedeparturevector["S_Predators"]
  S_Albatross <- statedeparturevector["S_Albatross"]
  
  I_Gentoo    <- statedeparturevector["I_Gentoo"]
  I_Predators <- statedeparturevector["I_Predators"]
  I_Albatross <- statedeparturevector["I_Albatross"]
  
  R_Gentoo    <- statedeparturevector["R_Gentoo"]
  R_Predators <- statedeparturevector["R_Predators"]
  R_Albatross <- statedeparturevector["R_Albatross"]
  
  S_Gentoo_Migout    <- S_Gentoo    * (1 - popvar$ResPer_gentoo)    / L
  S_Predators_Migout <- S_Predators * (1 - popvar$ResPer_preds)     / L
  S_Albatross_Migout <- S_Albatross * (1 - popvar$ResPer_albatross) / L
  
  I_Gentoo_Migout    <- I_Gentoo    * (1 - popvar$ResPer_gentoo)    / L
  I_Predators_Migout <- I_Predators * (1 - popvar$ResPer_preds)     / L
  I_Albatross_Migout <- I_Albatross * (1 - popvar$ResPer_albatross) / L
  
  
  R_Gentoo_Migout    <- R_Gentoo    * (1 - popvar$ResPer_gentoo)    / L
  R_Predators_Migout <- R_Predators * (1 - popvar$ResPer_preds)     / L
  R_Albatross_Migout <- R_Albatross * (1 - popvar$ResPer_albatross) / L
  
  migout_eventdata <- data.frame(
    var    = rep(c("S_Gentoo", "S_Predators", "S_Albatross", 
                   "I_Gentoo", "I_Predators", "I_Albatross",
                   "R_Gentoo", "R_Predators", "R_Albatross"),  
                 each = length(migout_times)),
    time   = rep(migout_times, times = 9),
    value  = rep(c(- S_Gentoo_Migout, - S_Predators_Migout, - S_Albatross_Migout, 
                   - I_Gentoo_Migout, - I_Predators_Migout, - I_Albatross_Migout,
                   - R_Gentoo_Migout, - R_Predators_Migout, - R_Albatross_Migout), 
                 each = length(migout_times)),
    method = "add"
  )
  
  migout_eventdata$time <- round(migout_eventdata$time, 3)
  migout_eventdata      <- migout_eventdata[order(migout_eventdata$time), ]
  
  return(migout_eventdata)
}

#=======================================================================================================================#
#                                                  Model Functions:
#=======================================================================================================================#
seasonality_transmission_calc <- function(time, t_arrival, mig_day_arrival, breed_start, 
                                          hatch_start, t_departure){
  if(time < (t_arrival + (mig_day_arrival / 2))){
    1 # Before arrival
  } else{
    if(time < breed_start){
      2 # During Arrival/Mating
    } else{
      if(time < hatch_start){
        1 # After eggs are laid
      } else{
        if(time < t_departure){
          2 # After eggs are hatched
        } else{
          1 # After migration out
        }
      }
    }
  }
}
make_beta                     <- function(gentgent_trans, gentpred_trans, 
                                          predgent_trans, predpred_trans, predalbs_trans, 
                                          albspred_trans, albsalbs_trans){
  species = c("gentoo", "predators", "albatross")
  beta <- matrix(
    c(
      gentgent_trans, gentpred_trans, 0,
      predgent_trans, predpred_trans, predalbs_trans,
      0,              albspred_trans, albsalbs_trans
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(species, species)
  )
  beta
}
calc_R0                       <- function(beta, S0, gamma, muI){
  removal = gamma + muI
  if(any(removal <= 0)) stop("All removal (gamma + muI) must be > 0")
  
  K <- sweep(beta, 2, removal, "/")
  K <- sweep(K, 1, S0, "*")
  
  ev <- eigen(K, only.values = TRUE)$values
  list(R0 = max(Mod(ev)), eigenvalues = ev)
}
HPAI_dyn                      <- function(time, state, parameters){
  # Extract state
  S_Gentoo        <- state["S_Gentoo"]
  I_Gentoo        <- state["I_Gentoo"]
  R_Gentoo        <- state["R_Gentoo"]
  D_Gentoo        <- state["D_Gentoo"]
  
  S_Predators     <- state["S_Predators"]
  I_Predators     <- state["I_Predators"]
  R_Predators     <- state["R_Predators"]
  D_Predators     <- state["D_Predators"]
  
  S_Albatross     <- state["S_Albatross"]
  I_Albatross     <- state["I_Albatross"]
  R_Albatross     <- state["R_Albatross"]
  D_Albatross     <- state["D_Albatross"]
  
  timer           <- state["timer"]
  
  # Extract parameters [Time]
  seasons_calc     = parameters["seasons_calc"]
  t_arrival        = parameters["t_arrival"] 
  mig_day_arrival  = parameters["mig_day_arrival"]
  breed_start      = parameters["breed_start"] 
  hatch_start      = parameters["hatch_start"]
  t_departure      = parameters["t_departure"]
  
  # Extract parameters [Mortality Rates (Healthy (H) vs. Infected (I))]
  muH_Gentoo      <- parameters["muH_Gentoo"]
  muH_Preds       <- parameters["muH_Preds"]
  muH_Albs        <- parameters["muH_Albs"]
  muI_Gentoo      <- parameters["muI_Gentoo"]
  muI_Preds       <- parameters["muI_Preds"]
  muI_Albs        <- parameters["muI_Albs"]
  
  
  # Extract parameters [Recovery Rates (Gamma) and Immunity Loss Rates (Sigma)]
  gamma_gentoo    <- parameters["gamma_gentoo"]
  sigma_gentoo    <- parameters["sigma_gentoo"]
  gamma_preds     <- parameters["gamma_preds"]
  sigma_preds     <- parameters["sigma_preds"]
  gamma_albatross <- parameters["gamma_albatross"]
  sigma_albatross <- parameters["sigma_albatross"]
  
  # Extract parameters [Betas]
  gentgent_trans  <- parameters["gentgent_trans"]
  gentpred_trans  <- parameters["gentpred_trans"]
  predgent_trans  <- parameters["predgent_trans"]
  predpred_trans  <- parameters["predpred_trans"]
  predalbs_trans  <- parameters["predalbs_trans"]
  albspred_trans  <- parameters["albspred_trans"]
  albsalbs_trans  <- parameters["albsalbs_trans"]
  
  # guard names exist
  need_states <- c("S_Gentoo","I_Gentoo","R_Gentoo","D_Gentoo",
                   "S_Predators","I_Predators","R_Predators","D_Predators",
                   "S_Albatross","I_Albatross","R_Albatross","D_Albatross","timer")
  need_parms  <- c("t_arrival","mig_day_arrival","breed_start","hatch_start","t_departure", "seasons_calc",
                   "muH_Gentoo","muH_Preds","muH_Albs",
                   "muI_Gentoo","muI_Preds","muI_Albs",
                   "gamma_gentoo","sigma_gentoo","gamma_preds","sigma_preds","gamma_albatross","sigma_albatross",
                   "gentgent_trans","gentpred_trans","predgent_trans","predpred_trans",
                   "predalbs_trans","albspred_trans","albsalbs_trans")
  if (!all(need_states %in% names(state))) stop("Missing state(s): ", paste(setdiff(need_states,names(state)), collapse=","))
  if (!all(need_parms %in% names(parameters))) stop("Missing param(s): ", paste(setdiff(need_parms,names(parameters)), collapse=","))
  
  # optional: clamp tiny negatives
  if(S_Gentoo    < 0 || I_Gentoo < 0    || R_Gentoo < 0 ||
     S_Predators < 0 || I_Predators < 0 || R_Predators < 0 ||
     S_Albatross < 0 || I_Albatross < 0 || R_Albatross < 0){
    S_Gentoo    <- max(S_Gentoo, 0); I_Gentoo    <- max(I_Gentoo, 0)
    S_Predators <- max(S_Predators,0); I_Predators<- max(I_Predators,0)
    S_Albatross <- max(S_Albatross,0); I_Albatross<- max(I_Albatross,0)
  }
  
  # Set the seasonality modifier
  if(seasons_calc == TRUE){
    seasMod         <- seasonality_transmission_calc(time, t_arrival, mig_day_arrival, breed_start, 
                                                     hatch_start, t_departure)
  } else{
    seasMod <- 1
  }
  
  #===============#
  #    Gentoo     #
  #===============#
  dS_Gentoo = 
    - muH_Gentoo      * S_Gentoo -                   
    gentgent_trans  * S_Gentoo    * I_Gentoo    * seasMod -     
    predgent_trans  * S_Gentoo    * I_Predators * seasMod +
    sigma_gentoo    * R_Gentoo
  
  dI_Gentoo = 
    - muI_Gentoo      * I_Gentoo +                  
    gentgent_trans  * S_Gentoo    * I_Gentoo    * seasMod +
    predgent_trans  * S_Gentoo    * I_Predators * seasMod -
    gamma_gentoo    * I_Gentoo
  
  dR_Gentoo = 
    + gamma_gentoo    * I_Gentoo -
    sigma_gentoo    * R_Gentoo
  
  dD_Gentoo = 
    muI_Gentoo        * I_Gentoo
  
  #===============#
  #   Predators   #
  #===============#
  dS_Predators = 
    - muH_Preds       * S_Predators -              
    gentpred_trans  * S_Predators * I_Gentoo    * seasMod -
    predpred_trans  * S_Predators * I_Predators * seasMod -
    albspred_trans  * S_Predators * I_Albatross * seasMod +
    sigma_preds     * R_Predators
  
  dI_Predators = 
    - muI_Preds       * I_Predators +              
    gentpred_trans  * S_Predators * I_Gentoo    * seasMod +
    predpred_trans  * S_Predators * I_Predators * seasMod +
    albspred_trans  * S_Predators * I_Albatross * seasMod -
    gamma_preds     * I_Predators
  
  dR_Predators = 
    + gamma_preds     * I_Predators -
    sigma_preds     * R_Predators
  
  dD_Predators = 
    muI_Preds         * I_Predators
  
  #===============#
  #   Albatross   #
  #===============#
  if (time < t_departure){
    mort_S_alb <- muH_Albs * S_Albatross
    mort_I_alb <- muI_Albs * I_Albatross
  } else{
    mort_S_alb <- 0
    mort_I_alb <- 0
  } # Do not count mortality after migration out starts
  dS_Albatross = 
    - mort_S_alb -
    predalbs_trans  * S_Albatross * I_Predators * seasMod -
    albsalbs_trans  * S_Albatross * I_Albatross * seasMod +
    sigma_albatross * R_Albatross
  
  dI_Albatross = 
    - mort_I_alb +
    predalbs_trans  * S_Albatross * I_Predators * seasMod +
    albsalbs_trans  * S_Albatross * I_Albatross * seasMod -
    gamma_albatross * I_Albatross
  
  dR_Albatross = 
    + gamma_albatross * I_Albatross -
    sigma_albatross * R_Albatross
  
  dD_Albatross =
    mort_I_alb
  
  #===============#
  #    Timer      #
  #===============#
  dTimer = 1
  
  #Clean up and return values
  derivs <- c(
    dS_Gentoo,    dI_Gentoo,    dR_Gentoo,    dD_Gentoo,
    dS_Predators, dI_Predators, dR_Predators, dD_Predators,
    dS_Albatross, dI_Albatross, dR_Albatross, dD_Albatross,
    dTimer
  )
  
  # strip any accidental names
  names(derivs) <- NULL
  
  list(derivs)
}
run_Model                     <- function(timevar, popvar, Initial_values, parameters){
  # Set the three distinct time pieces that will be evaluated
  time1 <- seq(from = 0, to = timevar$hatch_start, by = 1)
  time2 <- seq(from = timevar$hatch_start, to = timevar$t_departure, by = 1)
  time3 <- seq(from = timevar$t_departure, to = timevar$horizon, by = 1)
  
  # Set the migration in eventdata
  migin_eventdata <- create_migin_eventdata(timevar, popvar)
  
  # Run the first chunk of the ODE (Model start -> Hatching)
  solved1 <- ode(
    y      = Initial_values,
    func   = HPAI_dyn,
    parms  = parameters,
    times  = time1,
    events = list(data = migin_eventdata),
    method = "lsoda"
  )
  
  # Clean data, then prepare for the next evolution
  solved1 <- as.data.frame(solved1)
  statedeparturevector<- as.numeric(solved1[ nrow(solved1), -1])
  names(statedeparturevector) <- names(Initial_values)
  
  # Set the hatching eventdata
  N_Gent <- as.numeric(statedeparturevector["S_Gentoo"])    + 
    as.numeric(statedeparturevector["I_Gentoo"]) +
    as.numeric(statedeparturevector["R_Gentoo"])
  N_Pred <- as.numeric(statedeparturevector["S_Predators"]) + 
    as.numeric(statedeparturevector["I_Predators"]) +
    as.numeric(statedeparturevector["R_Predators"])
  N_Albs <- as.numeric(statedeparturevector["S_Albatross"]) + 
    as.numeric(statedeparturevector["I_Albatross"]) +
    as.numeric(statedeparturevector["R_Albatross"])
  hatcheventdata <- create_hatch_eventdata(timevar, popvar, N_Gent, N_Pred, N_Albs)
  
  # Run the second chunk of the ODE (Hatching -> Departure)
  solved2 <- ode(
    y      = statedeparturevector,
    func   = HPAI_dyn,
    parms  = parameters,
    times  = time2,
    events = list(data = hatcheventdata), 
    method = "lsoda"
  )
  
  # Clean data then prepare for the next evolution
  solved2 <- as.data.frame(solved2)
  statedeparturevector<- as.numeric(solved2[ nrow(solved2), -1])
  names(statedeparturevector) <- names(Initial_values)
  
  # Set the hatching eventdata
  N_Gent <- as.numeric(statedeparturevector["S_Gentoo"])    + 
    as.numeric(statedeparturevector["I_Gentoo"]) +
    as.numeric(statedeparturevector["R_Gentoo"])
  N_Pred <- as.numeric(statedeparturevector["S_Predators"]) + 
    as.numeric(statedeparturevector["I_Predators"]) +
    as.numeric(statedeparturevector["R_Predators"])
  N_Albs <- as.numeric(statedeparturevector["S_Albatross"]) + 
    as.numeric(statedeparturevector["I_Albatross"]) +
    as.numeric(statedeparturevector["R_Albatross"])
  migout_eventdata <- create_migout_eventdata(timevar, popvar, statedeparturevector)
  
  # Run the second chunk of the ODE (Departure -> End of Model)
  solved3 <- ode(
    y      = statedeparturevector,
    func   = HPAI_dyn,
    parms  = parameters,
    times  = time3,
    events = list(data = migout_eventdata),
    method = "lsoda"
  )
  
  # Clean data
  solved3 <- as.data.frame(solved3)
  
  # Clean up variables, and bind the three dataframes together
  out <- NULL
  out <- rbind(solved1, solved2, solved3)
  names(out) <- c(
    "Time",
    "Susceptible_Gentoo",    "Infected_Gentoo",    "Recovered_Gentoo",    "HPAI_Killed_Gentoo",
    "Susceptible_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
    "Susceptible_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross",
    "Timer"
  )
  out <- out[!duplicated(out$Time),]
  
  # choose the five interspecies series you want to store
  trans_vars <- c("inc_gentpred",   # Gentoo → Predators
                  "inc_predgent",   # Predators → Gentoo
                  "inc_predalbs",   # Predators → Albatross
                  "inc_albspred",   # Albatross → Predators
                  "inc_albsalbs")   # Albatross → Albatross  (replace if you meant another)
  
  # Connect pre- and post- migration outs and then calculate cross-species transmission events
  inc_gentgent <- parameters["gentgent_trans"] * out$Susceptible_Gentoo    * out$Infected_Gentoo
  inc_gentpred <- parameters["gentpred_trans"] * out$Susceptible_Predators * out$Infected_Gentoo
  inc_predgent <- parameters["predgent_trans"] * out$Susceptible_Gentoo    * out$Infected_Predators
  inc_predpred <- parameters["predpred_trans"] * out$Susceptible_Predators * out$Infected_Predators
  inc_predalbs <- parameters["predalbs_trans"] * out$Susceptible_Albatross * out$Infected_Predators
  inc_albspred <- parameters["albspred_trans"] * out$Susceptible_Predators * out$Infected_Albatross
  inc_albsalbs <- parameters["albsalbs_trans"] * out$Susceptible_Albatross * out$Infected_Albatross
  
  # Create a dataframe of the series
  run_tbl <- data.frame(
    Time          = out$Time,
    inc_gentpred  = inc_gentpred,
    inc_predgent  = inc_predgent,
    inc_predalbs  = inc_predalbs,
    inc_albspred  = inc_albspred,
    inc_albsalbs  = inc_albsalbs
  )
  
  # Convert to long format
  run_long <- pivot_longer(
    run_tbl,
    cols      = all_of(trans_vars),
    names_to  = "series",
    values_to = "value"
  )
  
  # Calculate the cumulative infections for each transmission type, compartment, and in total
  dt     <- diff(out$Time)[1]
  cum_gg <- sum(inc_gentgent) * dt
  cum_gp <- sum(inc_gentpred) * dt
  cum_pg <- sum(inc_predgent) * dt
  cum_pp <- sum(inc_predpred) * dt
  cum_pa <- sum(inc_predalbs) * dt
  cum_ap <- sum(inc_albspred) * dt
  cum_aa <- sum(inc_albsalbs) * dt
  CI_Gentoo    <- cum_gg + cum_pg
  CI_Predators <- cum_gp + cum_pp + cum_ap
  CI_Albatross <- cum_pa + cum_aa
  CI_Total     <- CI_Gentoo + CI_Predators + CI_Albatross
  
  # Calculate R0 (Next Generation Matrix Method)
  if(popvar$prev_infection_albs != 0  || popvar$prev_infection_preds != 0 || 
     popvar$ms_pred_day         != -1 || popvar$ms_albs_day          != -1){
    beta <- make_beta(
      parameters["gentgent_trans"], parameters["gentpred_trans"],
      parameters["predgent_trans"], parameters["predpred_trans"], parameters["predalbs_trans"],
      parameters["albspred_trans"], parameters["albsalbs_trans"])
    if(popvar$prev_infection_albs != 0){
      t_eval = timevar$t_arrival + (2 * timevar$mig_day_arrival)
    } else if(popvar$prev_infection_preds != 0){
      t_eval = timevar$t_arrival +      timevar$mig_day_arrival
    } else if(popvar$ms_pred_day != -1){
      t_eval = popvar$ms_pred_day + 2
    } else if(popvar$ms_albs_day != -1){
      t_eval = popvar$ms_albs_day + 2
    }
    S0 <- c(gentoo    = out[ t_eval, "Susceptible_Gentoo"],
            predators = out[ t_eval, "Susceptible_Predators"],
            albatross = out[ t_eval, "Susceptible_Albatross"])
    gamma <- c(parameters["gamma_gentoo"], parameters["gamma_preds"], parameters["gamma_albatross"])
    muI   <- c(parameters["muI_Gentoo"], parameters["muI_Preds"], parameters["muI_Albs"])
    eigen <- calc_R0(beta, S0, gamma, muI)
    R0    <- eigen$R0
  } else{
    R0 <- NA
  }
  
  # Create the return variable with all the data of interest
  df                     <- data.frame(prev_infection_preds = popvar$prev_infection_preds)
  df$prev_infection_albs <- popvar$prev_infection_albs
  df$midseason_pred_day  <- popvar$ms_pred_day
  df$midseason_albs_day  <- popvar$ms_albs_day
  df$gamma_gentoo        <- parameters["gamma_gentoo"]
  df$muI_gentoo          <- parameters["muI_Gentoo"]
  df$gamma_predators     <- parameters["gamma_preds"]
  df$muI_predators       <- parameters["muI_Preds"]
  df$gamma_albatross     <- parameters["gamma_albatross"]
  df$muI_albatross       <- parameters["muI_Albs"]
  df$Beta_GentGent       <- parameters["gentgent_trans"]
  df$Beta_GentPred       <- parameters["gentpred_trans"]
  df$Beta_PredGent       <- parameters["predgent_trans"]
  df$Beta_PredPred       <- parameters["predpred_trans"]
  df$Beta_PredAlbs       <- parameters["predalbs_trans"]
  df$Beta_AlbsPred       <- parameters["albspred_trans"]
  df$Beta_AlbsAlbs       <- parameters["albsalbs_trans"]
  df$R0                  <- R0
  df$maxI_Gentoo         <- max(out["Infected_Gentoo"])
  df$time_maxI_Gentoo    <- out[ which.max(out[["Infected_Gentoo"]]), "Time"]
  df$maxI_Predators      <- max(out["Infected_Predators"])
  df$time_maxI_Predators <- out[ which.max(out[["Infected_Predators"]]), "Time"]
  df$maxI_Albatross      <- max(out["Infected_Albatross"])
  df$time_maxI_Albatross <- out[ which.max(out[["Infected_Albatross"]]), "Time"]
  df$terI_Gentoo         <- round(out[ nrow(out), "Infected_Gentoo"])
  df$terI_Predators      <- round(out[ nrow(out), "Infected_Predators"])
  df$overwintering       <- ifelse(df$terI_Gentoo >= 1 || df$terI_Predators >= 1, 1, 0)
  df$outbreakdur         <- length(which(out$Infected_Gentoo + out$Infected_Predators + out$Infected_Albatross > 2))
  df$outbreakspan        <- ifelse(df$outbreakdur > (timevar$t_departure - timevar$t_arrival - 5), 1, 0)
  df$infperiod_Gentoo    <- length(which(out$Infected_Gentoo    > 1))
  df$infperiod_Predators <- length(which(out$Infected_Predators > 1))
  df$infperiod_Albatross <- length(which(out$Infected_Albatross > 1))
  df$deadGentoo          <- out[ nrow(out), "HPAI_Killed_Gentoo"]
  df$deadPredators       <- out[ nrow(out), "HPAI_Killed_Predators"]
  df$deadAlbatross       <- out[ nrow(out), "HPAI_Killed_Albatross"]
  df$gentpred_transevents<- cum_gp;       df$predgent_transevents<- cum_pg
  df$predalbs_transevents<- cum_pa;       df$albspred_transevents<- cum_ap
  df$cuminf_Total        <- CI_Total;     df$cuminf_Gentoo       <- CI_Gentoo
  df$cuminf_Predators    <- CI_Predators; df$cuminf_Albatross    <- CI_Albatross
  
  # Return the dataframe
  return(list(df       = df, 
              values   = out, 
              run_long = run_long))
}
build_histogram_features      <- function(valarray, bins = 52, popvar = NULL, 
                                          relativize_to_species_totals = TRUE){
  # Note: If this variable is changed, the species_for_comp variable must be too
  compartments_use = c(
    "Infected_Gentoo", "Infected_Predators", "Infected_Albatross",
    "HPAI_Killed_Gentoo", "HPAI_Killed_Predators", "HPAI_Killed_Albatross"
  )
  dn <- dimnames(valarray)
  if (is.null(dn) || is.null(dn$compartment)) {
    stop("valarray must have dimnames with a 'compartment' entry.")
  }
  if (!all(compartments_use %in% dn$compartment)) {
    missing <- setdiff(compartments_use, dn$compartment)
    stop("Missing compartments in valarray: ", paste(missing, collapse = ", "))
  }
  
  C <- length(compartments_use)
  T <- dim(valarray)[2]
  N <- dim(valarray)[3]
  
  # Optional per-capita relativization
  if (isTRUE(relativize_to_species_totals)) {
    if (is.null(popvar)) stop("relativize_to_species_totals=TRUE requires popvar.")
    total_gentoo    <- popvar$total_gentoo    * (1 + popvar$gentoo_chicksurvival)
    total_predators <- popvar$total_preds * (1 + popvar$preds_chicksurvival)
    total_albatross <- popvar$total_albatross * (1 + popvar$albs_chicksurvival)
    if (any(!is.finite(c(total_gentoo, total_predators, total_albatross))) ||
        any(c(total_gentoo, total_predators, total_albatross) <= 0)) {
      stop("Species totals in popvar must be positive, finite numbers.")
    }
    # Note: If the compartments_use variable is changed, this one must be too
    species_for_comp <- c("Gentoo", "Predators", "Albatross",
                          "Gentoo", "Predators", "Albatross")  
    totals <- list(Gentoo=total_gentoo, Predators=total_predators, Albatross=total_albatross)
    for (j in seq_along(compartments_use)) {
      comp <- compartments_use[j]
      denom <- totals[[ species_for_comp[j] ]]
      valarray[comp, , ] <- valarray[comp, , ] / denom
    }
  }
  
  # Shared breaks per compartment
  comp_breaks <- vector("list", C)
  names(comp_breaks) <- compartments_use
  for (j in seq_along(compartments_use)) {
    comp <- compartments_use[j]
    x_all <- as.numeric(valarray[comp, , ])
    x_all <- x_all[is.finite(x_all)]
    if (!length(x_all)) stop("All values non-finite for compartment: ", comp)
    rng <- range(x_all, na.rm = TRUE)
    if (rng[1] == rng[2]) {
      eps <- if (rng[1] == 0) 1e-8 else abs(rng[1]) * 1e-6
      rng <- rng + c(-eps, +eps)
    }
    comp_breaks[[j]] <- seq(rng[1], rng[2], length.out = bins + 1)
  }
  
  # Build feature matrix (density histograms)
  feature_mat <- matrix(NA_real_, nrow = dim(valarray)[3], ncol = 0)
  rn <- dimnames(valarray)$sample %||% paste0("Sample_", seq_len(dim(valarray)[3]))
  rownames(feature_mat) <- rn
  
  # Bin the values from the compartments
  for (j in seq_along(compartments_use)) {
    comp <- compartments_use[j]
    brks <- comp_breaks[[j]]
    block <- matrix(NA_real_, nrow = length(rn), ncol = bins)
    colnames(block) <- paste0(comp, "_bin", seq_len(bins))
    for (i in seq_along(rn)) {
      ts_ij <- as.numeric(valarray[comp, , i])
      h <- hist(ts_ij, breaks = brks, plot = FALSE, right = TRUE, include.lowest = TRUE)
      dens <- h$density
      if (any(!is.finite(dens))) {
        cnt <- h$counts; s <- sum(cnt); dens <- if (s > 0) cnt / s else rep(0, bins)
      }
      block[i, ] <- dens
    }
    feature_mat <- cbind(feature_mat, block)
  }
  
  # Drop zero-variance features (across runs)
  keep_cols <- which(apply(feature_mat, 2, function(v) stats::sd(v, na.rm=TRUE) > 0))
  feature_mat <- feature_mat[, keep_cols, drop = FALSE]
  
  # Return the relevant variables
  list(
    features         = feature_mat,
    breaks           = comp_breaks,
    compartments_use = compartments_use
  )
}
evaluate_kmeans_k             <- function(feature_mat, k_min = 2, k_max = 20, nstart = 50, 
                                          iter.max = 100, seed = 1L){
  if(nrow(feature_mat) <= k_max){
    stop("Number of clusters must be less than the number of samples.")
  }
  if(nrow(feature_mat) <= 2 * k_max){
    warning("Number of clusters should be less than half the number of samples.")
  }
  set.seed(seed)
  
  X <- feature_mat
  N <- nrow(X)
  p <- ncol(X)
  
  # Precompute pairwise distances once for silhouette; we recompute within each k for cluster labels
  D <- as.matrix(dist(X, method = "euclidean"))
  
  # total sum of squares (about the grand mean)
  grand_mean <- colMeans(X)
  TSS <- sum(rowSums((X - matrix(grand_mean, nrow=N, ncol=p, byrow=TRUE))^2))
  
  ks <- k_min:k_max
  out <- data.frame(
    k          = ks,
    WSS        = NA_real_, # within-cluster sum of squares
    CH         = NA_real_, # Calinski–Harabasz (higher is better)
    DB         = NA_real_, # Davies–Bouldin (lower is better)
    Silhouette = NA_real_  # mean silhouette width (higher is better)
  )
  
  # helper: compute cluster stats used for CH and DB
  cluster_stats <- function(X, cl, centers) {
    k <- length(unique(cl))
    # WSS
    WSS <- 0
    S_i <- numeric(k)     # average within-cluster distance to center (for DB)
    for (j in seq_len(k)) {
      idx <- which(cl == j)
      if (length(idx) == 0) next
      Cj <- centers[j, , drop = FALSE]
      diffs <- X[idx, , drop = FALSE] - matrix(Cj, nrow=length(idx), ncol=ncol(X), byrow=TRUE)
      sq <- rowSums(diffs^2)
      WSS <- WSS + sum(sq)
      # S_i: mean Euclidean distance to center
      S_i[j] <- mean(sqrt(sq))
    }
    list(WSS = WSS, S_i = S_i)
  }
  
  # helper: Davies–Bouldin
  davies_bouldin <- function(centers, S_i) {
    k <- nrow(centers)
    if (k == 1) return(NA_real_)
    # inter-center distances
    Cdist <- as.matrix(dist(centers))
    # for each i, compute max_j ( (S_i + S_j) / d(c_i,c_j) )
    R_i <- numeric(k)
    for (i in seq_len(k)) {
      vals <- (S_i[i] + S_i[-i]) / Cdist[i, -i]
      R_i[i] <- max(vals[is.finite(vals)])
    }
    mean(R_i)
  }
  
  mean_silhouette <- function(D, cl) {
    # D: full NxN distance matrix (Euclidean)
    # cl: integer cluster labels (length N)
    labs <- sort(unique(cl))
    k <- length(labs)
    if (k < 2) stop("Mean_silhouette cannot be calculated if k < 2")  # silhouette undefined for k=1
    
    N <- nrow(D)
    a <- numeric(N)  # mean intra-cluster distance (exclude self)
    b <- numeric(N)  # mean distance to the nearest *other* cluster
    
    for (lab in labs) {
      idx <- which(cl == lab)
      # a_i: mean distance to points in own cluster (exclude self)
      if (length(idx) > 1) {
        Di <- D[idx, idx, drop = FALSE]
        diag(Di) <- NA_real_                   # exclude self
        a[idx] <- rowMeans(Di, na.rm = TRUE)
      } else {
        a[idx] <- 0                            # singleton cluster: define a_i = 0
      }
      
      # b_i: for each other cluster, the mean distance to that cluster; take min
      other_labs <- setdiff(labs, lab)
      if (length(other_labs) == 0) {
        b[idx] <- 0
      } else {
        other_means <- lapply(other_labs, function(l2) {
          jdx <- which(cl == l2)
          rowMeans(D[idx, jdx, drop = FALSE], na.rm = TRUE)
        })
        Bmat <- do.call(cbind, other_means)
        if (is.null(dim(Bmat))) {
          # Only one other cluster -> Bmat is a vector
          b[idx] <- as.numeric(Bmat)
        } else {
          b[idx] <- apply(Bmat, 1, min)
        }
      }
    }
    
    s <- (b - a) / pmax(a, b)
    mean(s[is.finite(s)])
  }
  

  for (idx in seq_along(ks)) {
    k <- ks[idx]
    km <- kmeans(X, centers = k, nstart = nstart, iter.max = iter.max)
    stats <- cluster_stats(X, km$cluster, km$centers)
    WSS <- stats$WSS
    # Calinski–Harabasz
    CH <- ((TSS - WSS) / (k - 1)) / (WSS / (N - k))
    # Davies–Bouldin
    DB <- davies_bouldin(km$centers, stats$S_i)
    # Silhouette (mean)
    Sil <- mean_silhouette(D, km$cluster)
    out[idx, c("WSS","CH","DB","Silhouette")] <- c(WSS, CH, DB, Sil)
  }
  
  # Simple recommendations by criterion:
  # - CH: maximize
  # - Silhouette: maximize
  # - DB: minimize
  rec <- list()
  rec$k_CH  <- out$k[ which.max(out$CH) ]
  rec$k_Sil <- out$k[ which.max(out$Silhouette) ]
  rec$k_DB  <- out$k[ which.min(out$DB) ]

  list(metrics = out, recommended_k = rec, scaled_features = X)
}
fit_kmeans_at_k               <- function(X_scaled, k, nstart = 100, iter.max = 200, seed = 1L) {
  set.seed(seed)
  kmeans(X_scaled, centers = k, nstart = nstart, iter.max = iter.max)
}

#=======================================================================================================================#
#                                                 Generate Graphs:
#=======================================================================================================================#
vardistplots                  <- function(dfMaster, cnum){
  inputdf <- dfMaster %>%
    transmute(
      prev_infection_preds = ifelse(prev_infection_preds != 0, prev_infection_preds, NA),
      prev_infection_albs  = ifelse(prev_infection_albs  != 0, prev_infection_albs,  NA),
      midseason_pred_day   = ifelse(midseason_pred_day   != -1, midseason_pred_day,   NA),
      midseason_albs_day   = ifelse(midseason_albs_day   != -1, midseason_albs_day,   NA),
      
      g_InfDur  = 1 / gamma_gentoo,                    
      g_InfMort = (1 / gamma_gentoo) * muI_gentoo,
      
      p_InfDur  = 1 / gamma_predators,
      p_InfMort = (1 / gamma_predators) * muI_predators,
      
      a_InfDur  = 1 / gamma_albatross,
      a_InfMort = (1 / gamma_albatross) * muI_albatross,
      
      B_GentGent = Beta_GentGent,
      B_GentPred = Beta_GentPred,
      B_PredGent = Beta_PredGent,
      B_PredPred = Beta_PredPred,
      B_PredAlbs = Beta_PredAlbs,
      B_AlbsPred = Beta_AlbsPred,
      B_AlbsAlbs = Beta_AlbsAlbs,
      
      R0 = R0,
      clusternum = Cluster
    )
  
  df_long <- inputdf %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(
      cols = -c(row_id, clusternum),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(value = suppressWarnings(as.numeric(value))) %>%
    filter(is.finite(value))
  
  df_red <- df_long %>% filter(clusternum == cnum)
  
  # ---- plot: constant x so each facet's box fills the width ----
  ggplot(df_long, aes(x = 0, y = value)) +
    geom_violin(
      width = 0.9, trim = FALSE,
      fill = "lightblue", color = "darkblue"
    ) +
    geom_jitter(
      data  = df_red,
      aes(x = 0, y = value),
      width = 0.25, height = 0,
      color = "red", size = 1.5, alpha = 0.8
    ) +
    facet_wrap(~ variable, scales = "free_y", ncol = 4) +
    scale_x_continuous(breaks = NULL, labels = NULL, expand = expansion(mult = c(0, 0))) +
    labs(
      title = paste0("Violin Plot of Variable Space — Cluster ", cnum, " in red"),
      x = NULL, y = "Value"
    ) +
    theme_minimal() +
    theme(
      strip.text  = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x= element_blank()
    )
}
outdistgraph                  <- function(valarray, clusternum, cnum, compartments, popvar) {
  dn <- dimnames(valarray)
  time_labels <- dn$time
  
  to_int_time <- function(x) {
    x2 <- gsub("[^0-9\\-]+", "", x)         # strip non-digits
    suppressWarnings(as.integer(x2))
  }
  
  speccomp     <- length(compartments)/3
  yvals        <- c(popvar$total_gentoo * 1.5, popvar$total_preds * 1.5, popvar$total_albatross * 1.5)
  colorsr      <- c("red", "blue", "orange", "darkgreen")
  plots        <- vector("list", length(compartments))
  names(plots) <- compartments
  i = 1
  for (comp in compartments) {
    # slice: [time x sample] for this compartment
    mat_ts <- drop(valarray[comp, , , drop = FALSE])    # [T x N]
    if (nrow(mat_ts) != length(time_labels)) {
      stop("Row count of mat_ts does not match length of time_labels for compartment ", comp)
    }
    
    # keep only the chosen cluster
    keep <- which(clusternum == cnum)
    if (length(keep) == 0) {
      warning("No samples in cluster ", cnum, " for compartment ", comp, ". Skipping.")
      next
    }
    mat_ts_cluster <- mat_ts[, keep, drop = FALSE]      # [T x K]
    
    # build df with correct axes: rows=time, cols=runs
    df <- as.data.frame(mat_ts_cluster, check.names = FALSE)
    rownames(df) <- time_labels
    
    # long format: time column from rownames; runs come from column names
    df_all <- df %>%
      rownames_to_column(var = "time") %>%
      pivot_longer(cols = -time, names_to = "run", values_to = "value") %>%
      mutate(time = to_int_time(time)) %>%
      filter(is.finite(time), !is.na(value))
    
    # average over runs for each time
    df_avg <- df_all %>%
      group_by(time) %>%
      summarise(avg_value = mean(value, na.rm = TRUE), .groups = "drop")
    
    # plot
    p <- ggplot(df_all, aes(x = time, y = value, group = run)) +
      geom_line(color = "gray60", linetype = "dotted") +
      geom_line(data = df_avg, aes(x = time, y = avg_value),
                inherit.aes = FALSE, color = colorsr[(i %% speccomp) + 1], linewidth = 1.2) +
      scale_y_continuous(
        limits = c(-10, yvals[((i-1) %/% speccomp) + 1])
      ) +
      labs(
        title = paste0(comp, " (Cluster ", cnum, ")"),
        x = "Time (days)",
        y = "Number of Individuals"
      ) +
      theme_minimal()
    
    plots[[i]] <- p 
    i = i + 1
  }
  return(plots)
}

#=======================================================================================================================#
#                                              Validation Functions:
#=======================================================================================================================#
validate_variablecombinations <- function(iterate, cluster, check_iGraph, 
                                          num_samples, pred_Introduce, albs_Introduce, 
                                          midseason_pred, midseason_albs, check_inputDistrib){
  # Combination variable check
  if(iterate == FALSE && cluster == TRUE){
    stop("Error: Cannot do clustering without iteration. Set 'iterate' to TRUE.")
  }
  if(check_iGraph == TRUE && iterate == FALSE){
    stop("Error: Cannot check iGraph without iteration. Set 'iterate' to TRUE.")
  }
  if(check_iGraph == TRUE && num_samples > 30){
    warning("Warning: iGraph check will produce an unreasonable number of graphs.")
  }
  if(cluster == TRUE && num_samples < 1700){
    warning("Warning: Clustering may not be effective with fewer than 1700 samples.")
  }
  if(check_inputDistrib == TRUE && iterate == FALSE){
    stop("Error: Cannot check input distributions without iteration. Set 'iterate' to TRUE.")
  }
  if(check_inputDistrib == TRUE && cluster == FALSE){
    stop("Error: Cannot check input distributions without clustering. Set 'cluster' to TRUE.")
  }
  if(cluster == TRUE && pred_Introduce == FALSE && albs_Introduce == FALSE && 
     midseason_pred == FALSE && midseason_albs == FALSE){
    stop("Error: Cannot do clustering without introducing infection. Set one of the introduction variables to TRUE.")
  }
}
validate_Seasonality          <- function(func, time, timevar){
  # Calculate values
  values <- vapply(time,
                   function(x) func(x, timevar$t_arrival, timevar$mig_day_arrival, timevar$breed_start, 
                                    timevar$hatch_start, timevar$t_departure),
                   FUN.VALUE = numeric(1))
  
  df <- data.frame(
    time = time,
    value = values
  )
  
  # Generate Plot for Values
  p <- ggplot(df, aes(x = time, y = value)) +
    geom_line(color = "steelblue", linewidth = 1) +
    labs(
      title = "Seasonality Transmission Over Time",
      x     = "Time Interval (days)",
      y     = "Computed Value"
    ) +
    theme_minimal()
  
  # Return the values calculated
  list(data = df, plot = p)
}                  
validate_Births               <- function(timevar, popvar, Initial_values, muH_List){
  # Set parameters for no infection
  parmsBirth <- c(
    muH_Gentoo       = muH_List$muH_Gentoo,
    muH_Preds        = muH_List$muH_Preds,
    muH_Albs         = muH_List$muH_Albs,
    
    muI_Gentoo       = 0,
    muI_Preds        = 0,
    muI_Albs         = 0,
    
    gentgent_trans   = 0,
    gentpred_trans   = 0,
    predgent_trans   = 0,
    predpred_trans   = 0,
    predalbs_trans   = 0,
    albspred_trans   = 0,
    albsalbs_trans   = 0,
    
    gamma_gentoo     = 0,
    sigma_gentoo     = 0,
    gamma_preds      = 0,
    sigma_preds      = 0,
    gamma_albatross  = 0,
    sigma_albatross  = 0,
    
    seasons_calc     = FALSE,
    t_arrival        = timevar$t_arrival,
    mig_day_arrival  = timevar$mig_day_arrival,
    breed_start      = timevar$breed_start,
    hatch_start      = timevar$hatch_start,
    t_departure      = timevar$t_departure
  )
  
  popvar$prev_infection_preds <-  0
  popvar$prev_infection_albs  <-  0
  popvar$ms_pred_day          <- -1
  popvar$ms_albs_day          <- -1
  
  # Run the model to get values for the whole time range
  output <- run_Model(timevar, popvar, Initial_values, parmsBirth)
  values <- output$values
  
  # Set what columns to graph
  sus_cols <- c("Susceptible_Gentoo",
                "Susceptible_Predators",
                "Susceptible_Albatross")
  
  # Rearrange values data frame and only have the "Susceptible" compartment
  df_long <- as.data.frame(values) %>%
    select(Time, all_of(sus_cols)) %>%
    pivot_longer(-Time, names_to = "compartment", values_to = "value") %>%
    mutate(
      Species = sub("^Susceptible_", "", compartment),
      Species = factor(Species, levels = c("Gentoo", "Predators", "Albatross"))
    )
  
  # Plot df_long
  p <- ggplot(df_long, aes(x = Time, y = value)) +
    geom_line(color = "steelblue") +
    facet_wrap(~ Species, ncol = 1, scales = "free_y") +
    scale_y_continuous(
      limits = c(-1, NA),                 # force min=0, let max be computed per facet
      expand = expansion(mult = c(0, 0.05)) # small headroom above max
    ) +
    labs(
      title = "Migration, Birth, and Death Dynamics without Infection",
      x = "Time (days)",
      y = "Susceptible"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      panel.spacing.y = unit(0.8, "lines")
    )
  
  # Return p for visualization or saving
  return(p)
}
validate_iGraphs              <- function(values, popvar){
  inf_cols <- c("Infected_Gentoo","Infected_Predators","Infected_Albatross")
  
  df_long <- as.data.frame(values) %>%
    select(Time, all_of(inf_cols)) %>%
    pivot_longer(-Time, names_to = "compartment", values_to = "value") %>%
    mutate(
      Species = sub("^Infected_", "", compartment),
      Species = factor(Species, levels = c("Gentoo","Predators","Albatross"))
    )
  
  caps <- c(
    Gentoo    = as.numeric(popvar$total_gentoo)    * 1.3,
    Predators = as.numeric(popvar$total_preds)     * 1.3,
    Albatross = as.numeric(popvar$total_albatross) * 1.3
  )
  
  ggplot(df_long, aes(x = Time, y = value)) +
    geom_line(color = "orange") +
    facet_wrap(~ Species, ncol = 1, scales = "free_y") +   # <-- make y free
    ggh4x::facetted_pos_scales(
      y = list(
        Species == "Gentoo"    ~ scale_y_continuous(limits = c(0, caps[["Gentoo"]]),    expand = expansion(mult = c(0, 0.05))),
        Species == "Predators" ~ scale_y_continuous(limits = c(0, caps[["Predators"]]), expand = expansion(mult = c(0, 0.05))),
        Species == "Albatross" ~ scale_y_continuous(limits = c(0, caps[["Albatross"]]), expand = expansion(mult = c(0, 0.05)))
      )
    ) +
    labs(title = "Infected Dynamics Over Time", x = "Time (days)", y = "Infected") +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"),
          panel.spacing.y = unit(0.8, "lines"))
}
validate_clustering           <- function(scan, km_final, chosen_k, features_scaled = NULL, 
                                          file = "k_scan_log.txt"){
  # ---- basic validations ----
  stopifnot(is.list(scan),
            !is.null(scan$metrics),
            !is.null(scan$recommended_k),
            is.list(km_final),
            is.numeric(chosen_k), length(chosen_k) == 1L, chosen_k >= 1L)
  
  if (!is.null(km_final$centers)) {
    stopifnot(nrow(km_final$centers) == as.integer(chosen_k))
  }
  
  # prefer provided features; else use what's inside scan (if present)
  if (is.null(features_scaled) && !is.null(scan$scaled_features)) {
    features_scaled <- scan$scaled_features
  }
  
  N <- if (!is.null(features_scaled)) nrow(features_scaled) else length(km_final$cluster)
  p <- if (!is.null(features_scaled)) ncol(features_scaled) else if (!is.null(km_final$centers)) ncol(km_final$centers) else NA_integer_
  
  # extract metrics at chosen_k (if present)
  met <- scan$metrics
  has_k <- "k" %in% names(met)
  row_k <- if (has_k) which(met$k == as.integer(chosen_k)) else integer(0)
  
  metric_at_k <- if (length(row_k) == 1L) {
    as.list(met[row_k, , drop = FALSE])
  } else list()
  
  # cluster sizes from the provided fit
  cl_sizes <- sort(table(km_final$cluster))
  cluster_df <- data.frame(
    cluster = as.integer(names(cl_sizes)),
    size    = as.integer(cl_sizes)
  )
  
  # open connection and write
  con <- file(file, open = "wt"); on.exit(close(con), add = TRUE)
  w <- function(...) writeLines(paste0(...), con)
  
  w("==== K-SCAN AND FINAL CLUSTERING LOG (from existing fit) ====")
  w(sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
  if (!is.null(N)) w(sprintf("N (rows): %s", N))
  if (!is.null(p) && !is.na(p)) w(sprintf("p (features): %s", p))
  w("")
  
  # Recommendations
  w("== Recommended k values (by criterion) ==")
  for (nm in names(scan$recommended_k)) {
    w(sprintf("  %-12s : %s", nm, as.character(scan$recommended_k[[nm]])))
  }
  w("")
  
  # Chosen k and summary from the provided fit
  w("== Final choice and provided kmeans fit ==")
  w(sprintf("  chosen_k        : %d", chosen_k))
  # km_final diagnostics (present in base kmeans objects)
  if (!is.null(km_final$tot.withinss)) w(sprintf("  tot.withinss    : %.6f", km_final$tot.withinss))
  if (!is.null(km_final$betweenss))    w(sprintf("  betweenss       : %.6f", km_final$betweenss))
  if (!is.null(km_final$totss))        w(sprintf("  totss           : %.6f", km_final$totss))
  if (!is.null(km_final$iter))         w(sprintf("  iterations used : %s", km_final$iter))
  w("")
  
  # ---- Nicely aligned metrics table ----
  w("== Metrics across k ==")
  tmp <- scan$metrics
  
  # Format all columns to the same width
  col_widths <- vapply(tmp, function(col) {
    max(nchar(c(names(tmp)[1], as.character(col))), na.rm = TRUE)
  }, numeric(1))
  
  # Create a formatted text matrix
  formatted_lines <- apply(tmp, 1, function(row) {
    paste(
      vapply(seq_along(row), function(i)
        format(row[i], width = col_widths[i], justify = "right"), ""),
      collapse = "  "
    )
  })
  
  # Write header and rows
  header <- paste(
    vapply(names(tmp), function(nm)
      format(nm, width = col_widths[nm], justify = "right"), ""),
    collapse = "  "
  )
  
  writeLines(header, con)
  writeLines(formatted_lines, con)
  w("")
  
  # Metrics at chosen_k (if available)
  if (length(metric_at_k)) {
    w("== Metrics at chosen_k ==")
    for (nm in names(metric_at_k)) {
      val <- metric_at_k[[nm]]
      if (is.numeric(val)) {
        w(sprintf("  %-12s : %s", nm, sprintf("%.6f", val)))
      } else {
        w(sprintf("  %-12s : %s", nm, as.character(val)))
      }
    }
    w("")
  }
  
  # Cluster sizes from the provided fit
  w("== Final cluster size distribution (chosen_k) ==")
  utils::write.table(cluster_df, con, sep = "\t", row.names = FALSE, quote = FALSE)
  w("")
  
  invisible(list(file = normalizePath(file), cluster_sizes = cluster_df))
}

#=======================================================================================================================#
#                                                     Main:
#=======================================================================================================================#
# Validate variable combination used
validate_variablecombinations(iterate, cluster, check_iGraph, 
                              num_samples, pred_Introduce, albs_Introduce, 
                              midseason_pred, midseason_albs, check_inputDistrib)

# Generate popvar and timevar
popvar          <- generate_popvar()
timevar         <- generate_timevar()


# Generate Latin Hypercube parameters
U <- randomLHS(n = num_samples, k = 17)

# Set virus introduction conditions
ifelse(pred_Introduce == TRUE, 
       prev_infections_preds <-       qunif(U[, 1], min = 0.01, max = 0.5),
       prev_infections_preds <- rep(0, num_samples))
ifelse(albs_Introduce == TRUE, 
       prev_infections_albs  <-       qunif(U[, 2], min = 0.01, max = 0.5),
       prev_infections_albs  <- rep(0, num_samples))
ifelse(midseason_pred == TRUE, 
       midseason_pred_days   <- round(qunif(U[, 3], min = timevar$t_arrival, max = timevar$t_departure)), 
       midseason_pred_days   <- rep(-1, num_samples))
ifelse(midseason_albs == TRUE, 
       midseason_albs_days   <- round(qunif(U[, 4], min = timevar$t_arrival, max = timevar$t_departure)), 
       midseason_albs_days   <- rep(-1, num_samples))

# Generate global variables
Initial_values  <- generate_Intitial_values(popvar = popvar)
muH_List        <- generate_muH()
gammasigma_List <- generate_gamma_sigma(U)
gentgents_trans <- qunif(U[, 11], min = 0, max = 0.00005)
gentpreds_trans <- qunif(U[, 12], min = 0, max = 0.00005)
predgents_trans <- qunif(U[, 13], min = 0, max = 0.00005)
predpreds_trans <- qunif(U[, 14], min = 0, max = 0.00005)
predalbss_trans <- qunif(U[, 15], min = 0, max = 0.00005)
albspreds_trans <- qunif(U[, 16], min = 0, max = 0.00005)
albsalbss_trans <- qunif(U[, 17], min = 0, max = 0.00005)
time            <- seq(from = 0, to =  timevar$horizon)

if(iterate == TRUE){
  dfMaster         <- tibble()
  all_series_long  <- list()
  
  # Establish the valarray variable that will hold all the values for the compartments 
  #       for all the different iterations of the model
  compartments <- c(
    "Susceptible_Gentoo",    "Infected_Gentoo",    "Recovered_Gentoo",    "HPAI_Killed_Gentoo",
    "Susceptible_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
    "Susceptible_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross"
  )
  valarray     <- array(
    NA_real_,
    dim      = c(length(compartments), length(0:timevar$horizon), num_samples),
    dimnames = list(compartment = compartments,
                    time        = as.character(0:timevar$horizon),
                    sample      = paste0("sample_", seq_len(num_samples)))
  )
  
  # Create a progress bar
  pb <- progress_bar$new(
    format = "Sample Generation [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
    total = num_samples, clear = FALSE, width = 90
  )
  
  # Iterate over the number of samples
  for(i in 1:num_samples){
    # Set parameters for the current iteration
    parameters <- c(
      t_arrival        = timevar$t_arrival, 
      mig_day_arrival  = timevar$mig_day_arrival,
      breed_start      = timevar$breed_start, 
      hatch_start      = timevar$hatch_start,
      t_departure      = timevar$t_departure,
      
      muH_Gentoo       = muH_List$muH_Gentoo,
      muH_Preds        = muH_List$muH_Preds,
      muH_Albs         = muH_List$muH_Albs,
      
      muI_Gentoo       = gammasigma_List$muIs_gentoo[i],
      muI_Preds        = gammasigma_List$muIs_Preds[i],
      muI_Albs         = gammasigma_List$muIs_Albs[i],
      
      gamma_gentoo     = gammasigma_List$gammas_gentoo[i],
      sigma_gentoo     = gammasigma_List$sigma_gentoo,
      gamma_preds      = gammasigma_List$gammas_preds[i],
      sigma_preds      = gammasigma_List$sigma_preds,
      gamma_albatross  = gammasigma_List$gammas_albatross[i],
      sigma_albatross  = gammasigma_List$sigma_albatross,
      
      gentgent_trans   = gentgents_trans[i],
      gentpred_trans   = gentpreds_trans[i],
      predgent_trans   = predgents_trans[i],
      predpred_trans   = predpreds_trans[i],
      predalbs_trans   = predalbss_trans[i],
      albspred_trans   = albspreds_trans[i],
      albsalbs_trans   = albsalbss_trans[i],
      
      seasons_calc     = seasons_calc,
      t_departure      = timevar$t_departure
    )
    
    # Update the popvar with the current iteration's virus introduction conditions
    popvar$prev_infection_preds <- prev_infections_preds[i]
    popvar$prev_infection_albs  <- prev_infections_albs[i]
    popvar$ms_pred_day          <- midseason_pred_days[i]
    popvar$ms_albs_day          <- midseason_albs_days[i]
    
    # Run the model for the current iteration
    output <- run_Model(timevar, popvar, Initial_values, parameters)
    
    # Save values in the appropriate places
    run_long             <- output$run_long
    run_long$run_id      <- i
    values               <- output$values
    df                   <- output$df
    df$run_id            <- i
    dfMaster             <- bind_rows(dfMaster, df)
    all_series_long[[i]] <- run_long
    
    # Put the values into the valarray
    for(k in 1:length(compartments)){
      valarray[k, , i] <- values[[compartments[k]]]
    }
    # If checking iGraph, generate the plots
    if(check_iGraph == TRUE){
      p <- validate_iGraphs(values, popvar)
      print(p)
    }
    
    pb$tick()
  }
  
  
  # Measure time of execution
  iteration_time <- Sys.time()
  execution_time <- as.numeric(difftime(iteration_time, start_time, units = "secs"))
  hs <- as.integer(execution_time %/% 3600)
  ms <- as.integer((execution_time %% 3600) %/% 60)
  ss <- execution_time %% 60
  print(cat(sprintf("Iteration Complete \n Time: %02d:%02d:%05.2f\n", hs, ms, ss)))
}
if(cluster == TRUE){
  feat             <- build_histogram_features(valarray, relativize_to_species_totals = TRUE, popvar = popvar)
  print("Histogram Features Built.")
  scan             <- evaluate_kmeans_k(feat$features, k_min = k_min, k_max = k_max, nstart = 50, 
                                        iter.max = 100, seed = 1L)
  print("Scan Completed.")
  chosen_k         <- scan$recommended_k$k_Sil
  print(paste("Chosen K:", chosen_k))
  print(paste("Silhouette:", with(scan$metrics, Silhouette[k == chosen_k])))
  km_final         <- fit_kmeans_at_k(scan$scaled_features, k = chosen_k, nstart = 200, iter.max = 300)
  clusters         <- km_final$cluster
  dfMaster$Cluster <- clusters
}

# Validation & Saving protocols
fileLoc = here("Output")
folderLoc = paste0(fileLoc, "\\", format(Sys.time(), "O_%Y-%m-%d_%H-%M-%S"))
if(check_Initial_Values == TRUE ||
   check_season         == TRUE || 
   check_births         == TRUE ||
   check_inputDistrib   == TRUE ||
   check_clustering     == TRUE ||
   write_output         == TRUE ||
   save_data            == TRUE ||
   save_plots           == TRUE ){
  dir.create(folderLoc)
}
if(check_Initial_Values == TRUE){
  write.csv(Initial_values, paste0(folderLoc, "\\Initial_Values.csv"))
}
if(check_season         == TRUE){
  Seas_Val <- validate_Seasonality(seasonality_transmission_calc, time, timevar)
  ggsave(paste0(folderLoc, "\\seasonality.png"), Seas_Val$p, width = 10, height = 5)
}
if(check_births         == TRUE){
  birthsplot <- validate_Births(timevar, popvar, Initial_values, muH_List)
  ggsave(paste0(folderLoc, "\\birthmigrate.png"), birthsplot, width = 10, height = 10)
}
if(check_inputDistrib   == TRUE){
  for (i in 1:chosen_k){
    p <- vardistplots(dfMaster, i)
    ggsave(paste0(folderLoc, "\\InputDistrib_Cluster", i, ".png"), p, width = 14, height = 8)
  }
}
if(check_clustering     == TRUE){
  validate_clustering(scan, km_final, chosen_k, file = paste0(folderLoc, "\\k_scan_log.txt"))
}
if(write_output         == TRUE){
  dfScalar <- dfMaster[, !vapply(dfMaster, is.list, logical(1))]
  write.csv(dfScalar, paste0(folderLoc, "\\InputOutputTable.csv"))
} 
if(save_data            == TRUE){
  saveRDS(dfMaster,        paste0(folderLoc, "\\InputOutputTable.rds"))
  saveRDS(valarray,        paste0(folderLoc, "\\Valarray.rds"))
  saveRDS(all_series_long, paste0(folderLoc, "\\transevents.rds"))
}
if(save_plots           == TRUE){
  for(i in 1:chosen_k){
    plots <- outdistgraph(valarray, clusters, i, compartments, popvar)
    p     <- plot_grid(plotlist = plots, ncol = length(compartments)/3, align = "hv")
    ggsave(paste0(folderLoc, "\\CompartmentGraphs_Cluster_", i,".png"), p, width = 15, height = 10)
  }
}

#=======================================================================================================================#
#                                                   Clean Up:
#=======================================================================================================================#
# Measure time of execution
end_time <- Sys.time()
execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
hs <- as.integer(execution_time %/% 3600)
ms <- as.integer((execution_time %% 3600) %/% 60)
ss <- execution_time %% 60
print(cat(sprintf("Execution Time: %02d:%02d:%05.2f\n", hs, ms, ss)))

# Clear all environmental variables
rm(list = ls())

#=======================================================================================================================#
#                                                   Sources:
#=======================================================================================================================#
# [1] Robin W. Woods, The Birds of the Falkland Islands: An Annotated Checklist (British Ornithologists’ Club, 2017).
# [2] “Gentoo,” Falkland Conservation, n.d., accessed September 8, 2025, https://falklandsconservation.com/gentoo/#:~:
#         text=Breeding%20cycle,when%20they%20form%20small%20cr%C3%A9ches.
# [3] AnAge: The Animal Ageing and Longevity Database, “AnAge entry for Thalassarche melanophris” (Human Ageing Genomic 
#         Resources, March 7, 2023), https://genomics.senescence.info/species/entry.php?species=Thalassarche_melanophris.
# [4] AnAge: The Animal Ageing and Longevity Database, “AnAge entry for Stercorarius antarcticus” (Human Ageing Genomic 
#         Resources, March 7, 2023), https://genomics.senescence.info/species/entry.php?species=Stercorarius_antarcticus.
# [5] Viviane Hénaux and Michael D. Samuel, “AVIAN INFLUENZA SHEDDING PATTERNS IN WATERFOWL: IMPLICATIONS FOR SURVEILLANCE,
#         ENVIRONMENTAL TRANSMISSION, AND DISEASE SPREAD,” Journal of Wildlife Diseases 47, no. 3 (2011): 566–78, 
#         https://doi.org/10.7589/0090-3558-47.3.566.
#         Note: This is not a species specific source, and is an estimate
# [6] Justin D. Brown et al., “Susceptibility of North American Ducks and Gulls to H5N1 Highly Pathogenic Avian Influenza 
#         Viruses,” Emerging Infectious Diseases 12, no. 11 (2006): 1663–70, https://doi.org/10.3201/eid1211.060652.
#         Note: This study was done on laughing gulls, not kelp gulls, or Skuas, and is meant to provide an OOM estimation.
# [7] Williams, T. D. “Annual Variation in Breeding Biology of Gentoo Penguins, Pygoscelis Papua, at Bird Island, South 
#         Georgia.” Journal of Zoology 222, no. 2 (1990): 247–58. https://doi.org/10.1111/j.1469-7998.1990.tb05675.x.
# [8] Catry, Paulo, Jaume Forcada, and Ana Almeida. “Demographic Parameters of Black-Browed Albatrosses Thalassarche 
#         Melanophris from the Falkland Islands.” Polar Biology 34, no. 8 (2011): 1221–29. 
#         https://doi.org/10.1007/s00300-011-0984-3.
# [9] Tickell, W. L. N., and R. Pinder. “BREEDING BIOLOGY OF THE BLACK‐BROWED ALBATROSS DIOMEDEA MELANOPHRIS AND 
#         GREY‐HEADED ALBATROSS D. CHRYSOSTOMA AT BIRD ISLAND, SOUTH GEORGIA.” Ibis 117, no. 4 (1975): 433–51. 
#         https://doi.org/10.1111/j.1474-919X.1975.tb04237.x.
# [10] Catry, Paulo, Ana Almeida, Miguel Lecoq, José Pedro Granadeiro, and Rafael Matias. “Low Breeding Success and Sharp
#         Population Decline at the Largest Known Falkland Skua Colony.” Polar Biology 34, no. 8 (2011): 1239–41. 
#         https://doi.org/10.1007/s00300-011-0978-1.
