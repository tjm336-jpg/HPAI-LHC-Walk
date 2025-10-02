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

# Save Data variables
write_output         <- TRUE
save_plots           <- TRUE
save_data            <- TRUE

# Experimental Variables
seasons_calc         <- TRUE
iterate              <- TRUE
cluster              <- TRUE
num_samples          <- 1500
num_clusters         <- 30
pred_Introduce       <- TRUE
albs_Introduce       <- FALSE
midseason_pred       <- FALSE
midseason_albs       <- FALSE

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
  
  S_Gentoo_Migout    <- S_Gentoo    * (1 - popvar$ResPer_gentoo)    / L
  S_Predators_Migout <- S_Predators * (1 - popvar$ResPer_preds)     / L
  S_Albatross_Migout <- S_Albatross * (1 - popvar$ResPer_albatross) / L

  I_Gentoo_Migout    <- I_Gentoo    * (1 - popvar$ResPer_gentoo)    / L
  I_Predators_Migout <- I_Predators * (1 - popvar$ResPer_preds)     / L
  I_Albatross_Migout <- I_Albatross * (1 - popvar$ResPer_albatross) / L
  
  migout_eventdata <- data.frame(
    var    = rep(c("S_Gentoo","S_Predators", "S_Albatross", 
                   "I_Gentoo", "I_Predators", "I_Albatross"),  
                 each = length(migout_times)),
    time   = rep(migout_times, times = 6),
    value  = rep(c(- S_Gentoo_Migout, - S_Predators_Migout, - S_Albatross_Migout, 
                   - I_Gentoo_Migout, - I_Predators_Migout, - I_Albatross_Migout), 
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
  
  # First derivative of the time series
  dS_G <- diff(out$Susceptible_Gentoo)
  dI_G <- diff(out$Infected_Gentoo)
  dS_P <- diff(out$Susceptible_Predators)
  dI_P <- diff(out$Infected_Predators)
  dS_A <- diff(out$Susceptible_Albatross)
  dI_A <- diff(out$Infected_Albatross)

  # Find local maxima of the time series
  dS_G_Peaks <- which(diff(sign(dS_G)) == -2) + 1
  dI_G_Peaks <- which(diff(sign(dI_G)) == -2) + 1
  dS_P_Peaks <- which(diff(sign(dS_P)) == -2) + 1
  dI_P_Peaks <- which(diff(sign(dI_P)) == -2) + 1
  dS_A_Peaks <- which(diff(sign(dS_A)) == -2) + 1
  dI_A_Peaks <- which(diff(sign(dI_A)) == -2) + 1

  # Extract the times of the peaks
  tS_G <- out$Time[dS_G_Peaks]; tI_G <- out$Time[dI_G_Peaks]
  tS_P <- out$Time[dS_P_Peaks]; tI_P <- out$Time[dI_P_Peaks]
  tS_A <- out$Time[dS_A_Peaks]; tI_A <- out$Time[dI_A_Peaks]

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
  df$terI_Gentoo         <- out[ nrow(out), "Infected_Gentoo"]
  df$terI_Predators      <- out[ nrow(out), "Infected_Predators"]
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
  df$S_G_PeakNum         <- length(tS_G); df$I_G_PeakNum         <- length(tI_G)
  df$S_P_PeakNum         <- length(tS_P); df$I_P_PeakNum         <- length(tI_P)
  df$S_A_PeakNum         <- length(tS_A); df$I_A_PeakNum         <- length(tI_A)
  df$S_G_Peak            <- list(tS_G);   df$I_G_Peak            <- list(tI_G)
  df$S_P_Peak            <- list(tS_P);   df$I_P_Peak            <- list(tI_P)
  df$S_A_Peak            <- list(tS_A);   df$I_A_Peak            <- list(tI_A)

  
  # Return the dataframe
  return(list(df       = df, 
              values   = out, 
              run_long = run_long))
}
clusterSamples_euclidean      <- function(valarray, numclusters, linkage = "average", popvar,
                                          normalize_within_run = c("none", "by_compartment", "by_species", "global_row"),
                                          scale_across_runs = TRUE){
  # normalize_within_run:
  #  "none"           = original behavior (no within-run normalization)
  #  "by_compartment" = for each run and each compartment's T-time block, z-score within that run
  #  "global_row"     = z-score the entire row vector (all compartments & time) within that run
  normalize_within_run <- match.arg(normalize_within_run)
  
  compartments_use <- c(
    "Infected_Gentoo", "Infected_Predators", "Infected_Albatross",
    "HPAI_Killed_Gentoo", "HPAI_Killed_Predators", "HPAI_Killed_Albatross"
  )
  
  C <- length(compartments_use)
  T <- dim(valarray)[2]
  N <- dim(valarray)[3]
  
  dn <- dimnames(valarray)
  if (is.null(dn) || is.null(dn$compartment)) {
    stop("valarray must have dimnames with a 'compartment' entry.")
  }
  if (!all(compartments_use %in% dn$compartment)) {
    missing <- setdiff(compartments_use, dn$compartment)
    stop("These compartments are missing from valarray: ",
         paste(missing, collapse = ", "))
  }
  
  total_gentoo <- popvar$total_gentoo
  total_preds  <- popvar$total_preds
  total_albs   <- popvar$total_albatross
  
  if (is.null(total_gentoo) || is.null(total_preds) || is.null(total_albs)) {
    stop("popvar must contain total_gentoo, total_preds, and total_albatross.")
  }
  
  # Build design matrix X: rows = runs, columns = concatenated timepoints per compartment
  X <- matrix(NA_real_, nrow = N, ncol = C * T)
  for (j in seq_along(compartments_use)) {
    comp <- compartments_use[j]
    mat  <- valarray[comp, , ]                 # [time × sample]
    X[, ((j - 1) * T + 1):(j * T)] <- t(mat)   # [sample × time]
  }
  
  # ---- Within-run normalization (optional) ----
  if (normalize_within_run == "by_species") {
    # Map each compartment block to its species total
    species_for_comp <- c("Gentoo", "Predators", "Albatross",
                          "Gentoo", "Predators", "Albatross")
    
    totals <- list(
      Gentoo    = total_gentoo,
      Predators = total_preds,
      Albatross = total_albs
    )
    
    for (j in seq_along(compartments_use)) {
      idx <- ((j - 1) * T + 1):(j * T)
      denom <- totals[[ species_for_comp[j] ]]
      # Divide the entire block (all timepoints, all runs) by the species total
      X[, idx] <- X[, idx] / denom
    }
    # At this point, values are per-capita (or per-unit) trajectories by species.
  }
  
  # Helper: safe z-score (returns 0 if sd == 0 or all NA)
  safe_z <- function(v) {
    if (all(is.na(v))) return(v)
    m <- mean(v, na.rm = TRUE)
    s <- stats::sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      # flat block -> set to 0 so it doesn't influence distances
      return(rep(0, length(v)))
    } else {
      return((v - m) / s)
    }
  }
  
  if (normalize_within_run == "by_compartment") {
    # For each run (row) and each compartment block of length T, z-score within that run
    for (i in 1:N) {
      for (j in 1:C) {
        idx <- ((j - 1) * T + 1):(j * T)
        X[i, idx] <- safe_z(X[i, idx])
      }
    }
  } else if (normalize_within_run == "global_row") {
    # Z-score the entire row vector within the run (across all compartments & time)
    for (i in 1:N) {
      X[i, ] <- safe_z(X[i, ])
    }
  }
  # If "none", do nothing here.
  
  # ---- Across-run cleaning & scaling (as in original) ----
  # Drop columns with zero variance across runs
  keep_cols <- which(apply(X, 2, function(v) stats::sd(v, na.rm = TRUE) > 0))
  X <- X[, keep_cols, drop = FALSE]
  
  # Optional across-run standardization (column-wise z-score across runs)
  if (scale_across_runs) {
    # base::scale can propagate NA if an NA remains in a column; you can pre-impute if needed
    X <- scale(X)
  }
  
  # Distances & clustering
  d  <- dist(X, method = "euclidean")
  cat("Distance calculated: \n", summary(d), "\n")
  hc <- hclust(d, method = linkage)
  cl <- cutree(hc, k = numclusters)
  return(cl)
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
  if(num_clusters < 2 && cluster == TRUE){
    stop("Error: Number of clusters must be at least 2 when clustering is enabled.")
  }
  if(num_samples < 3 && cluster == TRUE){
    stop("Error: Number of samples must be at least 3 when clustering is enabled.")
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
  clusters <- clusterSamples_euclidean(valarray, num_clusters, linkage = "average", popvar,
                                       normalize_within_run = "by_species", scale_across_runs = TRUE)
  dfMaster$Cluster <- as.factor(clusters)
}

# Validation & Saving protocols
fileLoc = "C:\\Users\\tjm336\\OneDrive - Cornell University\\01 Lab Work\\01 Rotations\\01 Bento-Gamble\\04-TriSpecies-HPAI\\LHC Walking\\Output"
folderLoc = paste0(fileLoc, "\\", format(Sys.time(), "O_%Y-%m-%d_%H-%M-%S"))
if(check_Initial_Values == TRUE ||
   check_season         == TRUE || 
   check_births         == TRUE ||
   check_inputDistrib   == TRUE ||
   write_output         == TRUE ||
   save_data            == TRUE ||
   save_plots           == TRUE  ){
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
  for (i in 1:num_clusters){
    p <- vardistplots(dfMaster, i)
    ggsave(paste0(folderLoc, "\\InputDistrib_Cluster", i, ".png"), p, width = 14, height = 8)
  }
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
  for(i in 1:num_clusters){
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

