start_time           <- Sys.time()
#=======================================================================================================================#
#                                               Control Variables:
#=======================================================================================================================#
# Validation Variable
val                  <- data.frame(check_initVal = FALSE)
val$check_season     <- TRUE
val$check_eventdata1 <- FALSE
val$check_eventdata2 <- FALSE
val$check_births     <- TRUE

# Save Data variables
write_output         <- FALSE
save_plots           <- FALSE
save_data            <- FALSE

# Experimental Variables
seasons_calc         <- TRUE
num_samples          <- 3
numclusters          <- 2
pred_Introduce       <- FALSE
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

#=======================================================================================================================#
#                                                   Populations:
#=======================================================================================================================#
popvar                       <- data.frame(total_gentoo = 2000)                         # Total Gentoo population [Source: Dr. Gamble]
popvar$ResPer_gentoo         <- 1                                                       # Percent of Gentoo that are resident to the island [Source: 1]
popvar$Res_gentoo            <- popvar$total_gentoo * popvar$ResPer_gentoo              # Calculate the number of resident Gentoo
popvar$gentoo_migrate_in     <- popvar$total_gentoo * (1 - popvar$ResPer_gentoo)        # Calculate the number of migratory Gentoo

popvar$total_preds           <- 750                                                     # Total predator population [Source: Dr. Gamble]
popvar$ResPer_preds          <- 0.1                                                     # Percent of predator that are resident to the island [Source: 1]
popvar$Res_preds             <- popvar$total_preds * popvar$ResPer_preds                # Calculate the number of resident predators 
popvar$preds_migrate_in      <- popvar$total_preds * (1 - popvar$ResPer_preds)          # Calculate the number of migratory predators

popvar$total_albatross       <- 50000                                                   # Total Albatross population [Source: Dr. Gamble]
popvar$ResPer_albatross      <- 0                                                       # Percent of Albatross that are resident to the island [Source: 1]
popvar$Res_albatross         <- popvar$total_albatross * popvar$ResPer_albatross        # Resident Albatross 
popvar$albatross_migrate_in  <- popvar$total_albatross * (1 - popvar$ResPer_albatross)  # 100% of albatross are migratory

Initial_values <- c(
  S_Gentoo     = popvar$Res_gentoo ,     # Susceptible
  I_Gentoo     = 0,                      # Infected                   
  R_Gentoo     = 0,                      # Recovered
  D_Gentoo     = 0,                      # Dead
  
  S_Preds      = popvar$Res_preds,       # Susceptible
  I_Preds      = 0,                      # Infected                   
  R_Preds      = 0,                      # Recovered 
  D_Preds      = 0,                      # Dead
  
  S_Albatross  = popvar$Res_albatross,   # Susceptible
  I_Albatross  = 0,                      # Infected                   
  R_Albatross  = 0,                      # Recovered 
  D_Albatross  = 0,                      # Dead
  
  timer        = 0)

if(val$check_initVal == TRUE){
  print(Initial_values)
} 

U <- randomLHS(n = num_samples, k = 17)
ifelse(pred_Introduce == TRUE, prev_infections_preds  <- qunif(U[, 1], min = 0.01, max = 0.5), prev_infections_preds <- rep(0, num_samples))
ifelse(albs_Introduce == TRUE, prev_infections_albs   <- qunif(U[, 2], min = 0.01, max = 0.5), prev_infections_albs  <- rep(0, num_samples))
#=======================================================================================================================#
#                                                 Natural Mortality Rates:
#=======================================================================================================================#
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

#=======================================================================================================================#
#                                             Immunity & Infection Duration:
#=======================================================================================================================#
gentoo_inf_dur    <- qunif(U[, 3], min = 1, max = 15)
gentoo_inf_mort   <- qunif(U[, 4], min = 0.5, max = 0.9)
muIs_gentoo       <- gentoo_inf_mort / gentoo_inf_dur
gammas_gentoo     <- 1 / gentoo_inf_dur
gentoo_imm_dur    <- 182.5                                      # Estimate without a basis
#sigma_gentoo      <- 1 / gentoo_imm_dur                         # Rate of Gentoo immunity loss, >0 makes it an SIRS
sigma_gentoo      <- 0                                          # Rate of Gentoo immunity loss, >0 makes it an SIRS

preds_inf_dur     <- qunif(U[, 5], min = 1, max = 15)
preds_inf_mort    <- qunif(U[, 6], min = 0, max = 0.25)
muIs_Preds        <- preds_inf_mort / preds_inf_dur
gammas_preds      <- 1 / preds_inf_dur
preds_imm_dur     <- 182.5                                      # Estimate without a basis
#sigma_preds       <- 1 / preds_imm_dur                          # Rate of predator immunity loss, >0 makes it an SIRS
sigma_preds       <- 0                                          # Rate of predator immunity loss, >0 makes it an SIRS

Albs_inf_dur      <- qunif(U[, 7], min = 1, max = 15)
Albs_inf_mort     <- qunif(U[, 8], min = 0.5, max = 0.9)
muIs_Albs         <- Albs_inf_mort / Albs_inf_dur
gammas_albatross  <- 1 / Albs_inf_dur
albatross_imm_dur <- 182.5                                      # Estimate without a basis
#sigma_albatross   <- 1 / albatross_imm_dur                      # Rate of Albatross immunity loss, >0 makes it an SIRS
sigma_albatross   <- 0                                          # Rate of Albatross immunity loss, >0 makes it an SIRS

#=======================================================================================================================#
#                                                Temporal Values:
#=======================================================================================================================#
timevar                       <- data.frame(season_start = 100)                # First day of the season for seasonality tranmissibility calculation
timevar$season_end            <- 318                                           # Last day of the season for seasonality tranmissibility calculation
timevar$t_arrival             <- 62                                            # Day at which migratory birds begin arriving, unclear where the source is
timevar$mig_day_arrival       <- 31                                            # The duration of arrival time of migratory birds, no clear source
timevar$t_departure           <- 300                                           # Day at which migratory birds begin departing, unclear where the source is
timevar$mig_day_departure     <- 31                                            # The duration of departure time of migratory birds, no clear source
timevar$gentoo_breed_start    <- 31                                            # Gentoo breeding cycle start [Source: Unclear]
timevar$gentoo_breed_finish   <- 212                                           # Gentoo breeding cycle end   [Source: Unclear]
timevar$pred_breed_start      <- 93                                            # Predator breeding cycle start   [Source: Unclear]
timevar$pred_breed_finish     <- 225                                           # Predator breeding cycle finish  [Source: Unclear]
timevar$albatross_breed_start <- 31
timevar$albatross_breed_finish<- 245
timevar$horizon               <- 365                                                           # Number of days to be modeled
time1                         <- seq(from = 0, to = timevar$t_departure, by = 1)               # Set the front part of the time sequence
time2                         <- seq(from = timevar$t_departure, to = timevar$horizon, by = 1) # Set the back part of the time sequence

ifelse(midseason_pred == TRUE, midseason_pred_days   <- round(qunif(U[,  9], min = timevar$t_arrival + 10, max = timevar$t_departure - 10)), midseason_pred_days <- rep(0, num_samples))
ifelse(midseason_albs == TRUE, midseason_albs_days   <- round(qunif(U[, 10], min = timevar$t_arrival + 10, max = timevar$t_departure - 10)), midseason_albs_days <- rep(0, num_samples))
#=======================================================================================================================#
#                                                    Betas:
#=======================================================================================================================#
gentgents_trans <- qunif(U[, 11], min = 0, max = 0.00005)
gentpreds_trans <- qunif(U[, 12], min = 0, max = 0.00005)
predgents_trans <- qunif(U[, 13], min = 0, max = 0.00005)
predpreds_trans <- qunif(U[, 14], min = 0, max = 0.00005)
predalbss_trans <- qunif(U[, 15], min = 0, max = 0.00005)
albspreds_trans <- qunif(U[, 16], min = 0, max = 0.00005)
albsalbss_trans <- qunif(U[, 17], min = 0, max = 0.00005)

#=======================================================================================================================#
#                                                Event Functions:
#=======================================================================================================================#
create_eventdata1 <- function(timevar, popvar){
  # Duration (days) for each species' arrival window
  L <- timevar$mig_day_arrival
  if (L <= 0) stop("timevar$mig_day_arrival must be > 0")
  
  # --- Arrival windows (sequential) ---
  # Albatross arrive first
  t_alb_start <- timevar$t_arrival
  t_alb_end   <- t_alb_start + L - 1
  alb_times   <- seq(from = t_alb_start, to = t_alb_end, length.out = L)
  
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
    var    = rep(c("S_Albatross","I_Albatross"), each = length(alb_times)),
    time   = rep(alb_times, times = 2),
    value  = rep(c(alb_S_arrival, alb_I_arrival), each = length(alb_times)),
    method = "add"
  )
  
  df_pred <- data.frame(
    var    = rep(c("S_Preds","I_Preds"), each = length(pred_times)),
    time   = rep(pred_times, times = 2),
    value  = rep(c(pred_S_arrival, pred_I_arrival), each = length(pred_times)),
    method = "add"
  )
  
  eventdata1 <- rbind(df_alb, df_pred)
  eventdata1$time <- round(eventdata1$time, 3)
  eventdata1      <- eventdata1[order(eventdata1$time), ]
  
  return(eventdata1)
}
create_eventdata2  <- function(timevar, state){
  mig_departure_times   <- seq(
    from       = timevar$t_departure,
    to         = timevar$t_departure + timevar$mig_day_departure - 1,
    length.out = timevar$mig_day_departure)

  # Calculate the number of birds that depart on each day
  pred_S_departure <- - (0.9 / timevar$mig_day_departure) * state$S_Preds
  pred_I_departure <- - (0.9 / timevar$mig_day_departure) * state$I_Preds
  pred_R_departure <- - (0.9 / timevar$mig_day_departure) * state$R_Preds
  albs_S_departure <- - (1.0 / timevar$mig_day_departure) * state$S_Albatross
  albs_I_departure <- - (1.0 / timevar$mig_day_departure) * state$I_Albatross
  albs_R_departure <- - (1.0 / timevar$mig_day_departure) * state$R_Albatross
  
  # Build eventdata2
  eventdata2 <- data.frame(
    var    = rep(c("S_Preds","I_Preds","R_Preds", "S_Albatross","I_Albatross", "R_Albatross"),
                 each = length(mig_departure_times)),
    time   = rep(mig_departure_times, times = 6),
    value  = rep(c(pred_S_departure, pred_I_departure, pred_R_departure, albs_S_departure, albs_I_departure, albs_R_departure),
                 each = length(mig_departure_times)),
    method = rep("add", 6 * length(mig_departure_times))
  )
  
  # Sort and align precision
  eventdata2      <- eventdata2[order(eventdata2$time), ]
  eventdata2$time <- round(eventdata2$time, 3)
  
  return(eventdata2)
  
}

#=======================================================================================================================#
#                                                Model Functions:
#=======================================================================================================================#
seasonality_transmission_calc <- function(time, season_start, season_end){
  # Parameters for the first sigmoid curve
  k1   = 0.1
  x01  = 31
  a1   = 2
  b1   = 1
  
  # Parameters for the second sigmoid curve
  k2   = 0.1
  x02  = 193
  a2   = 2
  b2   = 1
  
  k3   = 0.1
  x03  = 326.5  # Midpoint of the third sigmoid
  a3   = 1.05    # Upper bound of the third sigmoid (previously `top`)
  b3   = 0.98      # Lower bound of the third sigmoid (previously `bottom`)
  
  ifelse(time < season_start, 
         b1 + (a1 - b1) / (1 + exp(-k1 * (time - x01))),
         ifelse(time <= season_end, 
                a2 - (a2 - b2) * (1 / (1 + exp(-k2 * (time - x02)))),
                b3 + (a3 - b3) / (1 + exp(-k3 * (time - x03)))
         )
  )
}                        # This function calculates the effect seasonality has on transmission rate
make_beta                     <- function(gentgent_trans, gentpred_trans, predgent_trans, predpred_trans, predalbs_trans, albspred_trans, albsalbs_trans, species = c("gentoo", "predators", "albatross")){
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
build_ngm_density             <- function(beta, S0, gamma, muI){
  removal = gamma + muI
  if(any(removal <= 0)) stop("All removal (gamma + muI) must be > 0")
  
  K <- sweep(beta, 2, removal, "/")
  K <- sweep(K, 1, S0, "*")
  K
}
R0_from_ngm                   <- function(K){
  ev <- eigen(K, only.values = TRUE)$values
  list(R0 = max(Mod(ev)), eigenvalues = ev)
}
HPAI_dyn                      <- function(time, state, parameters){
  # Stops analysis if time is not a numeric variable
  if(!is.numeric(time)){
    stop("Time must be a numeric variable")
  }
  
  with(as.list(c(state,parameters)),{
    # If seasonality is activated, use the seasonal infectivity modifier 
    if(seasons_calc == TRUE){
      seasonality_Transmission = seasonality_transmission_calc(time, timevar$season_start, timevar$season_end)
    } else{
      seasonality_Transmission = rep(1, length(time))
    }

    #===============#
    #    Gentoo     #
    #===============#
    dS_Gentoo = 
      ifelse(time >= timevar$gentoo_breed_start & time <= timevar$gentoo_breed_finish, 7 * muH_Gentoo * (S_Gentoo + I_Gentoo + R_Gentoo), 0) - #Gentoo birth rate
      
      seasonality_Transmission * gentgent_trans * I_Gentoo    * S_Gentoo - # Transmission from Gentoo to other Gentoo
      seasonality_Transmission * predgent_trans * I_Preds     * S_Gentoo - # Transmission from Predators to Gentoo to Gentoo

       muH_Gentoo * S_Gentoo  +                                            # Natural death rate
       sigma_gentoo * R_Gentoo                                             # Loss of Gentoo immunity

    dI_Gentoo = 
      seasonality_Transmission * gentgent_trans * I_Gentoo    * S_Gentoo + # Transmission from Gentoo to other Gentoo
      seasonality_Transmission * predgent_trans * I_Preds     * S_Gentoo - # Transmission from Predators to Gentoo
      
      muI_Gentoo * I_Gentoo -                                             # Infected death rate
      gamma_gentoo * I_Gentoo                                             # Recovery rate
    
    dR_Gentoo = 
      gamma_gentoo * I_Gentoo -                                            # Recovery rate
      muH_Gentoo * R_Gentoo -                                              # Natural death rate
      sigma_gentoo * R_Gentoo                                              # Loss of immunity rate

    dD_Gentoo = muI_Gentoo * I_Gentoo
    
    #===============#
    #   Predators   #
    #===============#
    dS_Preds = 
      ifelse(time >= timevar$pred_breed_start & time <= timevar$pred_breed_finish, 25 * muH_Preds * (S_Preds + I_Preds + R_Preds), 0) -  # Predator birth rate
      
      seasonality_Transmission * gentpred_trans * I_Gentoo    * S_Preds -  # Transmission from Gentoo to Predators
      seasonality_Transmission * predpred_trans * I_Preds     * S_Preds -  # Transmission from Predators to other Predators
      seasonality_Transmission * albspred_trans * I_Albatross * S_Preds -  # Transmission from Albatross to Predators
      
      muH_Preds * S_Preds +                                                # Natural death rate
      sigma_preds * R_Preds
    dI_Preds = 
      ifelse(time == midseason_pred_day, 0.01 * S_Preds, 0) +
      seasonality_Transmission * gentpred_trans * I_Gentoo    * S_Preds + # Transmission from Gentoo to Predators
      seasonality_Transmission * predpred_trans * I_Preds     * S_Preds + # Transmission from Predators to other Predators
      seasonality_Transmission * albspred_trans * I_Albatross * S_Preds - # Transmission from Albatross to Predators
      
      muI_Preds * I_Preds -                                               # Infected death rate
      gamma_preds * I_Preds                                               # Predator recovery rate
    dR_Preds =
      gamma_preds * I_Preds -                                             # Predator recovery rate
      muH_Preds * R_Preds -                                               # Natural death rate
      sigma_preds * R_Preds                                               # Loss of immunity
    
    dD_Preds = muI_Preds * I_Preds
    
    #===============#
    #   Albatross   #
    #===============#
    dS_Albatross = 
      ifelse(time >= timevar$albatross_breed_start & time <= timevar$albatross_breed_finish, 17 * muH_Albs * (S_Albatross + I_Albatross + R_Albatross), 0) -
      
      seasonality_Transmission * predalbs_trans * I_Preds     * S_Albatross - # Transmission from Predators to Albatross
      seasonality_Transmission * albsalbs_trans * I_Albatross * S_Albatross - # Transmission from Albatross to other Albatross
      
      muH_Albs * S_Albatross +                                                # Natural death rate
      sigma_albatross * R_Albatross                                           # Loss of Immunity
    
    dI_Albatross = 
      ifelse(time == midseason_albs_day, 0.01 * S_Albatross, 0) +
      seasonality_Transmission * predalbs_trans * I_Preds     * S_Albatross + # Transmission from Predators to Albatross
      seasonality_Transmission * albsalbs_trans * I_Albatross * S_Albatross - # Transmission from Albatross to other Albatross
      
      muI_Albs * I_Albatross -                                                # Infected death rate
      gamma_albatross * I_Albatross                                           # Recovery rate
    
    dR_Albatross = 
      gamma_albatross * I_Albatross -                                         # Recovery rate
      sigma_albatross * R_Albatross -                                         # Loss of immunity
      muH_Albs * R_Albatross                                                  # Natural death rate
    
    dD_Albatross = muI_Albs * I_Albatross
    
    dtimer = 1
    
    return(list(c(
      dS_Gentoo, dI_Gentoo, dR_Gentoo,dD_Gentoo, 
      dS_Preds, dI_Preds, dR_Preds, dD_Preds,
      dS_Albatross, dI_Albatross, dR_Albatross, dD_Albatross,
      
      dtimer)))
  })
}              # This function is the HPAI equations
run_Model                     <- function(timevar, popvar, Initial_values, HPAIdyn, parameters, time1, time2, R0){
  # Create the migration in event
  eventdata1 <- create_eventdata1(timevar, popvar)
  
  # Model Running before migration out
  out_pre            <- deSolve::ode(
    y        = Initial_values,
    func     = HPAI_dyn,
    parms    = parameters,
    times    = time1,
    events   = list(data = eventdata1),
    method   = "lsoda",
  )
  
  # Convert to data frame
  out_pre             <- as.data.frame(out_pre)
  
  # Extract last row of out_pre to get the state at day t_departure, and crease migration out event
  eventdata2 <- create_eventdata2(timevar, out_pre[ nrow(out_pre), ])
  
  # Get the vector of the state at t_departure, without the Time variable and name it so that it matches Intital_variables
  statedeparturevector<- as.numeric(out_pre[ nrow(out_pre), -1])
  names(statedeparturevector) <- names(Initial_values)
  
  # Model running during and after migration out
  out_post           <- deSolve::ode(
    y        = statedeparturevector,
    func     = HPAI_dyn,
    parms    = parameters,
    times    = time2,
    events   = list(data = eventdata2),
    method   = "lsoda",
  )
  
  # Connect pre- and post- migration outputs and then calculate cross-species transmission events
  output <- rbind(out_pre, out_post)
  inc_gentgent <- gentgent_trans * output$S_Gentoo    * output$I_Gentoo
  inc_gentpred <- gentpred_trans * output$S_Preds     * output$I_Gentoo
  inc_predgent <- predgent_trans * output$S_Gentoo    * output$I_Preds
  inc_predpred <- predpred_trans * output$S_Preds     * output$I_Preds
  inc_predalbs <- predalbs_trans * output$S_Albatross * output$I_Preds
  inc_albspred <- albspred_trans * output$S_Preds     * output$I_Albatross
  inc_albsalbs <- albsalbs_trans * output$S_Albatross * output$I_Albatross
  
  # Calculate the cumulative infections for each transmission type, compartment, and in total
  dt     <- diff(output$time)[1]
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
  
  # Label output with human-legible column names
  colnames(output) <- c("Time", 
                        "Susceptable_Gentoo", "Infected_Gentoo", "Recovered_Gentoo", "HPAI_Killed_Gentoo",
                        "Susceptable_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
                        "Susceptable_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross",
                        "Timer")
  output <- output[!duplicated(output$Time), ]
  # Create the return variable with all the data of interest
  #     Note: If the number of variables being returned changes, ensure that the dfMaster function NCol argument
  #           is adjusted accordingly, if not, the code will return an error. 
  df                     <- data.frame(prev_infection_preds = popvar$prev_infection_preds)
  df$midseason_pred_day  <- midseason_pred_day
  df$prev_infection_albs <- popvar$prev_infection_albs
  df$midseason_albs_day  <- midseason_albs_day
  df$gamma_gentoo        <- gamma_gentoo
  df$muI_gentoo          <- muI_Gentoo
  df$gamma_predators     <- gamma_preds
  df$muI_predators       <- muI_Preds
  df$gamma_albatross     <- gamma_albatross
  df$muI_albatross       <- muI_Albs
  df$Beta_GentGent       <- gentgent_trans
  df$Beta_GentPred       <- gentpred_trans
  df$Beta_PredGent       <- predgent_trans
  df$Beta_PredPred       <- predpred_trans
  df$Beta_PredAlbs       <- predalbs_trans
  df$Beta_AlbsPred       <- albspred_trans
  df$Beta_AlbsAlbs       <- albsalbs_trans
  df$R0                  <- R0
  df$maxI_Gentoo         <- max(output["Infected_Gentoo"])
  df$maxI_Predators      <- max(output["Infected_Predators"])
  df$maxI_Albatross      <- max(output["Infected_Albatross"])
  df$terI_Gentoo         <- output[ nrow(output), "Infected_Gentoo"]
  df$terI_Predators      <- output[ nrow(output), "Infected_Predators"]
  df$overwintering       <- ifelse(df$terI_Gentoo >= 1 || df$terI_Predators >= 1, 1, 0)
  df$outbreakdur         <- length(which(output$Infected_Gentoo + output$Infected_Predators + output$Infected_Albatross > 2))
  df$outbreakspan        <- ifelse(df$outbreakdur > (timevar$t_departure - timevar$t_arrival - 5), 1, 0)
  df$infperiod_Gentoo    <- length(which(output$Infected_Gentoo    > 1))
  df$infperiod_Predators <- length(which(output$Infected_Predators > 1))
  df$infperiod_Albatross <- length(which(output$Infected_Albatross > 1))
  df$deadGentoo          <- output[ nrow(output), "HPAI_Killed_Gentoo"]
  df$deadPredators       <- output[ nrow(output), "HPAI_Killed_Predators"]
  df$deadAlbatross       <- output[ nrow(output), "HPAI_Killed_Albatross"]
  df$gentpred_transevents<- cum_gp
  df$predgent_transevents<- cum_pg
  df$predalbs_transevents<- cum_pa
  df$albspred_transevents<- cum_ap
  df$cuminf_Total        <- CI_Total
  df$cuminf_Gentoo       <- CI_Gentoo
  df$cuminf_Predators    <- CI_Predators
  df$cuminf_Albatross    <- CI_Albatross
  #print(df)
  
  return(list(df = df, values = output))
  
} # This function runs a single iteration of the model with a given set of conditions
clusterSamples                <- function(valarray, numclusters, linkage = "average"){
  standardize <- "column"
  
  compartments_use <- c(
    "Infected_Gentoo", "Infected_Predators", "Infected_Albatross",
    "HPAI_Killed_Gentoo", "HPAI_Killed_Predators", "HPAI_Killed_Albatross"
  )
  
  C <- length(compartments_use)
  T <- dim(valarray)[2]
  N <- dim(valarray)[3]
  
  dn <- dimnames(valarray)
  if (!all(compartments_use %in% dn$compartment)) {
    missing <- setdiff(compartments_use, dn$compartment)
    stop("These compartments are missing from valarray: ",
         paste(missing, collapse = ", "))
  }
  
  X <- matrix(NA_real_, nrow = N, ncol = C * T)
  for (j in seq_along(compartments_use)) {
    comp <- compartments_use[j]
    mat  <- valarray[comp, , ]         # [time × sample]
    X[, ((j - 1) * T + 1):(j * T)] <- t(mat)  # transpose to [sample × time]
  }
  
  
  # Drop columns with zero variance (constant across samples)
  keep_cols <- which(apply(X, 2, function(v) sd(v, na.rm = TRUE) > 0))
  X <- X[, keep_cols, drop = FALSE]
  X <- scale(X)
  
  d  <- dist(X, method = "euclidean")
  hc <- hclust(d, method = linkage)
  cl <- cutree(hc, k = numclusters)
  return(cl)
  }

#=======================================================================================================================#
#                                              Validation Functions:
#=======================================================================================================================#
validate_Seasonality          <- function(func, time, timevar){
  # Calculate values
  values <- vapply(time,
                   function(x) func(x, as.numeric(timevar$season_start), as.numeric(timevar$season_end)),
                   FUN.VALUE = numeric(1))
  
  df <- data.frame(
    time = time,
    value = values
  )
  
  # Generate Plot for Values
  p <- ggplot(df, aes(x = time, y = value)) +
    geom_line(color = "steelblue", size = 1) +
    labs(
      title = "Seasonality Transmission Over Time",
      x     = "Time Interval (days)",
      y     = "Computed Value"
    ) +
    theme_minimal()
  
  # Return the values calculated
  list(data = df, plot = p)
}                  # This function validates the "seasonality_transmission_calc"
validate_Births               <- function(HPAI_dyn, time, timevar, popvar, Initial_Values, muH_Gentoo, muH_Preds, muH_Albs, time1, time2){
  muI_Gentoo             <- 0
  sigma_Gentoo           <- 0
  gamma_Gentoo           <- 0
  gentgent_trans         <- 0
  gentpred_trans         <- 0
  
  prev_infection_preds   <- 0
  midseason_pred_day     <- 0
  muI_Preds              <- 0
  sigma_preds            <- 0
  gamma_preds            <- 0
  predgent_trans         <- 0
  predpred_trans         <- 0
  predalbs_trans         <- 0
  
  prev_infection_albs    <- 0
  midseason_albs_day     <- 0
  muI_Albs               <- 0
  sigma_albatross        <- 0
  gamma_albatross        <- 0
  albspred_trans         <- 0
  albsalbs_trans         <- 0
  
  parameters <- c(
    muH_Gentoo,
    muI_Gentoo,
    sigma_Gentoo,
    gamma_Gentoo,
    gentgent_trans,
    gentpred_trans,
    
    prev_infection_preds,
    midseason_pred_day,
    muH_Preds,
    muI_Preds,
    sigma_preds,
    gamma_preds,
    predgent_trans,
    predpred_trans,
    predalbs_trans,
    
    prev_infection_albs,
    midseason_albs_day,
    muH_Albs,
    muI_Albs,
    sigma_albatross,
    gamma_albatross,
    albspred_trans,
    albsalbs_trans
  )
  popvar$prev_infection_preds <- prev_infection_preds
  popvar$prev_infection_albs  <- prev_infection_albs
  
  output <- run_Model(timevar, popvar, Initial_values, HPAIdyn, parameters, time1, time2, 0)
  values <- output$values
  
  print(paste0(output$df))
  
  sus_cols <- c("Susceptable_Gentoo",
                "Susceptable_Predators",
                "Susceptable_Albatross")
  
  df_long <- as.data.frame(values) %>%
    select(Time, all_of(sus_cols)) %>%
    pivot_longer(-Time, names_to = "compartment", values_to = "value") %>%
    mutate(
      Species = sub("^Susceptable_", "", compartment),
      Species = factor(Species, levels = c("Gentoo", "Predators", "Albatross"))
    )
  
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
  
  return(p)
}

#=======================================================================================================================#
#                                                 Generate Graphs:
#=======================================================================================================================#
vardistplots                  <- function(dfMaster, cnum){
  inputdf <- dfMaster %>%
    transmute(
      prev_infection_preds = prev_infection_preds,
      prev_infection_albs  = prev_infection_albs,
      
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
      clusternum = clusternum
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
#                                                    Execution:
#=======================================================================================================================#
# Establish the dfMaster variable that will eventually become the csv save file and then iterate through the model
dfMaster                      <- data.frame(matrix(NA, nrow = num_samples, ncol = 40))

# Establish the valarray variable that will hold all the values for the compartments for all the different iterations of the model
compartments <- c(
  "Susceptable_Gentoo", "Infected_Gentoo", "Recovered_Gentoo", "HPAI_Killed_Gentoo",
  "Susceptable_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
  "Susceptable_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross"
  )
valarray     <- array(
  NA_real_,
  dim      = c(length(compartments), length(0:timevar$horizon), num_samples),
  dimnames = list(compartment = compartments,
                  time        = as.character(0:timevar$horizon),
                  sample      = paste0("sample_", seq_len(num_samples)))
)

print(midseason_albs_days)
print(midseason_pred_days)

# Create a progress bar
pb <- progress_bar$new(
  format = "Sample Generation [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = num_samples, clear = FALSE, width = 90
)

# Run the model for each condition
for(i in 1:num_samples){
  # Set the unique variables for iteration i
  muI_Gentoo           <- muIs_gentoo[i]
  gamma_gentoo         <- gammas_gentoo[i]
  gentgent_trans       <- gentgents_trans[i]
  gentpred_trans       <- gentpreds_trans[i]
  
  prev_infection_preds <- prev_infections_preds[i]
  midseason_pred_day   <- midseason_pred_days[i]
  muI_Preds            <- muIs_Preds[i]
  gamma_preds          <- gammas_preds[i]
  predgent_trans       <- predgents_trans[i]
  predpred_trans       <- predpreds_trans[i]
  predalbs_trans       <- predalbss_trans[i]
  
  prev_infection_albs  <- prev_infections_albs[i]
  midseason_albs_day   <- midseason_albs_days[i]
  muI_Albs             <- muIs_Albs[i]
  gamma_albatross      <- gammas_albatross[i]
  albspred_trans       <- albspreds_trans[i]
  albsalbs_trans       <- albsalbss_trans[i]
  
  if(pred_Introduce == TRUE || albs_Introduce == TRUE){
    # Set the variables necessary for calculating R0
    beta    <- make_beta(gentgent_trans, gentpred_trans, predgent_trans, predpred_trans, predalbs_trans, albspred_trans, albsalbs_trans)
    S0      <- c(gentoo = 0,            predators = popvar$total_preds - (popvar$total_preds * prev_infection_preds), albatross = popvar$total_albatross - (popvar$total_albatross * prev_infection_albs))
    gamma   <- c(gentoo = gamma_gentoo, predators = gamma_preds, albatross = gamma_albatross) 
    muI     <- c(gentoo = muI_Gentoo,   predators = muI_Preds,   albatross = muI_Albs)
    
    # Calculate R0
    K    <- build_ngm_density(beta, S0, gamma, muI)
    spec <- R0_from_ngm(K)
    R0   <- spec$R0
  } else{
    R0   <- NA
  }

  # Set the parameters based on the variables pulled out of the LHC sets
  parameters <- c(
    muH_Gentoo,
    muI_Gentoo,
    sigma_gentoo,
    gamma_gentoo,
    gentgent_trans,
    gentpred_trans,
    
    prev_infection_preds,
    midseason_pred_day,
    muH_Preds,
    muI_Preds,
    sigma_preds,
    gamma_preds,
    predgent_trans,
    predpred_trans,
    predalbs_trans,
    
    prev_infection_albs,
    midseason_albs_day,
    muH_Albs,
    muI_Albs,
    sigma_albatross,
    gamma_albatross,
    albspred_trans,
    albsalbs_trans
  )
  popvar$prev_infection_preds <- prev_infection_preds
  popvar$prev_infection_albs  <- prev_infection_albs
  
  # Run model and save it appropriately
  output                      <- run_Model(timevar, popvar, Initial_values, HPAIdyn, parameters, time1, time2, R0)
  parameters                  <- NULL
  vec_i                       <- unlist(output$df, use.names = FALSE)
  dfMaster[i,]                <- vec_i
  
  for(k in 1:length(compartments)){
    valarray[compartments[k], , i] = output$values[[compartments[k]]]
  }
  
  pb$tick()
}

# Find the clusters, and tie them to the dfMaster dataframe
colnames(dfMaster) <- colnames(output$df)
clusternum <- clusterSamples(valarray, numclusters)
dfMaster <- cbind(dfMaster, clusternum)

# Save the dfMaster dataframe for record keeping
fileLoc = "C:\\Users\\tjm336\\OneDrive - Cornell University\\01 Lab Work\\01 Rotations\\01 Bento-Gamble\\04-TriSpecies-HPAI\\LHC Walking\\Output"
folderLoc = paste0(fileLoc, "\\", format(Sys.time(), "O_%Y-%m-%d_%H-%M-%S"))
if(write_output == TRUE || save_plots == TRUE || val$check_season == TRUE || val$check_births == TRUE){
  dir.create(folderLoc)
}

# Save the output file as a csv
if(write_output == TRUE){
  write.csv(dfMaster, paste0(folderLoc, "\\InputOutputTable.csv"))
} 

# Graph each cluster
if(save_plots == TRUE){
  for(i in 1:numclusters){
    pbw <- vardistplots(dfMaster, i)
    ggsave(paste0(folderLoc, "\\InputDist_Cluster_", i,".png"), pbw, width = 8, height = 10)
    plots <- outdistgraph(valarray, clusternum, i, compartments, popvar)
    p <- plot_grid(plotlist = plots, ncol = length(compartments)/3, align = "hv")
    ggsave(paste0(folderLoc, "\\CompartmentGraphs_Cluster_", i,".png"), p, width = 15, height = 10)
  } 
}

#Save valarray
if(save_data == TRUE){
  saveRDS(valarray, paste0(folderLoc, "\\valarray.rds"))
}

# Make and save seasonality transmission graph
if(val$check_season == TRUE){
  Seas_Val <- validate_Seasonality(seasonality_transmission_calc, seq(1,timevar$horizon), timevar)
  ggsave(paste0(folderLoc, "\\seasonality.png"), Seas_Val$plot)
}

# Make and save the validate births graph
if(val$check_birth == TRUE){
  pbirths <- validate_Births(HPAI_dyn, time, timevar, popvar, Initial_Values, muH_Gentoo, muH_Preds, muH_Albs, time1, time2)
  ggsave(paste0(folderLoc, "\\migrationrate.png"), pbirths, width = 8, height = 8)  
}

# Measure time of execution
end_time <- Sys.time()
execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
hs <- as.integer(execution_time %/% 3600)
ms <- as.integer((execution_time %% 3600) %/% 60)
ss <- execution_time %% 60
print(cat(sprintf("Execution Time: %02d:%02d:%05.2f\n", hs, ms, ss)))

#=======================================================================================================================#
#                                                         Sources:
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
