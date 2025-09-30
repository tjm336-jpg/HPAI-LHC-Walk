start_time           <- Sys.time()
#=======================================================================================================================#
#                                               Control Variables:
#=======================================================================================================================#
# Validation Variable
check_Initial_Values <- FALSE
check_season         <- FALSE
check_births         <- FALSE

# Save Data variables
#write_output         <- FALSE
#save_plots           <- FALSE
#save_data            <- FALSE

# Experimental Variables
#seasons_calc         <- TRUE
iterate              <- TRUE
num_samples          <- 3
#numclusters          <- 2
#doclusters           <- FALSE
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
library(tibble)
library(dplyr)

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
    
    S_Predators = popvar$Res_preds,       # Susceptible
    I_Predators = 0,                      # Infected
    
    S_Albatross = popvar$Res_albatross,   # Susceptible
    I_Albatross = 0,                      # Infected
    
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
  gentoo_inf_dur    <- qunif(U[,  5], min = 1, max = 15)
  gentoo_inf_mort   <- qunif(U[,  6], min = 0.5, max = 0.95)
  muIs_gentoo       <- gentoo_inf_mort / gentoo_inf_dur
  gammas_gentoo     <- 1 / gentoo_inf_dur
  gentoo_imm_dur    <- 182.5                                      # Estimate without a basis
  #sigma_gentoo      <- 1 / gentoo_imm_dur                         # Rate of Gentoo immunity loss, >0 makes it an SIRS
  sigma_gentoo      <- 0                                          # Rate of Gentoo immunity loss, >0 makes it an SIRS
  
  preds_inf_dur     <- qunif(U[,  7], min = 1, max = 15)
  preds_inf_mort    <- qunif(U[,  8], min = 0, max = 0.5)
  muIs_Preds        <- preds_inf_mort / preds_inf_dur
  gammas_preds      <- 1 / preds_inf_dur
  preds_imm_dur     <- 182.5                                      # Estimate without a basis
  #sigma_preds       <- 1 / preds_imm_dur                          # Rate of predator immunity loss, >0 makes it an SIRS
  sigma_preds       <- 0                                          # Rate of predator immunity loss, >0 makes it an SIRS
  
  Albs_inf_dur      <- qunif(U[,  9], min = 1, max = 15)
  Albs_inf_mort     <- qunif(U[, 10], min = 0.5, max = 0.95)
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
  
  migineventdata      <- rbind(df_alb, df_pred)
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
seasonality_transmission_calc <- function(time, timevar){
  if(time < (timevar$t_arrival + (timevar$mig_day_arrival / 2))){
    1 # Before arrival
  } else{
    if(time < timevar$breed_start){
      2 # During Arrival/Mating
    } else{
      if(time < timevar$hatch_start){
        1 # After eggs are laid
      } else{
        if(time < timevar$t_departure){
          2 # After eggs are hatched
        } else{
          1 # After migration out
        }
      }
    }
  }
}
HPAI_dyn                      <- function(time, state, parameters){
  # Extract state
  S_Gentoo    <- round(state["S_Gentoo"]   , 4)
  I_Gentoo    <- round(state["I_Gentoo"]   , 4)
  S_Predators <- round(state["S_Predators"], 4)
  I_Predators <- round(state["I_Predators"], 4)
  S_Albatross <- round(state["S_Albatross"], 4)
  I_Albatross <- round(state["I_Albatross"], 4)
  
  timer       <- state["timer"]
  
  # Extract parameters
  muH_Gentoo     <- parameters["muH_Gentoo"]
  muH_Preds      <- parameters["muH_Preds"]
  muH_Albs       <- parameters["muH_Albs"]
  
  muI_Gentoo     <- parameters["muI_Gentoo"]
  muI_Preds      <- parameters["muI_Preds"]
  muI_Albs       <- parameters["muI_Albs"]
  
  gentgent_trans <- parameters["gentgent_trans"]
  gentpred_trans <- parameters["gentpred_trans"]
  predgent_trans <- parameters["predgent_trans"]
  predpred_trans <- parameters["predpred_trans"]
  predalbs_trans <- parameters["predalbs_trans"]
  albspred_trans <- parameters["albspred_trans"]
  albsalbs_trans <- parameters["albsalbs_trans"]
  
  t_departure    <- parameters["t_departure"]
  
  #===============#
  #    Gentoo     #
  #===============#
  dS_Gentoo = 
    - muH_Gentoo     * S_Gentoo                  -
      gentgent_trans * S_Gentoo    * I_Gentoo    -
      predgent_trans * S_Gentoo    * I_Predators 
  
  dI_Gentoo = 
    - muI_Gentoo     * I_Gentoo                  +
      gentpred_trans * S_Gentoo    * I_Gentoo    +
      predgent_trans * S_Gentoo    * I_Predators 

  #===============#
  #   Predators   #
  #===============#
  dS_Predators = 
    - muH_Preds      * S_Predators               +
      gentpred_trans * S_Predators * I_Gentoo    -
      predpred_trans * S_Predators * I_Predators -
      albspred_trans * S_Predators * I_Albatross
  
  dI_Predators = 
    - muI_Preds      * S_Predators               +
      gentpred_trans * S_Predators * I_Gentoo    +
      predpred_trans * S_Predators * I_Predators +
      albspred_trans * S_Predators * I_Albatross
  
  #===============#
  #   Albatross   #
  #===============#
  dS_Albatross = 
    ifelse(time < t_departure, - muH_Albs * S_Albatross, 0) -
      predalbs_trans * S_Albatross * I_Predators -
      albsalbs_trans * S_Albatross * I_Albatross 
  
  dI_Albatross = 
    ifelse(time < t_departure, - muI_Albs * S_Albatross, 0) +
      predalbs_trans * S_Albatross * I_Predators +
      albsalbs_trans * S_Albatross * I_Albatross 
  
  dTimer = 1
  
  list(c(dS_Gentoo, dI_Gentoo, dS_Predators, dI_Predators, dS_Albatross, dI_Albatross, dTimer))
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
  N_Gent <- as.numeric(statedeparturevector["S_Gentoo"])    + as.numeric(statedeparturevector["I_Gentoo"])
  N_Pred <- as.numeric(statedeparturevector["S_Predators"]) + as.numeric(statedeparturevector["I_Predators"])
  N_Albs <- as.numeric(statedeparturevector["S_Albatross"]) + as.numeric(statedeparturevector["I_Albatross"])
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
  N_Gent <- as.numeric(statedeparturevector["S_Gentoo"])    + as.numeric(statedeparturevector["I_Gentoo"])
  N_Pred <- as.numeric(statedeparturevector["S_Predators"]) + as.numeric(statedeparturevector["I_Predators"])
  N_Albs <- as.numeric(statedeparturevector["S_Albatross"]) + as.numeric(statedeparturevector["I_Albatross"])
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
    "Susceptible_Gentoo", "Infected_Gentoo",
    "Susceptible_Predators", "Infected_Predators",
    "Susceptible_Albatross", "Infected_Albatross",
    "Timer"
  )
  out <- out[!duplicated(out$Time),]
  
  # Return the dataframe
  return(out)
}

#=======================================================================================================================#
#                                              Validation Functions:
#=======================================================================================================================#
validate_Seasonality          <- function(func, time, timevar){
  # Calculate values
  values <- vapply(time,
                   function(x) func(x, timevar),
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
    
    gammas_gentoo    = 0,
    sigma_gentoo     = 0,
    gammas_preds     = 0,
    sigma_preds      = 0,
    gammas_albatross = 0,
    sigma_albatross  = 0,
    
    gentgent_trans   = 0,
    gentpred_trans   = 0,
    predgent_trans   = 0,
    predpred_trans   = 0,
    predalbs_trans   = 0,
    albspred_trans   = 0,
    albsalbs_trans   = 0,
    
    t_departure      = timevar$t_departure
  )
  
  popvar$prev_infection_preds <- 0
  popvar$prev_infection_albs  <- 0
  
#  parameters <- NULL
  # Rune the model to get values for the whole spectrum
  values <- run_Model(timevar, popvar, Initial_values, parmsBirth)
  
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

#=======================================================================================================================#
#                                                     Main:
#=======================================================================================================================#
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
       midseason_pred_days   <- round(qunif(U[, 3], min = timevar$t_arrival + 10, max = timevar$t_departure - 10)), 
       midseason_pred_days   <- rep(-1, num_samples))
ifelse(midseason_albs == TRUE, 
       midseason_albs_days   <- round(qunif(U[, 4], min = timevar$t_arrival + 10, max = timevar$t_departure - 10)), 
       midseason_albs_days   <- rep(-1, num_samples))

# Generate global variables
popvar          <- generate_popvar()
Initial_values  <- generate_Intitial_values(popvar = popvar)
timevar         <- generate_timevar()
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

print(gammasigma_List$sigma_gentoo)

if(iterate == TRUE){
  for(i in 1:num_samples){
    popvar$prev_infection_preds <- prev_infections_preds[i]
    popvar$prev_infection_albs  <- prev_infections_albs[i]
    
    
    parameters <- c(
      muH_Gentoo       = muH_List$muH_Gentoo,
      muH_Preds        = muH_List$muH_Preds,
      muH_Albs         = muH_List$muH_Albs,
      
      muI_Gentoo       = gammasigma_List$muIs_gentoo[i],
      muI_Preds        = gammasigma_List$muIs_Preds[i],
      muI_Albs         = gammasigma_List$muIs_Albs[i],
      
      gentgent_trans   = gentgents_trans[i],
      gentpred_trans   = gentpreds_trans[i],
      predgent_trans   = predgents_trans[i],
      predpred_trans   = predpreds_trans[i],
      predalbs_trans   = predalbss_trans[i],
      albspred_trans   = albspreds_trans[i],
      albsalbs_trans   = albsalbss_trans[i],
      
      t_departure      = timevar$t_departure
    )
    
    values <- run_Model(timevar, popvar, Initial_values, parameters)
    print(values)
  }
}

# Validation protocols
if(check_Initial_Values == TRUE){
  print("Initial Values:")
  print(Initial_values)
}
if(check_season         == TRUE){
  Seas_Val <- validate_Seasonality(seasonality_transmission_calc, time, timevar)
  print(Seas_Val$p)
}
if(check_births         == TRUE){
  birthsplot <- validate_Births(timevar, popvar, Initial_values, muH_List)
  print(birthsplot)
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
#rm(list = ls())

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

