#=======================================================================================================================#
#                                               Compatability Notes:
#=======================================================================================================================#
# Compatible with LHC-Walk-v02.R

#=======================================================================================================================#
#                                               Control Variables:
#=======================================================================================================================#
# Validation Variable
val                  <- data.frame(check_initVal = FALSE)
val$check_eventdata1 <- FALSE
val$check_eventdata2 <- FALSE
start_time           <- Sys.time()

# Experimental Variables
seasons_calc     <- TRUE

#=======================================================================================================================#
#                                                   Libraries:
#=======================================================================================================================#
library(ggplot2)
library(deSolve)
library(tibble)
library(tidyr)
library(dplyr)
library(pheatmap)

#=======================================================================================================================#
#                                                   Populations:
#=======================================================================================================================#
popvar                       <- data.frame(colony_num = 20)                             # Number of colonies on the island

popvar$total_gentoo          <- 200000 / popvar$colony_num                              # Total Gentoo population [Source: 1]
popvar$ResPer_gentoo         <- 1                                                       # Percent of Gentoo that are resident to the island [Source: 1]
popvar$Res_gentoo            <- popvar$total_gentoo * popvar$ResPer_gentoo              # Calculate the number of resident Gentoo
popvar$gentoo_migrate_in     <- popvar$total_gentoo * (1 - popvar$ResPer_gentoo)        # Calculate the number of migratory Gentoo

popvar$total_preds           <- 74000 / popvar$colony_num                               # Total predator population [Source: 1]
popvar$ResPer_preds          <- 0.1                                                     # Percent of predator that are resident to the island [Source: 1]
popvar$Res_preds             <- popvar$total_preds * popvar$ResPer_preds                # Calculate the number of resident predators 
popvar$preds_migrate_in      <- popvar$total_preds * (1 - popvar$ResPer_preds)          # Calculate the number of migratory predators

popvar$total_albatross       <- 951000 / popvar$colony_num                              # Total Albatross population [Source: 1]
popvar$ResPer_albatross      <- 0                                                       # Percent of Albatross that are resident to the island [Source: 1]
popvar$Res_albatross         <- popvar$total_albatross * popvar$ResPer_albatross        # Resident Albatross 
popvar$albatross_migrate_in  <- popvar$total_albatross * (1 - popvar$ResPer_albatross) # 100% of albatross are migratory

Initial_values <- c(
  S_Gentoo = popvar$Res_gentoo ,         # Susceptible
  I_Gentoo = 0,                          # Infected                   
  R_Gentoo = 0,                          # Recovered
  D_Gentoo = 0,                          # Dead
  
  S_Preds      = popvar$Res_preds,       # Susceptible
  I_Preds      = 0,                      # Infected                   
  R_Preds      = 0,                      # Recovered 
  D_Preds      = 0,                     # Dead
  
  S_Albatross  = popvar$Res_albatross,   # Susceptible
  I_Albatross  = 0,                      # Infected                   
  R_Albatross  = 0,                      # Recovered 
  D_Albatross  = 0,                      # Dead
  
  timer        = 0)

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
gentoo_imm_dur    <- 182.5                                      # Estimate without a basis
sigma_gentoo      <- 1 / gentoo_imm_dur                         # Rate of Gentoo immunity loss

Preds_imm_dur     <- 182.5                                      # Estimate without a basis
sigma_preds       <- 1 / Preds_imm_dur                          # Rate of predator immunity loss

albatross_imm_dur <- 182.5                                      # Estimate without a basis
sigma_albatross   <- 1 / albatross_imm_dur                      # Rate of Albatross immunity loss

#=======================================================================================================================#
#                                                Temporal Values:
#=======================================================================================================================#
timevar                       <- data.frame(season_start = 100)                # First day of the season for seasonality tranmissibility calculation
timevar$season_end            <- 349                                           # Last day of the season for seasonality tranmissibility calculation
timevar$t_arrival             <- 62                                            # Day at which migratory birds begin arriving, unclear where the source is
timevar$mig_day_arrival       <- 31                                            # The duration of arrival time of migratory birds, no clear source
timevar$t_departure           <- 300                                           # Day at which migratory birds begin departing, unclear where the source is
timevar$mig_day_departure     <- 10                                            # The duration of departure time of migratory birds, no clear source
timevar$gentoo_breed_start    <- 31                                            # Gentoo breeding cycle start [Source: Unclear]
timevar$gentoo_breed_finish   <- 212                                           # Gentoo breeding cycle end   [Source: Unclear]
timevar$pred_breed_start      <- 93                                            # Predator breeding cycle start   [Source: Unclear]
timevar$pred_breed_finish     <- 225                                           # Predator breeding cycle finish  [Source: Unclear]
timevar$albatross_breed_start <- 31
timevar$albatross_breed_finish<- 245
timevar$horizon               <- 365                                                           # Number of days to be modeled
time1                         <- seq(from = 0, to = timevar$t_departure, by = 1)               # Set the front part of the time sequence
time1                         <- round(time1, 3)                                               # Round all the times
time2                         <- seq(from = timevar$t_departure, to = timevar$horizon, by = 1) # Set the back part of the time sequence

#=======================================================================================================================#
#                                                Event Functions:
#=======================================================================================================================#
create_eventdata1  <- function(timevar, popvar){
  # Create a list of the number of days where these arrival will happen
  mig_arrival_times     <- seq(
    from       = timevar$t_arrival, 
    to         = timevar$t_arrival + timevar$mig_day_arrival -1, 
    length.out = timevar$mig_day_arrival) 
  
  # Calculate the number of birds that arrive on each day
  pred_S_arrival <- (1 - popvar$prev_infection_preds) * popvar$preds_migrate_in / timevar$mig_day_arrival
  pred_I_arrival <-  popvar$prev_infection_preds       * popvar$preds_migrate_in / timevar$mig_day_arrival
  alb_S_arrival  <- (1 - popvar$prev_infection_albs)  * popvar$albatross_migrate_in / timevar$mig_day_arrival
  alb_I_arrival  <-  popvar$prev_infection_albs        * popvar$albatross_migrate_in / timevar$mig_day_arrival
  
  # Build eventdata1
  eventdata1 <- data.frame(
    var    = rep(c("S_Preds","I_Preds","S_Albatross","I_Albatross"),
                 each = length(mig_arrival_times)),
    time   = rep(mig_arrival_times, times = 4),
    value  = rep(c(pred_S_arrival, pred_I_arrival, alb_S_arrival, alb_I_arrival),
                 each = length(mig_arrival_times)),
    method = rep("add", 4 * length(mig_arrival_times))
  )
  
  # Sort and align precision
  eventdata1      <- eventdata1[order(eventdata1$time), ]
  eventdata1$time <- round(eventdata1$time, 3)
  
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
#                                                 Generate Graphs:
#=======================================================================================================================#
graphMasterInf          <- function(dfMasterInfGentoo, speciesname){
  df_all <- dfMasterInfGentoo %>%
    rownames_to_column(var = "run") %>%        # bring rownames in as a column
    pivot_longer(
      cols        = -run,                      # all except 'run'
      names_to    = "time",                    
      values_to   = "value"
    ) %>%
    mutate(time = as.integer(time))            # ensure 'time' is numeric
  df_avg <- df_all %>%
    group_by(time) %>%
    summarise(
      avg_value = mean(value, na.rm = TRUE),
      .groups   = "drop"
    )
  p <- ggplot(df_all, aes(x = time, y = value, group = run)) +
    geom_line(color = "gray60", linetype = "dotted") +
    geom_line(
      data        = df_avg,
      aes(x = time, y = avg_value),
      inherit.aes = FALSE,
      color       = "red",
      linewidth   = 1.2,
      linetype    = "solid"
    ) +
    labs(
      title = paste0(speciesname, " Dead (from HPAI) Trajectories with Average"),
      x     = "Time (days)",
      y     = "Number Infected"
    ) +
    theme_minimal()
  return(p)
}

#=======================================================================================================================#
#                                                Model Functions:
#=======================================================================================================================#
seasonality_transmission_calc <- function(time, timevar){
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
  x03  = 357.5  # Midpoint of the third sigmoid
  a3   = 1.05    # Upper bound of the third sigmoid (previously `top`)
  b3   = 0.98      # Lower bound of the third sigmoid (previously `bottom`)
  
  ifelse(time < timevar$season_start, 
         b1 + (a1 - b1) / (1 + exp(-k1 * (time - x01))),
         ifelse(time <= timevar$season_end, 
                a2 - (a2 - b2) * (1 / (1 + exp(-k2 * (time - x02)))),
                b3 + (a3 - b3) / (1 + exp(-k3 * (time - x03)))
         )
  )
}                        # This function calculates the effect seasonality has on transmission rate
HPAI_dyn                      <- function(time, state, parameters){
  # Stops analysis if time is not a numeric variable
  if(!is.numeric(time)){
    stop("Time must be a numeric variable")
  }
  
  with(as.list(c(state,parameters)),{
    # If seasonality is activated, use the seasonal infectivity modifier 
    if(seasons_calc == TRUE){
      seasonality_Transmission = seasonality_transmission_calc(time, timevar)
    } else{
      seasonality_Transmission = rep(1, length(time))
    }
    
    #===============#
    #    Gentoo     #
    #===============#
    dS_Gentoo = 
      ifelse(time >= timevar$gentoo_breed_start & time <= timevar$gentoo_breed_finish, 2 * muH_Gentoo * (S_Gentoo + I_Gentoo + R_Gentoo), 0) - #Gentoo birth rate
      
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
      ifelse(time >= timevar$pred_breed_start & time <= timevar$pred_breed_finish, 2.5 * muH_Preds * (S_Preds + I_Preds + R_Preds), 0) -  # Predator birth rate
      
      seasonality_Transmission * gentpred_trans * I_Gentoo    * S_Preds -  # Transmission from Gentoo to Predators
      seasonality_Transmission * predpred_trans * I_Preds     * S_Preds -  # Transmission from Predators to other Predators
      seasonality_Transmission * albspred_trans * I_Albatross * S_Preds -  # Transmission from Albatross to Predators
      
      muH_Preds * S_Preds +                                                # Natural death rate
      sigma_preds * R_Preds
    dI_Preds = 
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
      ifelse(time >= timevar$albatross_breed_start & time <= timevar$albatross_breed_finish, 1.7 * muH_Albs * (S_Albatross + I_Albatross + R_Albatross), 0) -
      
      seasonality_Transmission * predalbs_trans * I_Preds     * S_Albatross - # Transmission from Predators to Albatross
      seasonality_Transmission * albsalbs_trans * I_Albatross * S_Albatross - # Transmission from Albatross to other Albatross
      
      muH_Albs * S_Albatross +                                                # Natural death rate
      sigma_albatross * R_Albatross                                           # Loss of Immunity
    
    dI_Albatross = 
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
run_Model                     <- function(timevar, popvar, val, Initial_values, HPAIdyn, parameters, time1, time2, pres, samplenum){
  # Create the migration in event
  eventdata1 <- create_eventdata1(timevar, popvar)
  
  # Validation of no errors in eventdata1
  if(val$check_eventdata1 == TRUE){
    print(eventdata1)
  }
  
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
  
  # Validation of no errors in eventdata2
  if(val$check_eventdata2 == TRUE){
    print(eventdata2)
  }
  
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
  
  
  # Connect pre- and post- migration outputs and then use human intelligible labels
  output <- rbind(out_pre, out_post)
  colnames(output) <- c("Time", 
                        "Susceptable_Gentoo", "Infected_Gentoo", "Recovered_Gentoo", "HPAI_Killed_Gentoo",
                        "Susceptable_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
                        "Susceptable_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross",
                        "Timer")
  return(output)
} # This function runs the model with a given set of conditions

#=======================================================================================================================#
#                                                   Load data:
#=======================================================================================================================#
fileLoc   = "C:\\Users\\tjm336\\OneDrive - Cornell University\\01 Lab Work\\01 Rotations\\01 Bento-Gamble\\04-TriSpecies-HPAI\\LHC Walking"
inputFile = "Foundation.csv"
rawdf     <- as.data.frame(read.csv(paste0(fileLoc, "\\", inputFile)))

pheatmap(
  rawdf,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "column",
  treeheight_row = 0,
  main = "Heatmap of LHS Parameters",
  fontsize_row = 0.1,
  fontsize_col = 10,
  color = colorRampPalette(c("blue","white","red"))(100)
)


overwinterdf  <- rawdf[rawdf$overwintering == 1, ]
overwinterlen <- nrow(overwinterdf)

n_tp <- as.numeric(timevar$horizon) + 2
dfMasterInfGentoo    <- data.frame(matrix(NA, nrow = as.numeric(overwinterlen) + 1, ncol = n_tp))
dfMasterInfPredators <- data.frame(matrix(NA, nrow = as.numeric(overwinterlen) + 1, ncol = n_tp))
dfMasterInfAlbatross <- data.frame(matrix(NA, nrow = as.numeric(overwinterlen) + 1, ncol = n_tp))
for(i in 1:overwinterlen){
  samplenum       <- overwinterdf$X[i]
  muI_Gentoo      <- overwinterdf$muI_gentoo[i]
  gamma_gentoo    <- overwinterdf$gamma_genetoo[i]
  gentgent_trans  <- overwinterdf$Beta_GentGent[i]
  gentpred_trans  <- overwinterdf$Beta_GentPred[i]
  
  prev_infection_preds <- overwinterdf$prev_infection_preds[i]
  muI_Preds       <- overwinterdf$muI_predators[i]
  gamma_preds     <- overwinterdf$gamma_predators[i]
  predgent_trans  <- overwinterdf$Beta_PredGent[i]
  predpred_trans  <- overwinterdf$Beta_PredPred[i]
  predalbs_trans  <- overwinterdf$Beta_PredAlbs[i]
  
  prev_infection_albs <- overwinterdf$prev_infection_albs[i]
  muI_Albs        <- overwinterdf$muI_albatross[i]
  gamma_albatross <- overwinterdf$gamma_albatross[i]
  albspred_trans  <- overwinterdf$Beta_AlbsPred[i]
  albsalbs_trans  <- overwinterdf$Beta_AlbsAlbs[i]
  
  parameters <- c(
    muH_Gentoo,
    muI_Gentoo,
    sigma_gentoo,
    gamma_gentoo,
    gentgent_trans,
    gentpred_trans,
    
    prev_infection_preds,
    muH_Preds,
    muI_Preds,
    sigma_preds,
    gamma_preds,
    predgent_trans,
    predpred_trans,
    predalbs_trans,
    
    prev_infection_albs,
    muH_Albs,
    muI_Albs,
    sigma_albatross,
    gamma_albatross,
    albspred_trans,
    albsalbs_trans
  )
  popvar$prev_infection_preds <- prev_infection_preds
  popvar$prev_infection_albs  <- prev_infection_albs
  
  output <- run_Model(timevar, popvar, val, Initial_values, HPAIdyn, parameters, time1, time2, pres, samplenum)
#  print(length(output$Infected_Gentoo))
  dfMasterInfGentoo[i,] <- output$HPAI_Killed_Gentoo
  rownames(dfMasterInfGentoo)[i] <- paste0("S", samplenum)
  dfMasterInfPredators[i,] <- output$HPAI_Killed_Predators
  rownames(dfMasterInfPredators)[i] <- paste0("S", samplenum)
  dfMasterInfAlbatross[i,] <- output$HPAI_Killed_Albatross
  rownames(dfMasterInfAlbatross)[i] <- paste0("S", samplenum)
  
}

colnames(dfMasterInfGentoo) <- output$Time
dfMasterInfGentoo <- dfMasterInfGentoo[, !duplicated(names(dfMasterInfGentoo))]
colnames(dfMasterInfPredators) <- output$Time
dfMasterInfPredators <- dfMasterInfPredators[, !duplicated(names(dfMasterInfPredators))]
colnames(dfMasterInfAlbatross) <- output$Time
dfMasterInfAlbatross <- dfMasterInfAlbatross[, !duplicated(names(dfMasterInfAlbatross))]


p1 <- graphMasterInf(dfMasterInfGentoo, "Gentoo")
print(p1)
p2 <- graphMasterInf(dfMasterInfPredators, "Predators")
print(p2)
p3 <- graphMasterInf(dfMasterInfAlbatross, "Albatross")
print(p3)

#print(dfMasterInfGentoo)

