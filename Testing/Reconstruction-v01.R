#=======================================================================================================================#
#                                               Control Variables:
#=======================================================================================================================#
# Validation Variable
check_initVal<- FALSE
check_season <- TRUE
check_PFG    <- FALSE
check_PFA    <- FALSE
check_nPFG   <- FALSE
check_nPFA   <- FALSE
print_head   <- FALSE
print_tail   <- TRUE

# Experimental Variables
seasons_calc <- FALSE
migrate_out  <- TRUE

#=======================================================================================================================#
#                                                   Libraries:
#=======================================================================================================================#
library(ggplot2)
library(deSolve)

#=======================================================================================================================#
#                                                   Populations:
#=======================================================================================================================#
colony_num            <- 20                                         # Number of colonies on the island

total_gentoo          <- 200000 / colony_num                        # Total Gentoo population [Source: 1]
ResPer_gentoo         <- 1                                          # Percent of Gentoo that are resident to the island [Source: 1]
Res_gentoo            <- total_gentoo * ResPer_gentoo               # Calculate the number of resident Gentoo
gentoo_migrate_in     <- total_gentoo * (1 - ResPer_gentoo)         # Calculate the number of migratory Gentoo

total_preds           <- 74000 / colony_num                         # Total predator population [Source: 1]
ResPer_preds          <- 0.1                                        # Percent of predator that are resident to the island [Source: 1]
Res_preds             <- total_preds * ResPer_preds                 # Calculate the number of resident predators 
preds_migrate_in      <- total_preds * (1 - ResPer_preds)           # Calculate the number of migratory predators
prev_infection_preds  <- 0.00                                       # Rate at which arriving predators are infected


total_albatross       <- 951000 / colony_num                        # Total Albatross population [Source: 1]
ResPer_albatross      <- 0                                          # Percent of Albatross that are resident to the island [Source: 1]
Res_albatross         <- total_albatross * ResPer_albatross         # Resident Albatross 
albatross_migrate_in  <- total_albatross * (1 - ResPer_albatross)   # 100% of albatross are migratory
prev_infection_albs   <- 0.01                                       # Rate at which arriving Albatross are infected

Initial_values <- c(
  S_Preds      = Res_preds,       # Susceptible
  I_Preds      = 0,               # Infected                   
  R_Preds      = 0,               # Recovered 
  D_Preds      = 0,               # Dead
  
  S_Gentoo = Res_gentoo ,         # Susceptible
  I_Gentoo = 0,                   # Infected                   
  R_Gentoo = 0,                   # Recovered
  D_Gentoo = 0,                   # Dead
  
  S_Albatross  = Res_albatross,   # Susceptible
  I_Albatross  = 0,               # Infected                   
  R_Albatross  = 0,               # Recovered 
  D_Albatross  = 0,               # Dead
  
  timer        = 0)

if(check_initVal == TRUE){
  print(Initial_values)
} 

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
gentoo_inf_dur        <- 7                                          # Duration of infection in Gentoo penguins [Source: 5]
gamma_gentoo          <- 1 / gentoo_inf_dur                         # Calculate the recovery rate for Gentoo penguins
gentoo_inf_mort       <- 0.5                                        # Estimate without a basis
muI_Gentoo            <- gentoo_inf_mort / gentoo_inf_dur           # Calculation of the infected Gentoo death rate
gentoo_imm_dur        <- 182.5                                      # Estimate without a basis
sigma_gentoo          <- 1 / gentoo_imm_dur                         # Rate of Gentoo immunity loss

preds_inf_dur         <- 7                                          # Duration of infection in Predator species [Source: 5]
gamma_preds           <- 1 / preds_inf_dur                          # Calculate the recovery rate for Predator species 
preds_inf_mort        <- 0                                          # Estimate without a basis
muI_Preds             <- muH_Preds                                  # Calculation of the infected predator species death rate
Preds_imm_dur         <- 182.5                                      # Estimate without a basis
sigma_preds           <- 1 / Preds_imm_dur                          # Rate of predator immunity loss


Albs_inf_dur          <- 7                                          # Duration of infection in Albatross penguins [Source: 5]
gamma_albatross       <- 1 / Albs_inf_dur                           # Calculate the recovery rate for Albatross 
Albs_inf_mort         <- 0.5                                        # Estimate without a basis
muI_Albs              <- Albs_inf_mort / Albs_inf_dur               # Calculation of the infected Gentoo death rate
albatross_imm_dur     <- 182.5                                      # Estimate without a basis
sigma_albatross       <- 1 / albatross_imm_dur                      # Rate of Albatross immunity loss

#=======================================================================================================================#
#                                                Temporal Values:
#=======================================================================================================================#
horizon               <- 365                                        # Number of days to be modeled
time                  <- seq(from = 0, to = horizon, by = 1)        # Set the time sequence
season_start          <- 100                                        # First day of the season for seasonality tranmissibility calculation
season_end            <- 349                                        # Last day of the season for seasonality tranmissibility calculation
t_arrival             <- 62                                         # Day at which migratory birds begin arriving, unclear where the source is
mig_day_arrival       <- 31                                         # The duration of arrival time of migratory birds, no clear source
t_departure           <- 300                                        # Day at which migratory birds begin departing, unclear where the source is
mig_day_departure     <- 10                                         # The duration of departure time of migratory birds, no clear source
gentoo_breed_start    <- 31                                         # Gentoo breeding cycle start [Source: Unclear]
gentoo_breed_finish   <- 212                                        # Gentoo breeding cycle end   [Source: Unclear]
pred_breed_start      <- 93                                         # Predator breeding cycle start   [Source: Unclear]
pred_breed_finish     <- 225                                        # Predator breeding cycle finish  [Source: Unclear]
albatross_breed_start <- 31
albatross_breed_finish<- 245

#=======================================================================================================================#
#                                                    Betas:
#=======================================================================================================================#
# Transmission Rates From Gentoo to (other) Genntoos, Predators, or Albatrosses 
gentgent_trans <- 0.000003
gentpred_trans <- 0.000003
gentalbs_trans <- 0.000003

# Transmission Rates From Predators to Gentoos, (other) Predators, or Albatrosses 
predgent_trans <- 0.000003
predpred_trans <- 0.000003
predalbs_trans <- 0.000003

# Transmission Rates From Albatrosses to Gentoos, Predators, or (other) Albatrosses 
albsgent_trans <- 0.000003
albspred_trans <- 0.000003
albsalbs_trans <- 0.000003

#=======================================================================================================================#
#                                                Establish Parameters:
#=======================================================================================================================#
parameters <- c(
  muH_Gentoo,
  muI_Gentoo,
  sigma_gentoo,
  gamma_gentoo,
  gentgent_trans,
  gentpred_trans,
  gentalbs_trans,
  
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
  albsgent_trans,
  albspred_trans,
  albsalbs_trans
)
#print(parameters)

#=======================================================================================================================#
#                                                Model Functions:
#=======================================================================================================================#
seasonality_transmission_calc <- function(time){
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
  
  ifelse(time < season_start, 
         b1 + (a1 - b1) / (1 + exp(-k1 * (time - x01))),
         ifelse(time <= season_end, 
                a2 - (a2 - b2) * (1 / (1 + exp(-k2 * (time - x02)))),
                b3 + (a3 - b3) / (1 + exp(-k3 * (time - x03)))
         )
  )
}                        # This function calculates the effect seasonality has on transmission rate
PiecewiseFlux                 <- function(time, t1, t2, b){
  if(time >= t1 && time <= t2){
    return(b)
  } else{
    return(0)
  }
}             # This function models the arrival of Predators onto the island
HPAI_dyn                      <- function(time, state, parameters){
  # Stops analysis if time is not a numeric variable
  if(!is.numeric(time)){
    stop("Time must be a numeric variable")
  }
  
  with(as.list(c(state,parameters)),{
    # If seasonality is activated, use the seasonal infectivity modifier 
    if(seasons_calc == TRUE){
      seasonality_Transmission = seasonality_transmission_calc(time)
    } else{
      seasonality_Transmission = rep(1, length(time))
    }
    
    eps <- 1e-8
    #===============#
    #    Gentoo     #
    #===============#
    dS_Gentoo = 
      ifelse(time >= gentoo_breed_start & time <= gentoo_breed_finish, 2 * muH_Gentoo * (S_Gentoo + I_Gentoo + R_Gentoo), 0) - #Gentoo birth rate
      
      seasonality_Transmission * gentgent_trans * I_Gentoo    * S_Gentoo - # Transmission from Gentoo to other Gentoo
      seasonality_Transmission * predgent_trans * I_Preds     * S_Gentoo - # Transmission from Predators to Gentoo to Gentoo

       muH_Gentoo * S_Gentoo  +                                            # Natural death rate
       sigma_gentoo * R_Gentoo                                             # Loss of Gentoo immunity

    dI_Gentoo = 
      seasonality_Transmission * gentgent_trans * I_Gentoo    * S_Gentoo + # Transmission from Gentoo to other Gentoo
      seasonality_Transmission * predgent_trans * I_Preds     * S_Gentoo - # Transmission from Predators to Gentoo

      muI_Gentoo * I_Gentoo -                                              # Infected death rate
      gamma_gentoo * I_Gentoo                                              # Recovery rate
    
    dR_Gentoo = 
      gamma_gentoo * I_Gentoo -                                            # Recovery rate
      muH_Gentoo * R_Gentoo -                                              # Natural death rate
      sigma_gentoo * R_Gentoo                                              # Loss of immunity rate
    
    dD_Gentoo = muI_Gentoo * I_Gentoo
    
    #===============#
    #   Predators   #
    #===============#
    dS_Preds = 
      (1 - prev_infection_preds) * PiecewiseFlux(time, t_arrival - eps, t_arrival + mig_day_arrival + eps, preds_migrate_in/mig_day_arrival) + # Arrival of susceptible birds via migration
      ifelse(time >= pred_breed_start & time <= pred_breed_finish, 2.5 * muH_Preds * (S_Preds + I_Preds + R_Preds), 0) -  # Predator birth rate
      
      seasonality_Transmission * gentpred_trans * I_Gentoo    * S_Preds -  # Transmission from Gentoo to Predators
      seasonality_Transmission * predpred_trans * I_Preds     * S_Preds -  # Transmission from Predators to other Predators
      seasonality_Transmission * albspred_trans * I_Albatross * S_Preds -  # Transmission from Albatross to Predators
      
      muH_Preds * S_Preds +                                                # Natural death rate
      sigma_preds * R_Preds
    dI_Preds = 
      prev_infection_preds * PiecewiseFlux(time, t_arrival - eps, t_arrival + mig_day_arrival + eps, preds_migrate_in/mig_day_arrival) + # Arrival of infected birds via migration
      
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
      (1 - prev_infection_albs) * PiecewiseFlux(time, t_arrival - eps, t_arrival + mig_day_arrival + eps, preds_migrate_in/mig_day_arrival) + # Arrival of susceptible birds via migration
      ifelse(time >= albatross_breed_start & time <= albatross_breed_finish, 1.7 * muH_Albs * (S_Albatross + I_Albatross + R_Albatross), 0) -
      
      seasonality_Transmission * predalbs_trans * I_Preds     * S_Albatross - # Transmission from Predators to Albatross
      seasonality_Transmission * albsalbs_trans * I_Albatross * S_Albatross - # Transmission from Albatross to other Albatross
      
      muH_Albs * S_Albatross +                                                # Natural death rate
      sigma_albatross * R_Albatross                                           # Loss of Immunity
    
    dI_Albatross = 
      prev_infection_albs * PiecewiseFlux(time, t_arrival - eps, t_arrival + mig_day_arrival + eps, preds_migrate_in/mig_day_arrival) + # Arrival of susceptible birds via migration
      
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
}

#=======================================================================================================================#
#                                              Validation Functions:
#=======================================================================================================================#
validate_Seasonality          <- function(func, time){
  # Calculate values
  values <- vapply(time, func, FUN.VALUE = numeric(1))
  #print(values)
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
    )+
    theme_minimal()
  
  # Return the values calculated
  list(data = df, plot = p)
}                  # This function validates the "seasonality_transmission_calc"
validate_PFlux                <- function(func, time, t1, t2, b, spec){
  # Calculate values
  values <- vapply(time, func, FUN.VALUE = numeric(1), t1 = t1, t2 = t2, b = b)
  #print(values)
  df <- data.frame(
    time = time,
    value = values
  )
  
  # Generate Plot for Values
  p <- ggplot(df, aes(x = time, y = value)) +
    geom_line(color = "steelblue", size = 1) +
    labs(
      title = paste0(spec, " Arrival Over Time"),
      x     = "Time Interval (Days)",
      y     = "Computed Value"
    ) +
    theme_minimal()
  
  # Return List
  list(data = df, plot = p)
} # This function validates the PiecewiseFlux function for migration in
validate_nPFlux               <- function(func, time, t1, t2, b, spec){
  # Calculate Values
  values <- vapply(time, func, FUN.VALUE = numeric(1), t1 = t1, t2 = t2, b = b)
  #print(values)
  df <- data.frame(
    time = time,
    value = values
  )
  
  # Generate Plot for Values
  p <- ggplot(df, aes(x = time, y = value)) +
    geom_line(color = "steelblue", size = 1) +
    labs(
      title = paste0(spec, " Departure Over Time"),
      x     = "Time Interval (Days)",
      y     = "Computed Value"
    ) +
    theme_minimal()
  
  # Return List
  list(data = df, plot = p)
} # This function validates the PiecewiseFlux function for mirgration out

#=======================================================================================================================#
#                                                      Execution:
#=======================================================================================================================#
# Validation checklist
if(check_season == TRUE){
  Seas_Val <- validate_Seasonality(seasonality_transmission_calc, time)
  Seas_Val$plot
}            # validate_Seasonality activation
if(check_PFG    == TRUE){
  pfg_Val <- validate_PFlux(PiecewiseFlux, time, t_arrival, t_arrival + mig_day_arrival, preds_migrate_in/mig_day_arrival, "Predators")
  pfg_Val$plot
}            # validate_PFLux activation for Predators
if(check_PFA    == TRUE){
  pfa_Val <- validate_PFlux(PiecewiseFlux, time, t_arrival, t_arrival + mig_day_arrival, albatross_migrate_in/mig_day_arrival, "Albatross")
  pfa_Val$plot
}            # validate_PFLux activation for Albatross
if(check_nPFG   == TRUE){
  preds_migrating_out = total_preds * (1 - ResPer_preds)
  npfg_Val <- validate_nPFlux(PiecewiseFlux, time, t_departure, t_departure + mig_day_departure,  preds_migrating_out/mig_day_departure, "Predators")
  npfg_Val$plot
}            # validate_nPFLux activation for Predators
if(check_nPFA   == TRUE){
  albatross_migrating_out = total_albatross * (1 - ResPer_albatross)
  npfa_Val <- validate_nPFlux(PiecewiseFlux, time, t_departure, t_departure + mig_day_departure, albatross_migrating_out/mig_day_departure, "Albatross")
  npfa_Val$plot
}            # validate_nPFLux activation for Albatross

# Model Running
output <- as.data.frame(lsoda(
  y     = Initial_values,
  func  = HPAI_dyn,
  parms = parameters,
  time  = time
))

colnames(output) <- c("Time", 
                      "Susceptable_Gentoo", "Infected_Gentoo", "Recovered_Gentoo", "HPAI_Killed_Gentoo",
                      "Susceptable_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
                      "Susceptable_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross",
                      "Timer")
if(print_head == TRUE){
  head(output, 100)
}

if(print_tail == TRUE){
  tail(output, 10)
}

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